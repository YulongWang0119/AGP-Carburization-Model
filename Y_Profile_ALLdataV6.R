# --- 0. 載入必要的套件 ---
library(ggplot2)
library(dplyr)   # 用於 bind_rows，處理 data.frame 合併
library(purrr)   
library(readr)  
library(readr)
# ===================================================================
#           PART 1: 核心模擬函式
# ===================================================================

#' 運行一次多階段滲碳模擬
#'
#' @param recipe_stages 一個 list，定義了每個製程階段的參數。
#' @param L 鋼材厚度 (m)。
#' @param C_o 初始碳濃度 (wt%)。
#' @param Nx 空間網格點數量。
#' @return 一個包含 'final_profile' 和 'plot_data' 的 list。


run_simulation <- function(recipe_stages, L=0.008, C_o=0.2, Nx=101) {
  # 阿瑞尼斯方程式參數
  D0 <- 2.3e-5      # m^2/s
  Qd <- 148000      # J/mol
  R_const <- 8.314  # J/(mol·K)
  
  # --- 1a. 物理參數 ---
  calculate_D <- function(temp_C) { D0 * exp(-Qd / (R_const * (temp_C + 273.15))) }
  
  # --- 1b. 數值計算參數 ---
  total_time_s <- sum(sapply(recipe_stages, function(s) s$duration)) * 3600
  max_temp <- max(sapply(recipe_stages, function(s) s$temp_C))
  D_max <- calculate_D(max_temp)
  x <- seq(0, L, length.out = Nx); dx <- x[2] - x[1]
  dt <- 0.9 * (dx^2) / (2 * D_max); Nt <- as.integer(total_time_s / dt)
  
  # --- 1c. 初始化 ---
  C <- rep(C_o, Nx)
  plot_df <- data.frame(
    depth_mm = x * 1000, concentration = C, time_hr = 0, stage = "Initial"
  )
  stage_end_times_s <- cumsum(sapply(recipe_stages, function(s) s$duration * 3600))
  
  # --- 1d. 時間推進主迴圈 ---
  current_stage_index <- 1
  
  for (n in 1:Nt) {
    # 核心 FDM 計算
    current_stage <- recipe_stages[[current_stage_index]]
    current_D <- calculate_D(current_stage$temp_C)
    current_Cs <- current_stage$surface_C
    C_new <- C; lam <- current_D * dt / dx^2
    C_new[2:(Nx - 1)] <- C[2:(Nx - 1)] + lam * (C[3:Nx] - 2 * C[2:(Nx - 1)] + C[1:(Nx - 2)])
    C_new[1] <- current_Cs; C_new[Nx] <- C_new[Nx - 1]
    C <- C_new
    
    # 在計算完成後，檢查是否到達階段末尾
    if (current_stage_index < length(recipe_stages)) { # 只檢查到倒數第二個階段
      if ((n * dt) >= stage_end_times_s[current_stage_index]) {
        stage_info <- recipe_stages[[current_stage_index]] #取出當前這個即將結束的階段的詳細說明名字
        temp_df <- data.frame( 
          depth_mm = x * 1000, concentration = C,
          time_hr = round(stage_end_times_s[current_stage_index] / 3600, 1),
          stage = sprintf("End of %s", stage_info$name)
        )
        plot_df <- rbind(plot_df, temp_df)
        
        current_stage_index <- current_stage_index + 1 # 切換階段
      }
    }
  }
  
  # --- 1e. 手動記錄最終階段的狀態 ---
  # 在迴圈結束後，C 向量儲存的是最終時刻的濃度分佈
  final_stage_info <- recipe_stages[[length(recipe_stages)]]
  final_temp_df <- data.frame(
    depth_mm = x * 1000, concentration = C,
    time_hr = round(total_time_s / 3600, 1),
    stage = sprintf("End of %s", final_stage_info$name)
  )
  
  plot_df <- rbind(plot_df, final_temp_df)
  
  
  # --- 1e. 準備返回值 ---
  final_profile <- data.frame(depth_mm = x * 1000, concentration_pct = C)
  
  return(list(final_profile = final_profile, plot_data = plot_df))
}


# ===================================================================
#           PART 2: 資料庫儲存函式
# ===================================================================

# 建立資料庫基本結構
setup_database <- function(base_path) {
  dir.create(base_path, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(base_path, "profiles"), showWarnings = FALSE)
  dir.create(file.path(base_path, "plot_data"), showWarnings = FALSE)
  
  recipes_path <- file.path(base_path, "recipes.csv")
  if (!file.exists(recipes_path)) {
    empty_df <- data.frame(Run_ID=character(), Stage_Type=character(),
                           Profile_File_Path=character(), Plot_Data_File_Path=character())
    write.csv(empty_df, recipes_path, row.names = FALSE)
  }
}

# 新增 stage_type 參數，寫入製程類型標籤
save_run_to_db <- function(base_path, run_id, recipe, simulation_results, stage_type) {
  
  final_profile <- simulation_results$final_profile
  plot_data <- simulation_results$plot_data
  
  # 1. 準備 Recipe 主表
  recipe_list <- list(Run_ID = run_id, Stage_Type = stage_type) # 直接寫入標籤
  
  for (i in 1:length(recipe)) {
    stage <- recipe[[i]]
    recipe_list[[paste0("Stage", i, "_Temp_C")]] <- stage$temp_C
    recipe_list[[paste0("Stage", i, "_Time_hr")]] <- stage$duration
    recipe_list[[paste0("Stage", i, "_Cs_pct")]] <- stage$surface_C
  }
  
  profile_filename <- paste0(run_id, "_profile.csv")
  plot_data_filename <- paste0(run_id, "_plot_data.csv")
  recipe_list$Profile_File_Path <- file.path("profiles", profile_filename)
  recipe_list$Plot_Data_File_Path <- file.path("plot_data", plot_data_filename)
  
  recipes_path <- file.path(base_path, "recipes.csv")
  existing_recipes <- read_csv(recipes_path, col_types = cols(.default = "c")) %>% mutate_all(as.character)
  new_recipe_df <- as.data.frame(lapply(recipe_list, as.character))
  updated_recipes <- bind_rows(existing_recipes, new_recipe_df)
  write_csv(updated_recipes, recipes_path, na = "")
  
  # 2. 儲存 Raw Profile
  write_csv(final_profile, file.path(base_path, "profiles", profile_filename))
  
  # 3. 儲存 Plotting Data
  if (!is.null(plot_data)) {
    write_csv(plot_data, file.path(base_path, "plot_data", plot_data_filename))
  }
}

log_iteration <- function(base_path, iter, af_value, selected_run_id) {
  log_df <- data.frame(
    Iteration = iter, Timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    Acquisition_Function = "EI", Max_AF_Value = af_value, Selected_Run_ID = selected_run_id
  )
  write.table(log_df, file = file.path(base_path, "optimization_log.csv"), append = TRUE,
              sep = ",", row.names = FALSE, col.names = FALSE)
}

# ===================================================================
#           PART 3 : 定義連續參數空間與生成多階段隨機配方
# ===================================================================

# --- 3a. 定義製程參數的連續範圍 ---
param_ranges_boost <- list( # 單數階 (Boost)
  temp_C    = c(850, 1000),
  surface_C = c(1.2, 2.4),
  duration  = c(8, 20)      # 時間較長
)
param_ranges_diffuse <- list( # 雙數階 (Diffuse)
  temp_C    = c(850, 1000),
  surface_C = c(0.4, 1.3),
  duration  = c(4, 10)       # 時間較短
)

# --- 3b. 生成隨機的單階段配方 ---
num_single_stage <- 1000
set.seed(2027666) 

single_stage_recipes_v4 <- data.frame(
  temp_C    = round(runif(num_single_stage, min = param_ranges_boost$temp_C[1], max = param_ranges_boost$temp_C[2])),
  surface_C = round(runif(num_single_stage, min = param_ranges_boost$surface_C[1], max = param_ranges_boost$surface_C[2]), 2),
  duration  = round(runif(num_single_stage, min = param_ranges_boost$duration[1], max = param_ranges_boost$duration[2]), 1)
)
cat(sprintf("已生成 %d 組隨機的【單階段】製程配方。\n", nrow(single_stage_recipes_v4)))

# --- 3c. 生成滿足物理約束的隨機雙階段配方 ---
num_dual_stage <- 4000
dual_stage_recipes_v4 <- list()

cat("正在生成【雙階段】配方...\n")
while(length(dual_stage_recipes_v4) < num_dual_stage) {
  s1_params <- list(
    temp_C    = round(runif(1, min = param_ranges_boost$temp_C[1], max = param_ranges_boost$temp_C[2])),
    surface_C = round(runif(1, min = param_ranges_boost$surface_C[1], max = param_ranges_boost$surface_C[2]), 2),
    duration  = round(runif(1, min = param_ranges_boost$duration[1], max = param_ranges_boost$duration[2]), 1)
  )
  s2_params <- list(
    temp_C    = round(runif(1, min = param_ranges_diffuse$temp_C[1], max = param_ranges_diffuse$temp_C[2])),
    surface_C = round(runif(1, min = param_ranges_diffuse$surface_C[1], max = param_ranges_diffuse$surface_C[2]), 2),
    duration  = round(runif(1, min = param_ranges_diffuse$duration[1], max = param_ranges_diffuse$duration[2]), 1)
  )
  
  if (s1_params$surface_C > s2_params$surface_C && s1_params$duration > s2_params$duration) {
    dual_stage_recipes_v4[[length(dual_stage_recipes_v4) + 1]] <- list(Stage1 = s1_params, Stage2 = s2_params)
  }
}
cat(sprintf("已成功生成 %d 組【雙階段】製程配方。\n", length(dual_stage_recipes_v4)))

# --- 3d. 生成三階段配方 (B-D-B 和 B-D-D 嚴格各半) ---
num_triple_stage_total <- 7000
num_per_type <- floor(num_triple_stage_total / 2) # 每種類型生成的數量
triple_stage_recipes_v4 <- list()

cat("正在生成【三階段】配方 (B-D-B 和 B-D-D 各半)...\n")

# 生成三階段配方
generate_one_triple_recipe <- function(s3_type) {
  repeat {
    s1_params <- list(
      temp_C = round(runif(1, min = param_ranges_boost$temp_C[1], max = param_ranges_boost$temp_C[2])),
      surface_C = round(runif(1, min = param_ranges_boost$surface_C[1], max = param_ranges_boost$surface_C[2]), 2),
      duration = round(runif(1, min = param_ranges_boost$duration[1], max = param_ranges_boost$duration[2]), 1)
    )
    s2_params <- list(
      temp_C = round(runif(1, min = param_ranges_diffuse$temp_C[1], max = param_ranges_diffuse$temp_C[2])),
      surface_C = round(runif(1, min = param_ranges_diffuse$surface_C[1], max = param_ranges_diffuse$surface_C[2]), 2),
      duration = round(runif(1, min = param_ranges_diffuse$duration[1], max = param_ranges_diffuse$duration[2]), 1)
    )
    
    # 根據傳入的類型決定第三階段的範圍
    param_ranges_s3 <- if (s3_type == "Boost") param_ranges_boost else param_ranges_diffuse
    
    s3_params <- list(
      temp_C = round(runif(1, min = param_ranges_s3$temp_C[1], max = param_ranges_s3$temp_C[2])),
      surface_C = round(runif(1, min = param_ranges_s3$surface_C[1], max = param_ranges_s3$surface_C[2]), 2),
      duration = round(runif(1, min = param_ranges_s3$duration[1], max = param_ranges_s3$duration[2]), 1)
    )
    
    # 檢查核心約束
    if (s1_params$surface_C > s2_params$surface_C && s1_params$duration > s2_params$duration) {
      return(list(Stage1 = s1_params, Stage2 = s2_params, Stage3 = s3_params))
    }
  }
}

# 生成第一種類型: B-D-B
for(i in 1:num_per_type) {
  triple_stage_recipes_v4[[length(triple_stage_recipes_v4) + 1]] <- generate_one_triple_recipe("Boost")
}
cat(sprintf("已生成 %d / %d 筆 B-D-B 配方...\n", length(triple_stage_recipes_v4), num_triple_stage_total))

# 生成第二種類型: B-D-D
while(length(triple_stage_recipes_v4) < num_triple_stage_total) { # 使用 while 補足剩餘數量
  triple_stage_recipes_v4[[length(triple_stage_recipes_v4) + 1]] <- generate_one_triple_recipe("Diffuse")
}

cat(sprintf("已成功生成 %d 組【三階段】製程配方。\n", length(triple_stage_recipes_v4)))

# ===================================================================
#           PART 4: 執行模擬並建立 V6 資料庫
# ===================================================================

# --- 4a. 定義 V6 資料庫的主路徑 ---
db_main_path_v6 <- "C:/Users/USER/Desktop/PYCProfessor/CarburizationDB_ProfileV6"
cat(sprintf("所有 V6 權威資料庫將存放於: %s\n", db_main_path_v6))

# 初始化 V6 資料庫結構
setup_database(db_main_path_v6)

# --- 4b . 執行並儲存單階段數據 ---
cat(sprintf("\n--- 開始生成 %d 筆單階段數據 (V6) ---\n", nrow(single_stage_recipes_v4)))
for (i in 1:nrow(single_stage_recipes_v4)) {
  run_id <- sprintf("S1_R%04d", i)
  params <- single_stage_recipes_v4[i, ]
  recipe <- list(list(name = "Stage1", temp_C = params$temp_C, surface_C = params$surface_C, duration = params$duration))
  
  if (i %% 50 == 0 || i == nrow(single_stage_recipes_v4)) {
    cat(sprintf("\rRunning Single-Stage %s (%d/%d)...", run_id, i, nrow(single_stage_recipes_v4)))
    flush.console()
  }
  
  sim_results <- run_simulation(recipe_stages = recipe)
  # 直接傳入 stage_type 標籤
  save_run_to_db(base_path = db_main_path_v6, run_id = run_id, recipe = recipe, 
                 simulation_results = sim_results, stage_type = "Single_Stage")
}
cat("\n--- 單階段 V6 數據生成完畢！ ---\n")


# --- 4c. 執行並儲存雙階段數據 ---
cat(sprintf("\n--- 開始生成 %d 筆雙階段數據 (V6) ---\n", length(dual_stage_recipes_v4)))
for (i in 1:length(dual_stage_recipes_v4)) {
  run_id <- sprintf("S2_R%04d", i)
  s1_params <- dual_stage_recipes_v4[[i]]$Stage1
  s2_params <- dual_stage_recipes_v4[[i]]$Stage2
  recipe <- list(
    list(name = "Stage1_Boost", temp_C = s1_params$temp_C, surface_C = s1_params$surface_C, duration = s1_params$duration),
    list(name = "Stage2_Diffuse", temp_C = s2_params$temp_C, surface_C = s2_params$surface_C, duration = s2_params$duration)
  )
  
  if (i %% 50 == 0 || i == length(dual_stage_recipes_v4)) {
    cat(sprintf("\rRunning Dual-Stage %s (%d/%d)...", run_id, i, length(dual_stage_recipes_v4)))
    flush.console()
  }
  
  sim_results <- run_simulation(recipe_stages = recipe)
  # 傳入 stage_type 標籤
  save_run_to_db(base_path = db_main_path_v6, run_id = run_id, recipe = recipe, 
                 simulation_results = sim_results, stage_type = "Dual_BD")
}
cat("\n--- 雙階段 V6 數據生成完畢！ ---\n")


# --- 4d. 執行並儲存三階段數據 (B-D-B 和 B-D-D) ---
cat(sprintf("\n--- 開始生成 %d 筆三階段數據 (V6) ---\n", length(triple_stage_recipes_v4)))
num_bdb_type <- num_per_type 

for (i in 1:length(triple_stage_recipes_v4)) {
  run_id <- sprintf("S3_R%04d", i)
  
  # 根據生成時的順序，賦予正確的標籤
  type_label <- if (i <= num_bdb_type) "Triple_BDB" else "Triple_BDD"
  
  s1_params <- triple_stage_recipes_v4[[i]]$Stage1
  s2_params <- triple_stage_recipes_v4[[i]]$Stage2
  s3_params <- triple_stage_recipes_v4[[i]]$Stage3
  recipe <- list(
    list(name = "Stage1_Boost", temp_C = s1_params$temp_C, surface_C = s1_params$surface_C, duration = s1_params$duration),
    list(name = "Stage2_Diffuse", temp_C = s2_params$temp_C, surface_C = s2_params$surface_C, duration = s2_params$duration),
    list(name = "Stage3", temp_C = s3_params$temp_C, surface_C = s3_params$surface_C, duration = s3_params$duration)
  )
  
  if (i %% 50 == 0 || i == length(triple_stage_recipes_v4)) {
    cat(sprintf("\rRunning Triple-Stage %s (%s) (%d/%d)...", run_id, type_label, i, length(triple_stage_recipes_v4)))
    flush.console()
  }
  
  sim_results <- run_simulation(recipe_stages = recipe)
  # 傳入 stage_type 標籤
  save_run_to_db(base_path = db_main_path_v6, run_id = run_id, recipe = recipe, 
                 simulation_results = sim_results, stage_type = type_label)
}
cat("\n--- 三階段 V6 數據生成完畢！ ---\n")

cat("\n\n=============== 所有 V6 權威資料庫已生成完畢！ ===============\n")
# =========================================================================
#       V6 資料庫分類曲線視覺化腳本
# =========================================================================
# 目的：讀取 V6 統一資料庫，並根據 recipes.csv 中的 "Stage_Type" 欄位，
#       為每一種類型的製程曲線，各自生成一張視覺圖。
db_main_path_v6 <- "C:/Users/USER/Desktop/PYCProfessor/CarburizationDB_ProfileV6"

# 建立一個新的資料夾來儲存本次視覺化的所有圖檔
output_dir <- file.path(db_main_path_v6, "Visualization_V6_Output")
dir.create(output_dir, showWarnings = FALSE)
cat(sprintf("所有 V6 視覺化圖檔將儲存於: %s\n", output_dir))


# --- PART 2: 載入 V6 數據 ---
cat("--- 正在讀取 V6 資料庫...\n")

# 2a. 直接讀取 V6 資料庫根目錄下的 recipes.csv
recipes_path_v6 <- file.path(db_main_path_v6, "recipes.csv")
if (!file.exists(recipes_path_v6)) {
  stop("在指定的 V6 路徑下找不到 recipes.csv！請檢查路徑是否正確。")
}
recipes_all_v6 <- read_csv(recipes_path_v6, show_col_types = FALSE)

# 2b. 根據 recipes.csv 中的路徑，讀取所有 profile 檔案
all_profiles_list <- list()
for (i in 1:nrow(recipes_all_v6)) {
  current_recipe <- recipes_all_v6[i, ]
  profile_full_path <- file.path(db_main_path_v6, current_recipe$Profile_File_Path)
  
  if(file.exists(profile_full_path)) {
    profile_data <- read_csv(profile_full_path, show_col_types = FALSE) %>%
      mutate(
        Run_ID = current_recipe$Run_ID,
        Stage_Type = current_recipe$Stage_Type 
      )
    all_profiles_list[[i]] <- profile_data
  }
  
  if (i %% 200 == 0 || i == nrow(recipes_all_v6)) {
    cat(sprintf("\r已處理 %d / %d 個 profile 檔案...", i, nrow(recipes_all_v6)))
    flush.console()
  }
}

# 高效合併成一個大的 data.frame
all_profiles_df_v6 <- bind_rows(all_profiles_list)
cat("\n--- 所有 V6 profile 數據已整合完畢！ ---\n")


# =========================================================================
#           PART 3: 根據 Stage_Type 分類繪圖
# =========================================================================

# 函式：用於生成標準化的繪圖
create_spaghetti_plot <- function(data, title, subtitle, color) {
  ggplot(data, aes(x = depth_mm, y = concentration_pct, group = Run_ID)) +
    geom_line(alpha = 0.1, color = color) +
    labs(
      title = title,
      subtitle = paste(n_distinct(data$Run_ID), "條曲線"),
      x = "深度 (mm)", y = "碳濃度 (%)"
    ) +
    theme_bw(base_size = 16) 
}

# 3a. 繪製【單階段】
data_single <- filter(all_profiles_df_v6, Stage_Type == "Single_Stage")
plot_single <- create_spaghetti_plot(data_single, "V6 資料庫：所有【單階段】製程曲線", "類型: Single_Stage", "dodgerblue")
ggsave(file.path(output_dir, "1_V6_Single_Stage.png"), plot_single, width = 8, height = 6)
cat("已儲存：單階段曲線圖\n")

# 3b. 繪製【雙階段】
data_dual <- filter(all_profiles_df_v6, Stage_Type == "Dual_BD")
plot_dual <- create_spaghetti_plot(data_dual, "V6 資料庫：所有【雙階段 B-D】製程曲線", "類型: Dual_BD", "forestgreen")
ggsave(file.path(output_dir, "2_V6_Dual_BD.png"), plot_dual, width = 8, height = 6)
cat("已儲存：雙階段曲線圖\n")

# 3c. 繪製【三階段 B-D-B】
data_bdb <- filter(all_profiles_df_v6, Stage_Type == "Triple_BDB")
plot_bdb <- create_spaghetti_plot(data_bdb, "V6 資料庫：所有【三階段 B-D-B】製程曲線", "類型: Triple_BDB", "darkorange")
ggsave(file.path(output_dir, "3_V6_Triple_BDB.png"), plot_bdb, width = 8, height = 6)
cat("已儲存：三階段 B-D-B 類型曲線圖\n")

# 3d. 繪製【三階段 B-D-D】
data_bdd <- filter(all_profiles_df_v6, Stage_Type == "Triple_BDD")
plot_bdd <- create_spaghetti_plot(data_bdd, "V6 資料庫：所有【三階段 B-D-D】製程曲線", "類型: Triple_BDD", "firebrick")
ggsave(file.path(output_dir, "4_V6_Triple_BDD.png"), plot_bdd, width = 8, height = 6)
cat("已儲存：三階段 B-D-D 類型曲線圖\n")


cat("\n\n=============== 所有 V6 分類視覺化圖檔已生成完畢！ ===============\n")