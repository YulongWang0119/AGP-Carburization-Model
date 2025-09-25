# --- 0. 載入必要的套件 ---
library(ggplot2)
library(dplyr)   # 用於 bind_rows，處理 data.frame 合併
library(purrr)   
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

setup_database <- function(base_path) {
  dir.create(base_path, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(base_path, "profiles"), showWarnings = FALSE)
  dir.create(file.path(base_path, "plot_data"), showWarnings = FALSE)
  
  recipes_path <- file.path(base_path, "recipes.csv")
  if (!file.exists(recipes_path)) {
    # 創建一個包含所有可能欄位的空檔案頭
    empty_df <- data.frame(Run_ID=character(),
                           Profile_File_Path=character(), Plot_Data_File_Path=character())
    write.csv(empty_df, recipes_path, row.names = FALSE)
  }
  
  opt_log_path <- file.path(base_path, "optimization_log.csv")
  if (!file.exists(opt_log_path)) {
    write.csv(data.frame(Iteration=integer(), Timestamp=character(),
                         Acquisition_Function=character(), Max_AF_Value=double(),
                         Selected_Run_ID=character()), opt_log_path, row.names = FALSE)
  }
}

save_run_to_db <- function(base_path, run_id, recipe, simulation_results) {
  
  # 從 list 中提取數據
  final_profile <- simulation_results$final_profile
  plot_data <- simulation_results$plot_data
  
  # 1. 準備並儲存 Recipe 主表
  recipe_list <- list(Run_ID = run_id)
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
  # 讀取現有數據，使用 bind_rows 進行合併，可以處理欄位不匹配的情況
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
  cat(sprintf("Saved all data for Run_ID: %s\n", run_id))
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
#           PART 3 (V3): 定義連續參數空間與生成隨機配方
# ===================================================================

# --- 3a (V3). 定義製程參數的【連續】範圍 ---
param_ranges <- list(
  temp_C    = c(850, 1000), # 溫度範圍
  surface_C = c(0.4, 2.0),  # 碳勢範圍
  duration  = c(3, 15)      # 時間範圍
)

# --- 3b (V3). 生成隨機的單階段配方 ---
# 目標：生成約 1000 筆單階段數據
num_single_stage <- 1000
set.seed(202408) # 使用新的種子

# 使用 data.frame 和 replicate 快速生成隨機參數
single_stage_recipes_v3 <- data.frame(
  temp_C    = round(runif(num_single_stage, min = param_ranges$temp_C[1], max = param_ranges$temp_C[2])), # 溫度取整
  surface_C = round(runif(num_single_stage, min = param_ranges$surface_C[1], max = param_ranges$surface_C[2]), 2), # 碳勢保留2位小數
  duration  = round(runif(num_single_stage, min = param_ranges$duration[1], max = param_ranges$duration[2]), 1)  # 時間保留1位小數
)
cat(sprintf("已生成 %d 組隨機的單階段製程配方。\n", nrow(single_stage_recipes_v3)))


# --- 3c (V3). 生成滿足物理約束的隨機雙階段配方 ---
# 目標：生成約 5000-6000 筆雙階段數據
num_dual_stage <- 5500 # 我們可以設定一個目標值
dual_stage_recipes_v3 <- list() # 用一個 list 來收集合格的配方

cat("正在生成雙階段配方，直到滿足物理約束的數量達標...\n")
# 使用 while 迴圈，直到收集到足夠的合格配方
while(length(dual_stage_recipes_v3) < num_dual_stage) {
  
  # 隨機生成第一階段參數
  s1_temp_C    = round(runif(1, min = param_ranges$temp_C[1], max = param_ranges$temp_C[2]))
  s1_surface_C = round(runif(1, min = param_ranges$surface_C[1], max = param_ranges$surface_C[2]), 2)
  s1_duration  = round(runif(1, min = param_ranges$duration[1], max = param_ranges$duration[2]), 1)
  
  # 隨機生成第二階段參數
  s2_temp_C    = round(runif(1, min = param_ranges$temp_C[1], max = param_ranges$temp_C[2]))
  s2_surface_C = round(runif(1, min = param_ranges$surface_C[1], max = param_ranges$surface_C[2]), 2)
  s2_duration  = round(runif(1, min = param_ranges$duration[1], max = param_ranges$duration[2]), 1)
  
  # 【關鍵】檢查物理約束
  if (s1_surface_C > s2_surface_C) {
    # 如果滿足條件，就將這組配方加入到 list 中
    recipe <- list(
      Stage1 = list(temp_C = s1_temp_C, surface_C = s1_surface_C, duration = s1_duration),
      Stage2 = list(temp_C = s2_temp_C, surface_C = s2_surface_C, duration = s2_duration)
    )
    dual_stage_recipes_v3[[length(dual_stage_recipes_v3) + 1]] <- recipe
  }
  
  # 打印進度
  if (length(dual_stage_recipes_v3) %% 100 == 0) {
    cat(sprintf("\r已收集 %d / %d 筆合格配方...", length(dual_stage_recipes_v3), num_dual_stage))
    flush.console()
  }
}
cat(sprintf("\n已成功生成 %d 組滿足物理約束的隨機雙階段製程配方。\n", length(dual_stage_recipes_v3)))


# ===================================================================
#           PART 4 (V3): 執行模擬並建立 V3 資料庫
# ===================================================================

# --- 4a (V3). 建立單階段資料庫 V3 ---
db_path_single_v3 <- "C:/Users/USER/Desktop/PYCProfessor/CarburizationDB_SingleStage_ProfileV3"
setup_database(db_path_single_v3)

cat(sprintf("\n--- 開始生成 %d 筆單階段數據 (V3) ---\n", nrow(single_stage_recipes_v3)))
for (i in 1:nrow(single_stage_recipes_v3)) {
  run_id <- sprintf("S%04d", i)
  params <- single_stage_recipes_v3[i, ]
  
  recipe <- list(
    list(name = "Process", temp_C = params$temp_C, surface_C = params$surface_C, duration = params$duration)
  )
  
  cat(sprintf("\rRunning %s (%d/%d)...", run_id, i, nrow(single_stage_recipes_v3)))
  flush.console()
  
  sim_results <- run_simulation(recipe_stages = recipe)
  save_run_to_db(base_path = db_path_single_v3, run_id = run_id, recipe = recipe, simulation_results = sim_results)
}
cat("\n--- 單階段資料庫 V3 生成完畢！ ---\n")


# --- 4b (V3). 建立雙階段資料庫 V3 ---
db_path_dual_v3 <- "C:/Users/USER/Desktop/PYCProfessor/CarburizationDB_DualStage_ProfileV3"
setup_database(db_path_dual_v3)

cat(sprintf("\n--- 開始生成 %d 筆雙階段數據 (V3) ---\n", length(dual_stage_recipes_v3)))

for (i in 1:length(dual_stage_recipes_v3)) {
  run_id <- sprintf("D%05d", i)
  
  # 從 list 中提取配方
  s1_params <- dual_stage_recipes_v3[[i]]$Stage1
  s2_params <- dual_stage_recipes_v3[[i]]$Stage2
  
  recipe <- list(
    list(name = "Boost", temp_C = s1_params$temp_C, surface_C = s1_params$surface_C, duration = s1_params$duration),
    list(name = "Diffuse", temp_C = s2_params$temp_C, surface_C = s2_params$surface_C, duration = s2_params$duration)
  )
  
  cat(sprintf("\rRunning %s (%d/%d)...", run_id, i, length(dual_stage_recipes_v3)))
  flush.console()
  
  sim_results <- run_simulation(recipe_stages = recipe)
  save_run_to_db(base_path = db_path_dual_v3, run_id = run_id, recipe = recipe, simulation_results = sim_results)
}
cat("\n--- 雙階段資料庫 V3 生成完畢！ ---\n")