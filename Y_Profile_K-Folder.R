# =========================================================================
#             Y-PROFILE 混合維度模型最嚴謹交叉驗證腳本
# =========================================================================

# --- Step 0: 環境設定與載入函式庫 ---
library(readr); library(ggplot2); library(dplyr); library(caret); library(tidyr)
# 載入您包含所有函式的 Y_Profile_Function_V1.R
source("C:/Users/USER/Desktop/PYCProfessor/定期紀錄/Y_Profile_Function_V1.R") 

# =========================================================================
#           PART 1: 準備完整的 200 筆混合真實數據集
# =========================================================================
cat("--- PART 1: 準備 200 筆混合真實數據集 ---\n")
set.seed(202406)

# --- 1a. 抽樣雙階段數據 ---
db_full_path <- "C:/Users/USER/Desktop/PYCProfessor/CarburizationDB_DualStage_ProfileV2"
recipes_full_df <- read_csv(file.path(db_full_path, "recipes.csv"), show_col_types = FALSE)
num_2stage_sample <- 150
indices_2stage <- sample(1:nrow(recipes_full_df), num_2stage_sample)
recipes_2stage_sampled <- recipes_full_df[indices_2stage, ]
cat(sprintf("已從雙階段資料庫中隨機抽取 %d 筆數據。\n", num_2stage_sample))

# --- 1b. 抽樣單階段數據 ---
db_single_path <- "C:/Users/USER/Desktop/PYCProfessor/CarburizationDB_SingleStage_ProfileV2"
recipes_single_df <- read_csv(file.path(db_single_path, "recipes.csv"), show_col_types = FALSE)
num_1stage_sample <- 50
indices_1stage <- sample(1:nrow(recipes_single_df), num_1stage_sample)
recipes_1stage_sampled <- recipes_single_df[indices_1stage, ]
cat(sprintf("已從單階段資料庫中隨機抽取 %d 筆數據。\n", num_1stage_sample))

# --- 1c. 整合 X_raw 和 Y_matrix ---
X_2stage <- as.matrix(recipes_2stage_sampled[, c("Stage1_Temp_C", "Stage1_Time_hr", "Stage1_Cs_pct", "Stage2_Temp_C", "Stage2_Time_hr", "Stage2_Cs_pct")])
X_1stage_3col <- as.matrix(recipes_1stage_sampled[, c("Stage1_Temp_C", "Stage1_Time_hr", "Stage1_Cs_pct")])
X_1stage <- matrix(NA, nrow = num_1stage_sample, ncol = 6)
X_1stage[, 1:3] <- X_1stage_3col
X_raw_combined <- rbind(X_2stage, X_1stage)

paths_2stage <- file.path(db_full_path, recipes_2stage_sampled$Profile_File_Path)
paths_1stage <- file.path(db_single_path, recipes_1stage_sampled$Profile_File_Path)
all_profile_paths <- c(paths_2stage, paths_1stage)
first_profile <- read_csv(all_profile_paths[1], show_col_types = FALSE)
N <- nrow(first_profile)
Y_matrix_combined <- matrix(NA, nrow = N, ncol = nrow(X_raw_combined))
for (i in 1:nrow(X_raw_combined)) {
  profile_data <- read_csv(all_profile_paths[i], col_names = c("depth_mm", "concentration_pct"), skip = 1, show_col_types = FALSE)
  if (nrow(profile_data) == N) { Y_matrix_combined[, i] <- profile_data$concentration_pct }
}
cat("--- 完整的 200 筆 X 和 Y_matrix 準備完成！---\n")


# =========================================================================
#           PART 2: 僅供探索與報告說明用的全局降維與視覺化
# =========================================================================
cat("\n--- PART 2: 正在執行僅供探索用的全局降維 ---\n")
cat("--- (此部分的結果將不會用於後續的嚴謹驗證) ---\n")

# 執行一次全局降維，以產生用於解釋方法論的圖表
global_decomposition <- decompose_Y(Y_matrix_combined, variance_threshold = 0.999)
J_global <- global_decomposition$num_components

# 繪製全局的變化模式圖
depth_vector <- first_profile$depth_mm
plot_df_base_global <- data.frame(depth = depth_vector, mean_curve = global_decomposition$Y_mean_curve)
plot_df_wide_global <- plot_df_base_global %>% bind_cols(as.data.frame(global_decomposition$eigenvectors))
colnames(plot_df_wide_global)[-c(1,2)] <- paste0("模式", 1:J_global)
plot_df_tidy_global <- plot_df_wide_global %>% pivot_longer(cols = starts_with("模式"), names_to = "模式", values_to = "變化量")
plot_all_modes_global <- ggplot() +
  geom_line(data = plot_df_base_global, aes(x = depth, y = mean_curve, color = "平均曲線"), linewidth = 1.5) +
  geom_line(data = plot_df_tidy_global, aes(x = depth, y = 變化量, color = 模式), linetype = "dashed", linewidth = 1) +
  labs(title = "平均滲碳曲線與所有主要變化模式 (基於全部200筆數據)", 
       x = "深度 (mm)", y = "濃度 / 變化量", color = "曲線") +
  scale_color_manual(name = "曲線",
                     values = c("平均曲線" = "#E41A1C", "模式1" = "#377EB8", "模式2" = "#4DAF4A", 
                                "模式3" = "#984EA3", "模式4" = "#FF7F00", "模式5" = "#FFFF33"),
                     breaks = c("平均曲線", paste0("模式", 1:J_global))) +
  theme_bw(base_size = 14)

print(plot_all_modes_global)
cat("--- 全局探索性分析與繪圖完成 ---\n")

# =========================================================================
#           PART 3: 運行最嚴謹的交叉驗證
# =========================================================================
cat("\n--- PART 3: 開始對各權重模型進行最嚴謹的5-Fold 交叉驗證 ---\n")
# --- 3a. 設定參數 ---
# 不再使用 J_global，因為每一折的 J 可能都不同。
# 我們設定一個合理的驗證上限，例如，我們只關心前 4 個權重的表現。
J_to_validate <- 4 
strict_cv_results_list <- list()

data_source_labels <- c(rep("2-stage", num_2stage_sample), 
                        rep("1-stage", num_1stage_sample))

agp_bounds <- list(
  lower = c(rep(1e-2, 6), rep(1e-2, 3)), 
  upper = c(rep(100, 6), rep(20, 3))
)

# --- 3b. 迴圈，對每一個權重進行最嚴謹的交叉驗證 ---
for (j in 1:J_to_validate) {
  start_time_cv <- Sys.time()
  
  result <- validate_gp_model(
    X_raw                = X_raw_combined, 
    Y_matrix_raw         = Y_matrix_combined,
    j_weight_to_validate = j,
    model_train_func     = train_AGP_model,
    model_predict_func   = predict_AGP,
    hyper_params_bounds  = agp_bounds,
    k_folds              = 5,
    seed                 = 2024,
    strata_labels        = data_source_labels,
    decompose_func       = decompose_Y
  )
  
  # 加入一個檢查，如果 result 是 NULL (代表出錯或所有折都跳過)，就處理一下
  if(!is.null(result)) {
    strict_cv_results_list[[j]] <- result
    end_time_cv <- Sys.time()
    
    cat(sprintf("\n--- 【權重 %d】最嚴謹交叉驗證效能 ---", j))
    print(paste("R-squared:", round(result$r_squared, 4)))
    print(paste("RMSE:", round(result$rmse, 4)))
    print(paste("交叉驗證總耗時:", format(end_time_cv - start_time_cv)))
  } else {
    cat(sprintf("\n--- 無法完成【權重 %d】的驗證，可能是因為在某些折中主成分數量不足。 ---\n", j))
  }
}

cat("\n\n=============== 所有權重模型【最嚴謹】驗證完成 ===============\n")

#---------------------------------------------------------------------------
# install.packages("gridExtra")
# =========================================================================
#           PART 4: 實驗結果彙總與儲存
# =========================================================================
cat("\n--- PART 4: 正在生成並儲存本次實驗的彙總圖檔 ---\n")

# --- 4a. 建立一個儲存結果的資料夾 ---
output_dir <- paste0("HoldOut_Results_", format(Sys.time(), "%Y%m%d_%H%M%S"))
dir.create(output_dir)
cat(sprintf("所有結果將儲存於資料夾: %s\n", output_dir))

# --- 4b. 繪製並儲存全局變化模式】 ---
# 這個圖基於訓練集，展示了模型的基礎如何建立
depth_vector <- train_data$depth
plot_df_base_train <- data.frame(depth = depth_vector, mean_curve = decomposition_train$Y_mean_curve)
plot_df_wide_train <- plot_df_base_train %>%
  bind_cols(as.data.frame(decomposition_train$eigenvectors))
colnames(plot_df_wide_train)[-c(1,2)] <- paste0("模式", 1:J)

plot_df_tidy_train <- plot_df_wide_train %>%
  pivot_longer(cols = starts_with("模式"), names_to = "模式", values_to = "變化量")

plot_train_modes <- ggplot() +
  geom_line(data = plot_df_base_train, aes(x = depth, y = mean_curve, color = "平均曲線"), linewidth = 1.5) +
  geom_line(data = plot_df_tidy_train, aes(x = depth, y = 變化量, color = 模式), linetype = "dashed", linewidth = 1) +
  labs(title = "模型基礎：基於200筆訓練集的平均曲線與變化模式",
       subtitle = "這是後續所有預測與權重計算的基準",
       x = "深度 (mm)", y = "濃度 / 變化量", color = "曲線") +
  scale_color_manual(name = "曲線",
                     values = c("平均曲線" = "#E41A1C", "模式1" = "#377EB8", "模式2" = "#4DAF4A", 
                                "模式3" = "#984EA3", "模式4" = "#FF7F00"),
                     breaks = c("平均曲線", paste0("模式", 1:J))) +
  theme_bw(base_size = 14)

print(plot_train_modes)
ggsave(file.path(output_dir, "A_Training_Decomposition_Modes.png"), plot_train_modes, width = 10, height = 6)
cat("已儲存：模型基礎變化模式圖\n")


# --- 4c. 繪製並儲存所有權重的預測 vs. 真實圖---
plots_list <- list()
for (j in 1:J) {
  # 為 subtitle 準備文字
  subtitle_text <- sprintf("R² = %.4f | RMSE = %.4f", results_df$R2[j], results_df$RMSE[j])
  
  # 動態選擇 x 和 y 的欄位
  plot_j <- ggplot(plot_data, aes_string(x = paste0("true_W", j), y = paste0("pred_W", j))) +
    geom_point(alpha = 0.7, color = scales::hue_pal()(J)[j], size = 3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(title = paste("最終模型性能：權重", j, "(在50筆獨立測試集上)"),
         x = paste("真實權重", j), y = paste("預測權重", j),
         subtitle = subtitle_text) +
    theme_bw(base_size = 14)
  
  plots_list[[j]] <- plot_j # 存入 list
}

# 儲存個別權重的圖檔
for (j in 1:J) {
  file_name <- sprintf("B_Weight_%d_Performance.png", j)
  ggsave(file.path(output_dir, file_name), plots_list[[j]], width = 8, height = 6)
  cat(sprintf("已儲存：權重 %d 性能圖\n", j))
}

# 
# # --- 4d. 使用 gridExtra 套件將所有權重圖合併在一張大圖上 ---
# if (requireNamespace("gridExtra", quietly = TRUE)) {
#   plot_all_weights_grid <- gridExtra::grid.arrange(grobs = plots_list, ncol = 2) # 排成 2x2 的網格
#   
#   # 為合併後的大圖添加一個總標題
#   final_grid_with_title <- gridExtra::grid.arrange(
#     top = grid::textGrob("所有權重模型在獨立測試集上的最終性能彙總", gp=grid::gpar(fontsize=20)),
#     plot_all_weights_grid
#   )
#   
#   ggsave(file.path(output_dir, "C_All_Weights_Performance_Grid.png"), final_grid_with_title, width = 12, height = 10)
#   cat("已儲存：所有權重性能合併圖\n")
# } else {
#   cat("\n提示：若要生成所有權重性能的合併圖，請先安裝 'gridExtra' 套件 ( install.packages('gridExtra') )。\n")
# }
# 
# cat("\n--- 所有圖檔生成與儲存完畢！ ---\n")