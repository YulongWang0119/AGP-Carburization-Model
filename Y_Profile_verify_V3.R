# ========================================================================
#       Y-PROFILE 混合維度模型Hold-Out Validation
# =========================================================================
# 目的：訓練一個基於特定訓練集的最終模型，並在一個從未見過的、
#       完全獨立的測試集上評估其最終性能。

# --- Step 0: 環境設定與載入函式庫 ---
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
# 載入 project_Y_to_weights的函式庫
source("C:/Users/USER/Desktop/PYCProfessor/定期紀錄/Y_Profile_Function_V1.R") 

# =========================================================================
#           PART 1: 數據分割 (Train/Test Split) 
# =========================================================================
cat("--- PART 1: 正在準備完全獨立的訓練集與最終測試集 ---\n")
set.seed(20251111)

db_full_path <- "C:/Users/USER/Desktop/PYCProfessor/CarburizationDB_DualStage_ProfileV3"
recipes_2stage_all <- read_csv(file.path(db_full_path, "recipes.csv"), show_col_types = FALSE)

db_single_path <- "C:/Users/USER/Desktop/PYCProfessor/CarburizationDB_SingleStage_ProfileV3"
recipes_1stage_all <- read_csv(file.path(db_single_path, "recipes.csv"), show_col_types = FALSE)

# --- 抽樣訓練集 (200筆) ---
num_train_2stage <- 100
num_train_1stage <- 100
train_indices_2stage <- sample(1:nrow(recipes_2stage_all), num_train_2stage)
train_indices_1stage <- sample(1:nrow(recipes_1stage_all), num_train_1stage)
recipes_train_2stage <- recipes_2stage_all[train_indices_2stage, ]
recipes_train_1stage <- recipes_1stage_all[train_indices_1stage, ]
#sample(...)：從這個索引序列中，不重複地隨機抽取 num_train_2stage（100）個索引
#recipes_2stage_all[train_indices_2stage, ]：使用這些行號，從完整的主索引表中提取出對應的 100 行，是雙階段訓練集配方

# --- 抽樣最終測試集 (50筆) ---
remaining_indices_2stage <- setdiff(1:nrow(recipes_2stage_all), train_indices_2stage)
remaining_indices_1stage <- setdiff(1:nrow(recipes_1stage_all), train_indices_1stage)
num_test_2stage <- 25; num_test_1stage <- 25
test_indices_2stage <- sample(remaining_indices_2stage, num_test_2stage)
test_indices_1stage <- sample(remaining_indices_1stage, num_test_1stage)
recipes_test_2stage <- recipes_2stage_all[test_indices_2stage, ]
recipes_test_1stage <- recipes_1stage_all[test_indices_1stage, ]
#setdiff(A, B)：會回傳一個在集合 A 中，但不在集合 B 中的所有元素。

cat(sprintf("訓練集抽樣完成: %d 筆雙階段 + %d 筆單階段 = %d 總數\n", num_train_2stage, num_train_1stage, num_train_2stage + num_train_1stage))
cat(sprintf("最終測試集抽樣完成: %d 筆雙階段 + %d 筆單階段 = %d 總數\n", num_test_2stage, num_test_1stage, num_test_2stage + num_test_1stage))

# --- 輔助函式 ---
#recipes_2s：雙階段製程配方
#recipes_1s：單階段製程配方
#db_path_2s：雙階段數據庫的資料夾路徑
#db_path_1s：單階段數據庫的資料夾路徑
prepare_data <- function(recipes_2s, recipes_1s, db_path_2s, db_path_1s) {
  
  X_2stage <- as.matrix(recipes_2s[, c("Stage1_Temp_C", "Stage1_Time_hr", "Stage1_Cs_pct", "Stage2_Temp_C", "Stage2_Time_hr", "Stage2_Cs_pct")])
  X_1stage_3col <- as.matrix(recipes_1s[, c("Stage1_Temp_C", "Stage1_Time_hr", "Stage1_Cs_pct")])
  X_1stage <- matrix(NA, nrow = nrow(recipes_1s), ncol = 6)
  #matrix(NA, ...)：建立一個新的矩陣，行數為 100，欄數為 6。矩陣裡的所有元素都先用 NA 
  X_1stage[, 1:3] <- X_1stage_3col
  
  X_raw <- rbind(X_2stage, X_1stage)
  #rbind() 會將兩個或多個矩陣垂直地堆疊起來。
  #將 100 x 6 的 X_2stage 矩陣和 100 x 6 的 X_1stage 矩陣堆疊在一起，形成一個 200 x 6 的大矩陣 X_raw
  
  paths_2stage <- file.path(db_path_2s, recipes_2s$Profile_File_Path)
  paths_1stage <- file.path(db_path_1s, recipes_1s$Profile_File_Path)
  all_paths <- c(paths_2stage, paths_1stage)
  #recipes_2s$Profile_File_Path：從 recipes_2s 這個 data.frame 中，取出名為 Profile_File_Path 的那一欄。
  #這一欄儲存的是每個深度剖面的相對檔案路徑（例如 "profiles/D00001_profile.csv"）。
  #c(...)：將雙階段和單階段的所有完整檔案路徑合併成一個長的字串向量 all_paths。
  #這個向量的順序與 X_raw 的行順序是完全對應的。
  #all_paths 是一個包含了 200 個元素的一維向量。
  #all_paths[1] 存的是第一筆雙階段資料的完整檔案路徑。
  #all_paths[100] 存的是第一百筆雙階段資料的完整檔案路徑。
  
  first_profile <- read_csv(all_paths[1], show_col_types = FALSE)
  N <- nrow(first_profile)
  #N會是101，因為我設定101個網格點
  
  Y_matrix <- matrix(NA, nrow = N, ncol = length(all_paths))
  #length(all_paths)：計算總共有多少個檔案路徑，也就是總樣本數（例如 200）。
  #matrix(NA, nrow = N, ncol = ...)：根據剛剛得到的 N 和總樣本數，預先建立一個 N x n（例如 101 x 200）的空矩陣 Y_matrix，並填滿 NA。
  
  for (i in 1:length(all_paths)) {
    profile_data <- read_csv(all_paths[i], col_names = c("depth_mm", "concentration_pct"), skip = 1, show_col_types = FALSE)
    if (nrow(profile_data) == N) { Y_matrix[, i] <- profile_data$concentration_pct }
  }
  
  return(list(X_raw = X_raw, Y_matrix = Y_matrix, depth = first_profile$depth_mm))
}

train_data <- prepare_data(recipes_train_2stage, recipes_train_1stage, db_full_path, db_single_path)
X_train <- train_data$X_raw
Y_train <- train_data$Y_matrix
test_data <- prepare_data(recipes_test_2stage, recipes_test_1stage, db_full_path, db_single_path)
X_test <- test_data$X_raw
Y_test <- test_data$Y_matrix
cat("--- 訓練集與最終測試集的 X, Y 已準備完成！ ---\n")

#得到了四個物件：
#X_train (200 x 6 矩陣)：訓練模型的輸入特徵。
#Y_train (101 x 200 矩陣)：訓練模型的輸出目標（深度剖面）。
#X_test (50 x 6 矩陣)：測試模型的輸入特徵。
#Y_test (101 x 50 矩陣)：測試模型的標準答案。

# =========================================================================
#           PART 2: 模型建立 (只使用訓練集) 
# =========================================================================
cat("\n--- PART 2: 正在使用 200 筆訓練集建立最終模型 ---\n")

# --- 2a. 對訓練集進行 KL 降維 ---
decomposition_train <- decompose_Y(Y_train, variance_threshold = 0.999)
#輸入：Y_train (一個 101 x 200 的矩陣)，以及一個閾值 0.999。
#輸出：decomposition_train 是一個 list，裡面包含了 KL 降維的所有重要結果：平均曲線、eigenvectors、eigenvalues、降維後的權重，以及主成分數量num_components
J <- decomposition_train$num_components
weights_train <- decomposition_train$weights # 這是 J x n 矩陣
cat(sprintf("訓練集降維完成，選定 %d 個主成分來建立模型。\n", J))
#GaussianProcess這篇文獻的公式 14-16

# --- 2b. 訓練 J 個最終的 AGP 模型 ---
final_agp_models <- list()
agp_bounds <- list(lower = c(rep(1e-2, 6), rep(1e-2, 3)), upper = c(rep(100, 6), rep(20, 3)))
lower_b_log <- log(agp_bounds$lower)
upper_b_log <- log(agp_bounds$upper)
#log(...)：將邊界轉換到對數尺度。因為最佳化演算法在對數尺度上搜索正數參數（如 θ 和 σ²）希望能更高效

# 對 X 進行正規化，並儲存正規化參數供後續使用
X_min <- apply(X_train, 2, min, na.rm = TRUE)
X_max <- apply(X_train, 2, max, na.rm = TRUE)
X_scale_factors <- X_max - X_min
X_scale_factors[X_scale_factors == 0 | is.na(X_scale_factors)] <- 1 # 避免除以零
X_train_scaled <- scale(X_train, center = X_min, scale = X_scale_factors)
#將原始的 X_train 進行 Min-Max 正規化，結果 X_train_scaled 的每一欄數值都會落在 [0, 1] 區間內。

#迴圈將為每一個權重（每一個主成分）訓練一個專屬的 AGP 模型。
for (j in 1:J) {
  cat(sprintf("正在訓練權重 %d 的最終模型...\n", j))
  start_time <- Sys.time()
  
  # 提取第 j 個權重，並在此處進行正規化
  #weights_train[j, ]：從 J x 200 的權重矩陣中，提取出第 j 列。這一列包含了所有 200 筆訓練數據對應的第 j 個權重值。
  y_j_raw <- as.numeric(weights_train[j, ]) # 正確提取第 j 列權重
  y_j_mean <- mean(y_j_raw)
  y_j_sd <- sd(y_j_raw)
  y_j_scaled <- (y_j_raw - y_j_mean) / y_j_sd
  #進行了標準化，使得 y_j_scaled 這個新的目標向量的平均值為 0，標準差為 1。希望 GP 模型的數值穩定性
  
  model_j <- train_AGP_model(
    X_train = X_train_scaled,
    Y_train_scaled = matrix(y_j_scaled, ncol=1),
    lower_b_log = lower_b_log,
    upper_b_log = upper_b_log,
    y_mean = y_j_mean, # 傳入正規化參數
    y_sd = y_j_sd
  )
  #將剛剛計算出的權重的均值和標準差，一併傳入並儲存到模型物件中，這樣在後續預測時才能正確地將結果還原。
 
  final_agp_models[[j]] <- model_j
  end_time <- Sys.time()
  cat(sprintf("權重 %d 模型訓練完成，耗時: %s\n", j, format(end_time - start_time)))
}
cat("--- 所有 J 個最終模型已訓練完成！ ---\n")

# =========================================================================
#           PART 3: 在獨立測試集上評估)
# =========================================================================
cat("\n--- PART 3: 正在對 50 筆從未見過的測試集進行最終審判 ---\n")

# --- 3a. 預測權重 ---
# 預測前，用訓練集的參數來正規化 X_test
X_test_scaled <- scale(X_test, center = X_min, scale = X_scale_factors)
#使用訓練階段 (PART 2) 計算出的 X_min 和 X_scale_factors，來對測試集的輸入 X_test 進行正規化。
#X_test_scaled 是一個 50 x 6 的矩陣，它的尺度與訓練時使用的 X_train_scaled 一致

weights_pred <- matrix(NA, nrow = nrow(X_test), ncol = J)
#預先分配一個 50 x J 的空矩陣 weights_pred，用來存放模型對 50 筆測試樣本的 J 個權重的預測值
# 使用 for 迴圈逐點預測，而不是 lapply
for (i in 1:nrow(X_test_scaled)) {
  for (j in 1:J) {
    weights_pred[i, j] <- predict_AGP(X_test_scaled[i, ], final_agp_models[[j]])$mean
  }
}
# predict_AGP 是來自 Forrester的 Kriging 預測公式 (2.40), p. 60
#雙層迴圈。外層迴圈 i 遍歷 50 筆測試樣本，內層迴圈 j 遍歷 J 個權重模型。
#predict_AGP(X_test_scaled[i, ], final_agp_models[[j]])：在迴圈的核心，呼叫 predict_AGP 函式，用第 j 個訓練好的模型 (final_agp_models[[j]])，去預測第 i 筆測試樣本 (X_test_scaled[i, ]) 對應的第 j 個權重值。
#$mean：從預測結果（一個包含 mean 和 sd 的 list）中，只提取預測的平均值，也就是模型給出的最佳猜測。
#weights_pred[i, j] <- ...：將預測出的權重值，填入 weights_pred 矩陣的第 i 行、第 j 欄。

colnames(weights_pred) <- paste0("pred_W", 1:J)
#colnames()：為 weights_pred 矩陣的每一欄設定一個有意義的名稱，
#例如 "pred_W1", "pred_W2" 等。讓後續的資料處理和繪圖更方便。

# --- 3b. 計算真實權重 (標準答案) ---
# 確保呼叫的是修正後、帶有 eigenvalues 參數的 project_Y_to_weights
weights_true <- project_Y_to_weights(
  Y_matrix_new = Y_test,
  Y_mean_curve = decomposition_train$Y_mean_curve,
  eigenvectors = decomposition_train$eigenvectors,
  eigenvalues = decomposition_train$eigenvalues 
)
colnames(weights_true) <- paste0("true_W", 1:J)
#Y_matrix_new = Y_test：傳入測試集的真實深度剖面矩陣。

# --- 3c. 計算並報告各權重的效能  ---
results_df <- data.frame()
for (j in 1:J) {
  pred <- weights_pred[, j]
  true <- weights_true[, j]
  r_squared <- cor(pred, true)^2
  rmse <- sqrt(mean((pred - true)^2))
  cat(sprintf("\n--- 【權重 %d】最終測試效能 ---", j)); cat(sprintf("\nR-squared: %.4f", r_squared)); cat(sprintf("\nRMSE: %.4f\n", rmse))
  results_df <- rbind(results_df, data.frame(Weight = j, R2 = r_squared, RMSE = rmse))
}
#sqrt(mean((pred - true)^2))：計算均方根誤差 (RMSE)，衡量的是預測值和真實值的平均差距。越小越好。

# --- 3d. 視覺化比較 ---
plot_data <- as.data.frame(cbind(weights_true, weights_pred))
plot_w1 <- ggplot(plot_data, aes(x = true_W1, y = pred_W1)) +
  geom_point(alpha = 0.7, color = "#377EB8", size = 3) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(title = "最終模型性能：權重 1 (在獨立測試集上)", x = "真實權重 1", y = "預測權重 1", subtitle = paste0("R² = ", round(results_df$R2[1], 4))) +
  theme_bw(base_size = 14)
print(plot_w1)

if (J > 1) {
  plot_w2 <- ggplot(plot_data, aes(x = true_W2, y = pred_W2)) +
    geom_point(alpha = 0.7, color = "#4DAF4A", size = 3) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(title = "最終模型性能：權重 2 (在獨立測試集上)", x = "真實權重 2", y = "預測權重 2", subtitle = paste0("R² = ", round(results_df$R2[2], 4))) +
    theme_bw(base_size = 14)
  print(plot_w2)
}

cat("\n\n=============== 留出驗證 (Hold-Out Validation) 完成 ===============\n")

# =========================================================================
#           PART 4: 實驗結果彙總與儲存
# =========================================================================
cat("\n--- PART 4: 正在生成並儲存本次實驗的彙總圖檔 ---\n")

# --- 4a.建立一個指定路徑下的儲存結果資料夾 ---

# 1. 定義主儲存路徑
base_output_path <- "C:/Users/USER/Desktop/PYCProfessor/Picture/Y_Profile"

# 2. 確保這個主路徑存在，如果不存在就自動創建它
dir.create(base_output_path, showWarnings = FALSE, recursive = TRUE)
#recursive = TRUE：如果中間層的資料夾（例如 Picture）也不存在，請一併自動建立。

# 3. 結合主路徑和帶有時間戳的資料夾名稱，形成最終的儲存路徑
output_dir <- file.path(base_output_path, 
                        paste0("HoldOut_Results", format(Sys.time(), "%Y%m%d_%H%M%S")))

# 4. 創建這個最終的、帶有時間戳的資料夾
dir.create(output_dir)

cat(sprintf("所有結果將儲存於資料夾: %s\n", output_dir))

# --- 4b. 繪製並儲存全局變化模式】 ---
# 這個圖基於訓練集，展示了模型的基礎是如何建立的
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

# --- 4c (雙版本). 繪製並儲存【所有權重的預測 vs. 真實圖】 ---
# 這是本次實驗最核心的效能展示
# 我們將為每個權重生成兩種版本的圖：一個包含R²，一個只含RMSE

# 我們不再需要一個外部的 plots_list，可以在迴圈內直接處理
for (j in 1:J) {
  
  # --- 版本一：同時包含 R² 和 RMSE ---
  
  subtitle_v1 <- sprintf("R² = %.4f | RMSE = %.4f", results_df$R2[j], results_df$RMSE[j])
  
  plot_v1 <- ggplot(plot_data, aes_string(x = paste0("true_W", j), y = paste0("pred_W", j))) +
    geom_point(alpha = 0.7, color = scales::hue_pal()(J)[j], size = 3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(title = paste("最終模型性能：權重", j, "(在50筆獨立測試集上)"),
         x = paste("真實權重", j), y = paste("預測權重", j),
         subtitle = subtitle_v1) +
    theme_bw(base_size = 14)
  
  # 在 RStudio 中顯示版本一的圖
  print(plot_v1)
  
  # 儲存版本一的圖，並在檔名中標註
  file_name_v1 <- sprintf("B_Weight_%d_Performance_R2_RMSE.png", j)
  ggsave(file.path(output_dir, file_name_v1), plot_v1, width = 8, height = 6)
  
  
  # --- 版本二：只包含 RMSE ---
  
  subtitle_v2 <- sprintf("RMSE = %.4f", results_df$RMSE[j])
  
  plot_v2 <- ggplot(plot_data, aes_string(x = paste0("true_W", j), y = paste0("pred_W", j))) +
    geom_point(alpha = 0.7, color = scales::hue_pal()(J)[j], size = 3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(title = paste("最終模型性能：權重", j, "(在50筆獨立測試集上)"),
         x = paste("真實權重", j), y = paste("預測權重", j),
         subtitle = subtitle_v2) +
    theme_bw(base_size = 14)
  
  # 儲存版本二的圖，並在檔名中標註
  file_name_v2 <- sprintf("B_Weight_%d_Performance_RMSE_only.png", j)
  ggsave(file.path(output_dir, file_name_v2), plot_v2, width = 8, height = 6)
  
  
  # --- 打印完成訊息 ---
  cat(sprintf("已為權重 %d 生成並儲存 R2+RMSE 和 僅RMSE 兩種版本的性能圖。\n", j))
}

# --- PART 5: 單點預測、重構與貢獻度分析 ---
cat("\n--- PART 5: 正在進行單一樣本的完整重構分析 ---\n")

# 1. 挑選一個測試樣本 (例如，第 10 筆測試數據)
sample_idx <- 10
x_sample_raw <- X_test[sample_idx, ]
y_sample_true <- Y_test[, sample_idx]

cat("正在為以下輸入參數進行預測與分析：\n")
print(x_sample_raw)

# 2. 使用增強版模型進行預測
prediction_results <- predict_and_reconstruct_profile(
  x_new_raw = x_sample_raw,
  model_list = final_agp_models,
  decomposition_results = decomposition_train,
  X_min = X_min, 
  X_scale_factors = X_scale_factors
)

# 提取所有需要的資訊
y_sample_pred <- prediction_results$final_prediction
predicted_weights <- prediction_results$predicted_weights
mean_curve <- prediction_results$mean_curve
mode_contributions_df <- prediction_results$mode_contributions
colnames(mode_contributions_df) <- paste0("貢獻_模式", 1:J)


# 3. 視覺化完整重構過程
# 準備繪圖用的資料框
depth_vector <- train_data$depth 
comparison_df_long <- data.frame(
  depth = depth_vector,
  `真實 Profile` = y_sample_true,
  `最終預測 Profile` = y_sample_pred,
  `平均曲線` = mean_curve
) %>% 
  bind_cols(mode_contributions_df) %>%
  pivot_longer(
    cols = -depth,
    names_to = "曲線類型",
    values_to = "濃度"
  )

# 自定義顏色和線型
line_colors <- c("真實 Profile" = "black", 
                 "最終預測 Profile" = "#E41A1C", # 紅色
                 "平均曲線" = "#377EB8",       # 藍色
                 "貢獻_模式1" = "#4DAF4A",      # 綠色
                 "貢獻_模式2" = "#984EA3",      # 紫色
                 "貢獻_模式3" = "#FF7F00",      # 橘色
                 "貢獻_模式4" = "#F781BF")      # 粉色
line_types <- c("真實 Profile" = "dashed", 
                "最終預測 Profile" = "solid", 
                "平均曲線" = "dotted",
                "貢獻_模式1" = "longdash",
                "貢獻_模式2" = "longdash",
                "貢獻_模式3" = "longdash",
                "貢獻_模式4" = "longdash")
line_widths <- c("真實 Profile" = 1.2, 
                 "最終預測 Profile" = 1.5, 
                 "平均曲線" = 1.0,
                 "貢獻_模式1" = 0.8,
                 "貢獻_模式2" = 0.8,
                 "貢獻_模式3" = 0.8,
                 "貢獻_模式4" = 0.8)


# 建立一個標題來顯示預測的權重係數
weights_text <- paste(sprintf("ξ%d=%.2f", 1:J, predicted_weights), collapse = ", ")
plot_title <- paste("單一樣本預測重構全覽 (測試樣本 #", sample_idx, ")")
plot_subtitle <- paste("預測權重:", weights_text)

# 繪圖
plot_reconstruction <- ggplot(comparison_df_long, aes(x = depth, y = 濃度, color = 曲線類型)) +
  geom_line(aes(linetype = 曲線類型, linewidth = 曲線類型)) +
  labs(
    title = plot_title,
    subtitle = plot_subtitle,
    x = "深度 (mm)",
    y = "濃度 (%) / 貢獻量"
  ) +
  scale_color_manual(values = line_colors) +
  scale_linetype_manual(values = line_types) +
  scale_linewidth_manual(values = line_widths) +
  theme_bw(base_size = 14) +
  guides(linewidth = "none") # 不顯示線寬的圖例

print(plot_reconstruction)

# 4. 儲存圖片到指定的資料夾
if (!exists("output_dir")) {
  # 如果是獨立執行，提供一個預設路徑
  output_dir <- file.path("C:/Users/USER/Desktop/PYCProfessor/Picture/Y_Profile", 
                          paste0("Single_Prediction_Results", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}

# 儲存重構圖
reconstruction_filename <- file.path(output_dir, sprintf("C_Reconstruction_Sample_%d.png", sample_idx))
ggsave(reconstruction_filename, plot_reconstruction, width = 12, height = 7)
cat(sprintf("已儲存單一樣本的完整重構圖至: %s\n", reconstruction_filename))


# 計算並儲存文字檔，包含 RMSE 和權重
single_point_rmse <- sqrt(mean((y_sample_true - y_sample_pred)^2))
info_text <- c(
  paste("測試樣本編號:", sample_idx),
  paste("輸入參數:", paste(round(x_sample_raw, 4), collapse = ", ")),
  paste("預測曲線的 RMSE:", round(single_point_rmse, 5)),
  "\n--- 預測的權重係數 (ξ_hat) ---",
  weights_text
)
info_filename <- file.path(output_dir, sprintf("C_Info_Sample_%d.txt", sample_idx))
writeLines(info_text, info_filename)
cat(sprintf("已儲存單一樣本的詳細資訊至: %s\n", info_filename))

# =========================================================================
#           PART 6: 生成並儲存簡潔版預測 vs. 真實對比圖 (附RMSE)
# =========================================================================
cat("\n--- PART 6: 正在生成簡潔版預測對比圖 ---\n")

# 繼續使用 PART 5 中已計算好的 sample_idx, y_sample_true, y_sample_pred

# 計算這條曲線的 RMSE
single_point_rmse <- sqrt(mean((y_sample_true - y_sample_pred)^2))
#mean(...) 或 (1/N) * Σ[...]：計算誤差平方的平均值 (Mean Squared Error, MSE)。
#sqrt(...)：最後再開根號，將單位從「濃度的平方」還原回「濃度」
#這樣是誤差的單位會和原始數據的單位一致，更容易直觀地理解誤差的大小
# 準備繪圖用的資料框
comparison_df_simple <- data.frame(
  depth = depth_vector,
  `真實曲線` = y_sample_true,
  `模型預測` = y_sample_pred
) %>%
  tidyr::pivot_longer(
    cols = -depth,
    names_to = "圖例",
    values_to = "濃度"
  )

# 將 RMSE 加入圖表的副標題
plot_subtitle_simple <- sprintf("整條曲線的 RMSE = %.5f", single_point_rmse)

# 繪圖
plot_simple_comparison <- ggplot(comparison_df_simple, aes(x = depth, y = 濃度, color = 圖例)) +
  geom_line(aes(linetype = 圖例), linewidth = 1.2) + # 讓線型也跟著圖例變
  labs(
    title = paste("最終模型預測 vs. 真實深度剖面對比 (測試樣本 #", sample_idx, ")"),
    subtitle = plot_subtitle_simple, # 使用新的副標題
    x = "深度 (mm)",
    y = "濃度 (%)"
  ) +
  scale_color_manual(name = "圖例", values = c("真實曲線" = "black", "模型預測" = "red")) +
  scale_linetype_manual(name = "圖例", values = c("真實曲線" = "dashed", "模型預測" = "solid")) +
  theme_bw(base_size = 16) + # 放大字體讓圖更清晰
  theme(legend.position = "bottom") # 將圖例放在底部

print(plot_simple_comparison)

# 儲存圖片到與 RMSE 檔案相同的資料夾
simple_comparison_filename <- file.path(output_dir, sprintf("D_Simple_Comparison_Sample_%d.png", sample_idx))
ggsave(simple_comparison_filename, plot_simple_comparison, width = 10, height = 6)
cat(sprintf("已儲存簡潔版對比圖至: %s\n", simple_comparison_filename))

cat(sprintf("\n此單一樣本【整條預測曲線】的 RMSE 為: %.5f\n", single_point_rmse))