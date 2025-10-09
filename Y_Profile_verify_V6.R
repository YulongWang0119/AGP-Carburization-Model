# ========================================================================
#       Y-PROFILE 混合維度模型 Hold-Out Validation (V6)
# =========================================================================
# 目的：訓練一個基於特定訓練集的最終模型，並在一個從未見過的、
#       完全獨立的測試集上評估其最終性能。

# --- Step 0: 環境設定與載入函式庫 ---
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr) 

source("C:/Users/USER/Desktop/PYCProfessor/定期紀錄/Y_Profile_Function_V6.R") 

# =========================================================================
#           PART 1: 手動指定樣本數的分層抽樣 (V6 手動控制版)
# =========================================================================
cat("--- PART 1: 正在從 V6 資料庫，根據手動設定準備訓練/測試集 ---\n")

# --- 1a. 設定區---
db_main_path_v6 <- "C:/Users/USER/Desktop/PYCProfessor/CarburizationDB_ProfileV6"
RANDOM_SEED <- 20308965

# --- 指定抽樣數量---
n_train_s1    <- 75      # 單階段 (Single_Stage) 訓練數量
n_train_s2_bd <- 75     # 雙階段 (Dual_BD) 訓練數量
n_train_s3_bdb<- 75      # 三階段 (Triple_BDB) 訓練數量 像一階
n_train_s3_bdd<- 75      # 三階段 (Triple_BDD) 訓練數量

n_test_s1     <- 25      # 單階段 (Single_Stage) 測試數量
n_test_s2_bd  <- 25      # 雙階段 (Dual_BD) 測試數量
n_test_s3_bdb <- 25      # 三階段 (Triple_BDB) 測試數量 像一階
n_test_s3_bdd <- 25      # 三階段 (Triple_BDD) 測試數量

# --- 1b. 讀取 V6 recipes ---
set.seed(RANDOM_SEED)

recipes_path_v6 <- file.path(db_main_path_v6, "recipes.csv")
if (!file.exists(recipes_path_v6)) stop("找不到 recipes.csv！請檢查路徑。")
recipes_all <- read_csv(recipes_path_v6, show_col_types = FALSE)
# show_col_types = FALSE 讀取時不要顯示每個欄位的型別

# 建立一個包含使用者設定的 data.frame
sampling_plan <- tribble(
  ~Stage_Type,  ~n_train, ~n_test,
  "Single_Stage", n_train_s1,    n_test_s1,
  "Dual_BD",      n_train_s2_bd, n_test_s2_bd,
  "Triple_BDB",   n_train_s3_bdb,n_test_s3_bdb,
  "Triple_BDD",   n_train_s3_bdd,n_test_s3_bdd
)

# 驗證設定是否合理
db_counts <- recipes_all %>% count(Stage_Type, name = "db_total")
#把recipes_all表格放到流水線上
# count(Stage_Type)`: dplyr 去計算 `recipes_all` 裡面，每個 `Stage_Type` 各有多少筆資料
# `name = "db_total"`: 計算出來的數量欄位命名為 `db_total`
validation <- sampling_plan %>%
  left_join(db_counts, by = "Stage_Type") %>% # 把抽樣計畫表和資料庫總數表根據 Stage_Type 合併
  mutate(required = n_train + n_test) %>% ## 新增一個欄位 `required`，計算每個類別共需要多少筆資料
  filter(required > db_total)

if (nrow(validation) > 0) {
  print(validation)
  stop("錯誤：您指定的訓練+測試數量，超過了資料庫中該類別的總數！請檢查設定。")
}

cat("使用者設定驗證通過。\n")
cat("抽樣計畫如下：\n"); print(sampling_plan)


# --- 1c. 根據手動設定執行分層抽樣 ---
recipes_by_type <- split(recipes_all, recipes_all$Stage_Type)
#`recipes_by_type` 會變成一個list，裡面有 "Single_Stage" 的表格、"Dual_BD" 的表格...

# 為訓練集抽樣
# imap會對列表裡的每一個元素和它的名字做同一件事
train_recipes_list <- imap(recipes_by_type, function(df, type_name) {
  n_to_sample <- sampling_plan$n_train[sampling_plan$Stage_Type == type_name]
  # 去 `sampling_plan` 表格裡，找到對應類別 (`type_name`) 的訓練樣本數 (`n_train`)
  slice_sample(df, n = n_to_sample) # 從目前小表格 (`df`) 中，隨機抽出 `n_to_sample` 筆資料
})
recipes_train <- bind_rows(train_recipes_list)
# 把抽出來的好幾個小表格，再合併成一個大的訓練集表格 `recipes_train`

# 為測試集抽樣 (從剩餘數據中)
#  `filter`的條件是：`Run_ID` 不在 (`! ... %in% ...`)剛剛抽出的訓練集 `recipes_train` 裡
remaining_recipes <- filter(recipes_all, !Run_ID %in% recipes_train$Run_ID)
remaining_by_type <- split(remaining_recipes, remaining_recipes$Stage_Type)

test_recipes_list <- imap(remaining_by_type, function(df, type_name) {
  n_to_sample <- sampling_plan$n_test[sampling_plan$Stage_Type == type_name]
  slice_sample(df, n = n_to_sample)
})
recipes_test <- bind_rows(test_recipes_list)

cat("\n手動分層抽樣完成，最終結果：\n")
cat("訓練集實際分佈:\n"); print(count(recipes_train, Stage_Type))
cat("測試集實際分佈:\n"); print(count(recipes_test, Stage_Type))

# 最終驗證獨立性
intersection_count <- length(intersect(recipes_train$Run_ID, recipes_test$Run_ID))
if (intersection_count == 0) {
  cat("\nSUCCESS: 訓練集與測試集完全獨立！\n")
} else {
  stop("\nERROR: 訓練集與測試集存在重疊，抽樣邏輯有誤！")
}


# --- 1d. 通用函式，準備 X 和 Y ---
prepare_model_data_v6 <- function(recipes_df, base_path) {
  
  # --- 準備 X ---
  X_matrix <- recipes_df %>%
    select(starts_with("Stage")) %>% # 挑出所有以 "Stage" 開頭的欄位
    select(-"Stage_Type") %>%        # 再從中排除 "Stage_Type" 這個文字欄位
    mutate_all(as.numeric) %>%       # 把所有剩下的欄位都轉換成數值格式
    replace(is.na(.), 0) %>%         # 如果有 NA ，把它們都換成 0
    as.matrix()                      # 最後轉換成 matrix 格式
  rownames(X_matrix) <- recipes_df$Run_ID # 把矩陣的每一列，用 Run_ID 來命名
  
  # --- 準備 Y ---
  profile_paths <- file.path(base_path, recipes_df$Profile_File_Path) # 根據 recipes 表格裡的檔案路徑，組合成每一筆資料對應的濃度曲線 CSV 檔的完整路徑
  Y_list <- map(profile_paths, ~read_csv(.x, show_col_types = FALSE)$concentration_pct) 
  # `map 會對 `profile_paths` 裡的每一個路徑做同一件事
  # 依序讀取每一個濃度曲線 CSV 檔，只取出 `concentration_pct` 這一欄
  # `Y_list` 會是一個 list，裡面裝著很多條濃度曲線（每一條都是一個向量）
  
  Y_matrix_temp <- do.call(cbind, Y_list)
  # `cbind` 是把很多個向量按行合併成一個大矩陣。`do.call` 是一種呼叫函式的方式
  # 這個矩陣的每一行是一筆實驗，每一列是一個深度點。 (n x N)
  Y_matrix <- t(Y_matrix_temp) #(N x n)
  rownames(Y_matrix) <- recipes_df$Run_ID # 同樣用 Run_ID 命名
  
  depth_vector <- read_csv(profile_paths[1], show_col_types = FALSE)$depth_mm
  #  讀取任何一個檔案，可以知道深度的座標軸是什麼
  return(list(X = X_matrix, Y = Y_matrix, depth = depth_vector))
}


# --- 1e. 生成最終的訓練集與測試集 ---
train_data <- prepare_model_data_v6(recipes_train, db_main_path_v6)
X_train <- train_data$X
Y_train_rows <- train_data$Y # Y_train_rows 是 N x n 格式
Y_train <- t(Y_train_rows)   # 轉置回 n x N 格式，KL 函式 decompose_Y 需要這個格式

test_data <- prepare_model_data_v6(recipes_test, db_main_path_v6)
X_test <- test_data$X
Y_test_rows <- test_data$Y
Y_test <- t(Y_test_rows)

cat("\n--- 訓練集與最終測試集的 X, Y 已準備完成！ ---\n")
cat(sprintf("X_train: %d x %d | Y_train: %d x %d (N x n)\n", nrow(X_train), ncol(X_train), nrow(Y_train), ncol(Y_train)))
cat(sprintf("X_test: %d x %d | Y_test: %d x %d (N x n)\n", nrow(X_test), ncol(X_test), nrow(Y_test), ncol(Y_test)))

# =========================================================================
#           PART 2: 模型建立 (只使用訓練集) 
# =========================================================================
# 訓練集總數為 75+75+75 = 300 筆
cat(sprintf("\n--- PART 2: 正在使用 %d 筆訓練集建立最終模型 ---\n", nrow(X_train)))

# --- 2a. 對訓練集進行 KL 降維 ---
decomposition_train <- decompose_Y(Y_train, variance_threshold = 0.999)
# 找出幾條最關鍵的基本變化模式 eigenvectors），用來組合出所有原始的濃度曲線
# `variance_threshold = 0.999` 的意思是：選取足夠多的基本模式，直到它們能解釋原始資料中 99.9% 的變異。
J <- decomposition_train$num_components
weights_train <- decomposition_train$weights # 這是 J x n 矩陣
# `weights` 是每一條原始曲線，分別需要多少比例的模式1、多少比例的模式」...來組合
# 原本要預測一整條曲線 (很多個點)，現在問題被簡化成：只要預測 J 個權重值就好
cat(sprintf("訓練集降維完成，選定 %d 個主成分來建立模型。\n", J))

# --- 2b. 訓練 J 個最終的 AGP 模型 ---
final_agp_models <- list() # 準備一個空的列表，用來存放 J 個訓練好的模型

# 定義 16 個超參數的邊界 (9 thetas + 7 sigmas)
agp_bounds <- list(
  lower = c(rep(1e-2, 9), rep(1e-2, 7)),
  upper = c(rep(100, 9),  rep(20, 7))
)
lower_b_log <- log(agp_bounds$lower)
upper_b_log <- log(agp_bounds$upper)

# 對 X 進行正規化，並儲存正規化參數供後續使用
X_min <- apply(X_train, 2, min, na.rm = TRUE)
X_max <- apply(X_train, 2, max, na.rm = TRUE)
X_scale_factors <- X_max - X_min
X_scale_factors[X_scale_factors == 0 | is.na(X_scale_factors)] <- 1 # 避免除以零
# 把這些正規化參數 (min, max) 存起來，因為等等預測測試集時必須用同一套標準來正規化測試集
X_train_scaled <- scale(X_train, center = X_min, scale = X_scale_factors)
# `scale()` 函式會做正規化，公式是：(X - center) / scale。
# 就是把所有 X 特徵都縮放到 [0, 1] 的範圍內

# 迴圈將為每一個權重（每一個主成分）訓練一個專屬的 AGP 模型
for (j in 1:J) {
  cat(sprintf("正在訓練權重 %d 的最終模型...\n", j))
  start_time <- Sys.time()
  
  # 提取第 j 個權重，並在此處進行正規化
  y_j_raw <- as.numeric(weights_train[j, ])  # 從權重矩陣中，拿出第 j 列，也就是所有樣本的第 j 個權重
  y_j_mean <- mean(y_j_raw)
  y_j_sd <- sd(y_j_raw)
  y_j_scaled <- (y_j_raw - y_j_mean) / y_j_sd
  
  model_j <- train_AGP_model(
    X_train = X_train_scaled,
    Y_train_scaled = matrix(y_j_scaled, ncol=1), # 目標是正規化後的第 j 個權重
    lower_b_log = lower_b_log,
    upper_b_log = upper_b_log,
    y_mean = y_j_mean, # 傳入正規化參數
    y_sd = y_j_sd
  )
  
  final_agp_models[[j]] <- model_j # 把訓練好的模型，存到列表的第 j 個位置。
  end_time <- Sys.time()
  cat(sprintf("權重 %d 模型訓練完成，耗時: %s\n", j, format(end_time - start_time)))
}
cat("--- 所有 J 個最終模型已訓練完成！ ---\n")

# 儲存完整的模型物件
model_to_save <- list(
  models = final_agp_models,
  decomp = decomposition_train,
  X_min = X_min,
  X_scale_factors = X_scale_factors
)
# 建立儲存資料夾
dir.create("C:/Users/USER/Desktop/PYCProfessor/Picture/Y_Profile/Saved_Models", showWarnings = FALSE, recursive = TRUE)
# 儲存模型
saveRDS(model_to_save, "C:/Users/USER/Desktop/PYCProfessor/Picture/Y_Profile/Saved_Models/My_Model_V6.rds")
cat("--- 最終模型已儲存至 .rds 檔案！ ---\n")


# =========================================================================
#           PART 3: 在獨立測試集上評估
# =========================================================================
# 測試集總數為 25+25+25+25 = 100 筆
cat(sprintf("\n--- PART 3: 正在對 %d 筆從未見過的測試集進行最終評估 ---\n", nrow(X_test)))

# --- 3a. 預測權重 ---
# 預測前，用訓練集的參數來正規化 X_test
X_test_scaled <- scale(X_test, center = X_min, scale = X_scale_factors)

weights_pred <- matrix(NA, nrow = nrow(X_test), ncol = J) # 準備一個空的矩陣來放預測結果
# 使用 for 迴圈逐點預測
for (i in 1:nrow(X_test_scaled)) {
  for (j in 1:J) {
    weights_pred[i, j] <- predict_AGP(final_agp_models[[j]], X_test_scaled[i, ])$mean
  }
}
colnames(weights_pred) <- paste0("pred_W", 1:J)

# --- 3b. 計算真實權重 (標準答案) ---
weights_true <- project_Y_to_weights(
  Y_matrix_new = Y_test,
  Y_mean_curve = decomposition_train$Y_mean_curve,
  eigenvectors = decomposition_train$eigenvectors,
  eigenvalues = decomposition_train$eigenvalues 
)
colnames(weights_true) <- paste0("true_W", 1:J)

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

# --- 3e.對整個測試集進行曲線重構並計算 RMSE 和 MAPE ---
cat("\n--- 正在對整個測試集進行曲線重構以計算性能指標分佈 ---\n")

# 準備空的向量來儲存每一條曲線的 RMSE 和 MAPE
all_test_rmse <- numeric(nrow(X_test))
all_test_mape <- numeric(nrow(X_test))
epsilon <- 1e-6 # 用於防止除以零

for (i in 1:nrow(X_test)) {
  # 重構預測曲線
  reconstructed_curve <- decomposition_train$Y_mean_curve
  for (j in 1:J) {
    reconstructed_curve <- reconstructed_curve + 
      weights_pred[i, j] * decomposition_train$eigenvectors[, j] * sqrt(decomposition_train$eigenvalues[j])
  }
  
  # 獲取真實曲線
  true_curve <- Y_test[, i]
  
  # 計算並儲存 RMSE
  all_test_rmse[i] <- sqrt(mean((reconstructed_curve - true_curve)^2))
  
  # 計算並儲存 MAPE
  percentage_errors <- abs((reconstructed_curve - true_curve) / (true_curve + epsilon))
  all_test_mape[i] <- mean(percentage_errors) * 100 # 以百分比形式儲存
}

# 將 RMSE 和 MAPE 結果與 recipes_test 合併
results_test_df <- recipes_test %>%
  mutate(
    curve_rmse = all_test_rmse,
    curve_mape = all_test_mape
  )

cat("已計算完所有 100 筆測試樣本的曲線 RMSE 和 MAPE。\n")
cat("--- RMSE 統計摘要 ---\n"); print(summary(results_test_df$curve_rmse))
cat("--- MAPE 統計摘要 ---\n"); print(summary(results_test_df$curve_mape))

# =========================================================================
#           PART 4: 實驗結果彙總與儲存
# =========================================================================
cat("\n--- PART 4: 正在生成並儲存本次實驗的彙總圖檔 ---\n")

# --- 4a. 建立一個指定路徑下的儲存結果資料夾 ---
base_output_path <- "C:/Users/USER/Desktop/PYCProfessor/Picture/Y_Profile"
dir.create(base_output_path, showWarnings = FALSE, recursive = TRUE)
output_dir <- file.path(base_output_path, 
                        paste0("HoldOut_Results_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(output_dir)
cat(sprintf("所有結果將儲存於資料夾: %s\n", output_dir))

# --- 4b. 繪製並儲存全局變化模式 ---
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
  labs(title = sprintf("模型基礎：基於%d筆訓練集的平均曲線與變化模式", nrow(X_train)),
       subtitle = "這是後續所有預測與權重計算的基準",
       x = "深度 (mm)", y = "濃度 / 變化量", color = "曲線") +
  scale_color_manual(name = "曲線",
                     values = c("平均曲線" = "#E41A1C", "模式1" = "#377EB8", "模式2" = "#4DAF4A", 
                                "模式3" = "#984EA3", "模式4" = "#FF7F00", "模式5" = "#F781BF"),
                     breaks = c("平均曲線", paste0("模式", 1:J))) +
  theme_bw(base_size = 14)

print(plot_train_modes)
ggsave(file.path(output_dir, "A_Training_Decomposition_Modes.png"), plot_train_modes, width = 10, height = 6)
cat("已儲存：模型基礎變化模式圖\n")

# --- 4b-2. 計算並報告特徵值與變異貢獻率 ---
cat("\n--- KL 展開結果分析 ---\n")

# 從降維結果中提取特徵值
eigenvalues <- decomposition_train$eigenvalues

# 計算每個模式解釋的變異比例
variance_explained <- eigenvalues / sum(eigenvalues)
cumulative_variance <- cumsum(variance_explained)

# 建立一個清晰的報告表格
eigen_report <- data.frame(
  模式 = paste0("模式 ", 1:J),
  特徵值_lambda = eigenvalues,
  變異貢獻率 = variance_explained,
  累積貢獻率 = cumulative_variance
)

cat("各變化模式的重要性分析 (基於訓練集):\n")
print(eigen_report, row.names = FALSE)

# 將這份報告儲存為文字檔
eigen_report_filename <- file.path(output_dir, "A_Eigenvalue_Report.txt")
# 使用 sink() 將 print() 的輸出重導向到檔案
sink(eigen_report_filename)
cat("各變化模式的重要性分析 (基於訓練集):\n\n")
print(eigen_report, row.names = FALSE)
sink() # 結束重導向
cat(sprintf("已將特徵值分析報告儲存至: %s\n", eigen_report_filename))

# --- 4c. 繪製並儲存所有權重的預測 vs. 真實圖 ---
for (j in 1:J) {
  # 版本一：同時包含 R² 和 RMSE
  subtitle_v1 <- sprintf("R² = %.4f | RMSE = %.4f", results_df$R2[j], results_df$RMSE[j])
  plot_v1 <- ggplot(plot_data, aes_string(x = paste0("true_W", j), y = paste0("pred_W", j))) +
    geom_point(alpha = 0.7, color = scales::hue_pal()(J)[j], size = 3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(title = sprintf("最終模型性能：權重 %d (在%d筆獨立測試集上)", j, nrow(X_test)),
         x = paste("真實權重", j), y = paste("預測權重", j),
         subtitle = subtitle_v1) +
    theme_bw(base_size = 14)
  print(plot_v1)
  
  file_name_v1 <- sprintf("B_Weight_%d_Performance_R2_RMSE.png", j)
  ggsave(file.path(output_dir, file_name_v1), plot_v1, width = 8, height = 6)
  
  # 版本二：只包含 RMSE
  subtitle_v2 <- sprintf("RMSE = %.4f", results_df$RMSE[j])
  plot_v2 <- ggplot(plot_data, aes_string(x = paste0("true_W", j), y = paste0("pred_W", j))) +
    geom_point(alpha = 0.7, color = scales::hue_pal()(J)[j], size = 3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(title = sprintf("最終模型性能：權重 %d (在%d筆獨立測試集上)", j, nrow(X_test)),
         x = paste("真實權重", j), y = paste("預測權重", j),
         subtitle = subtitle_v2) +
    theme_bw(base_size = 14)
  
  file_name_v2 <- sprintf("B_Weight_%d_Performance_RMSE_only.png", j)
  ggsave(file.path(output_dir, file_name_v2), plot_v2, width = 8, height = 6)
  
  cat(sprintf("已為權重 %d 生成並儲存 R2+RMSE 和 僅RMSE 兩種版本的性能圖。\n", j))
}


# ---  找出 RMSE 或 MAPE 最高的幾個案例 ---
cat("\n--- 預測 RMSE 最差的 5 個樣本 ---\n")
worst_rmse_predictions <- results_test_df %>%
  arrange(desc(curve_rmse)) %>%
  head(5)
print(worst_rmse_predictions)

cat("\n--- 預測 MAPE 最差的 5 個樣本 ---\n")
worst_mape_predictions <- results_test_df %>%
  arrange(desc(curve_mape)) %>%
  head(5)
print(worst_mape_predictions)



# =========================================================================
#           PART 5 & 6: 迴圈處理多個單點預測
# =========================================================================
cat("\n--- PART 5 & 6: 正在對指定的多個樣本進行完整分析 ---\n")

#indices_to_analyze <- c(1, 9) # 可以換成任何想看的樣本編號
# 自動鎖定最差的 5 個樣本進行分析
worst_indices <- match(worst_rmse_predictions$Run_ID, recipes_test$Run_ID)
indices_to_analyze <- worst_indices

# 開始迴圈
for (sample_idx in indices_to_analyze) {
  
# =========================================================================
#           PART 5: 單點預測、重構與貢獻度分析
# =========================================================================
cat("\n--- PART 5: 正在進行單一樣本的完整重構分析 ---\n")

# 可以更改 sample_idx 來分析不同的測試樣本 (範圍: 1 到 nrow(X_test))

x_sample_raw <- X_test[sample_idx, ]
y_sample_true <- Y_test[, sample_idx]
# 獲取並顯示樣本的真實類型
sample_info <- recipes_test[sample_idx, ]
sample_type <- sample_info$Stage_Type
cat(sprintf("\n分析的樣本 #%d 屬於 [%s] 類型。\n", sample_idx, sample_type))

cat("正在為以下輸入參數進行預測與分析：\n"); print(na.omit(x_sample_raw))

# 使用模型進行預測
prediction_results <- predict_and_reconstruct_profile(
  x_new_raw = x_sample_raw,
  model_list = final_agp_models,
  decomposition_results = decomposition_train,
  X_min = X_min, 
  X_scale_factors = X_scale_factors
)

# <<< 請在這裡加入以下診斷碼 >>>
cat("\n\n=============== 診斷開始 ===============\n")
cat("--- 每個權重預測的標準差 (SD) ---\n")
print(prediction_results$predicted_sds)
cat("--- 最終預測曲線的標準差 (SD) 的前5個值 ---\n")
print(head(prediction_results$final_sd_curve, 5))
cat("=============== 診斷結束 ===============\n\n")
# <<< 診斷碼結束 >>>

# 提取所有需要的資訊
y_sample_pred <- prediction_results$final_prediction
predicted_weights <- prediction_results$predicted_weights
mean_curve <- prediction_results$mean_curve
mode_contributions_df <- prediction_results$mode_contributions
colnames(mode_contributions_df) <- paste0("貢獻_模式", 1:J)

y_sample_pred_sd <- prediction_results$final_sd_curve 

# 計算 95% 信賴區間
upper_ci_bound <- y_sample_pred + 1.96 * y_sample_pred_sd
lower_ci_bound <- y_sample_pred - 1.96 * y_sample_pred_sd

# --- 整合權重與特徵值的分析 ---
# 從降維結果中提取特徵值和變異貢獻率
eigenvalues <- decomposition_train$eigenvalues
variance_explained <- eigenvalues / sum(eigenvalues)

# 建立一個包含所有分析資訊的 data.frame
analysis_df <- data.frame(
  模式 = 1:J,
  預測權重_xi = predicted_weights,
  模式重要性_lambda = eigenvalues,
  變異貢獻率_pct = variance_explained * 100
)

# 計算每個模式對最終曲線 RMSE 的貢獻
# 是近似值，因為模式之間是正交的
mode_rmse_contribution <- numeric(J)
for(j in 1:J){
  # 這裡需要真實權重來做比較
  true_weights_sample <- project_Y_to_weights(
    matrix(y_sample_true, ncol=1),
    decomposition_train$Y_mean_curve,
    decomposition_train$eigenvectors,
    decomposition_train$eigenvalues
  )
  weight_error_j <- predicted_weights[j] - true_weights_sample[j]
  
  # 誤差傳遞: (權重誤差 * 模式範數 * sqrt(特徵值))^2
  mode_rmse_contribution[j] <- (weight_error_j * sqrt(eigenvalues[j]))^2
}
analysis_df$對RMSE平方的貢獻 <- mode_rmse_contribution / sum(mode_rmse_contribution) * 100


cat(sprintf("\n--- 對測試樣本 #%d 的詳細分析 ---\n", sample_idx))
print(analysis_df)

# =========================================================
#診斷程式碼
cat("\n--- 權重不確定性診斷 ---\n")
cat("每個權重預測的標準差 (SD):\n")
print(prediction_results$predicted_sds)
# =========================================================

# 將分析表格儲存
analysis_filename <- file.path(output_dir, sprintf("C_Detailed_Analysis_Sample_%d.txt", sample_idx))
sink(analysis_filename)
cat(sprintf("--- 對測試樣本 #%d 的詳細分析 ---\n\n", sample_idx))
cat("輸入參數:\n"); print(na.omit(x_sample_raw)); cat("\n")
print(analysis_df)
sink()
cat(sprintf("已將詳細分析儲存至: %s\n", analysis_filename))

# 使用 check.names = FALSE 阻止 R 自動更改欄位名
comparison_df_long <- data.frame(
  depth = depth_vector,
  `真實 Profile` = y_sample_true,
  `最終預測 Profile` = y_sample_pred,
  `平均曲線` = mean_curve,
  check.names = FALSE 
) %>% 
  bind_cols(as.data.frame(mode_contributions_df)) %>%
  pivot_longer(
    cols = -depth,
    names_to = "曲線類型",
    values_to = "濃度"
  )

# =========================================================
# 診斷程式碼
cat(sprintf("\n--- 診斷資訊 (樣本 #%d) ---\n", sample_idx))
cat("`comparison_df_long` 中的 `曲線類型` 包含以下幾種：\n")
print(unique(comparison_df_long$曲線類型))
cat("----------------------------------\n\n")
# =========================================================

# --- 動態生成顏色和線型 ---
# 基礎曲線的設定
base_colors <- c("真實 Profile" = "black", "最終預測 Profile" = "#E41A1C", "平均曲線" = "#377EB8")
base_types <- c("真實 Profile" = "dashed", "最終預測 Profile" = "solid", "平均曲線" = "dotted")
base_widths <- c("真實 Profile" = 1.2, "最終預測 Profile" = 1.5, "平均曲線" = 1.0)

# 動態生成模式曲線的設定
mode_names <- paste0("貢獻_模式", 1:J)
# 使用 ggplot 內建的調色盤，確保顏色夠用且不同
mode_colors <- scales::hue_pal()(J) 
names(mode_colors) <- mode_names

mode_types <- rep("longdash", J)
names(mode_types) <- mode_names

mode_widths <- rep(0.8, J)
names(mode_widths) <- mode_names

# 合併成最終的設定
line_colors <- c(base_colors, mode_colors)
line_types <- c(base_types, mode_types)
line_widths <- c(base_widths, mode_widths)

# 建立一個標題來顯示預測的權重係數
weights_text <- paste(sprintf("ξ%d=%.2f", 1:J, predicted_weights), collapse = ", ")
# 在標題中加入類型資訊
plot_title <- sprintf("單一樣本預測重構全覽 (測試樣本 #%d, 類型: %s)", sample_idx, sample_type)
plot_subtitle <- paste("預測權重:", weights_text)

plot_reconstruction <- ggplot(comparison_df_long, aes(x = depth, y = 濃度, color = 曲線類型)) +
  geom_line(aes(linetype = 曲線類型, linewidth = 曲線類型)) +
  labs(title = plot_title, subtitle = plot_subtitle, x = "深度 (mm)", y = "濃度 (%) / 貢獻量") +
  # 移除了 breaks 参数，让 ggplot 自动匹配
  scale_color_manual(values = line_colors) +
  scale_linetype_manual(values = line_types) +
  scale_linewidth_manual(values = line_widths) +
  theme_bw(base_size = 14) +
  guides(linewidth = "none")

print(plot_reconstruction)

# 儲存圖片到指定的資料夾
reconstruction_filename <- file.path(output_dir, sprintf("C_Reconstruction_Sample_%d.png", sample_idx))
ggsave(reconstruction_filename, plot_reconstruction, width = 12, height = 7)
cat(sprintf("已儲存單一樣本的完整重構圖至: %s\n", reconstruction_filename))

# 計算並儲存文字檔
single_point_rmse <- sqrt(mean((y_sample_true - y_sample_pred)^2))
info_text <- c(
  paste("測試樣本編號:", sample_idx),
  paste("輸入參數:", paste(na.omit(round(x_sample_raw, 4)), collapse = ", ")),
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

# 計算 RMSE
single_point_rmse <- sqrt(mean((y_sample_true - y_sample_pred)^2))

epsilon <- 1e-6
single_point_mape <- mean(abs((y_sample_true - y_sample_pred) / (y_sample_true + epsilon))) * 100
plot_subtitle_simple <- sprintf("整條曲線的 RMSE = %.5f | MAPE = %.2f%%", single_point_rmse, single_point_mape)

# 將信賴區間加入 data frame
comparison_df_simple_with_ci <- data.frame(
  depth = depth_vector,
  真實曲線 = y_sample_true,
  模型預測 = y_sample_pred,
  信賴上限 = upper_ci_bound,
  信賴下限 = lower_ci_bound
) %>%
  tidyr::pivot_longer(cols = c("真實曲線", "模型預測"), names_to = "圖例", values_to = "濃度")

# 繪圖 (加入 geom_ribbon)
plot_simple_comparison <- ggplot(comparison_df_simple_with_ci, aes(x = depth)) +
  # 繪製 95% 信賴區間的陰影
  geom_ribbon(aes(ymin = 信賴下限, ymax = 信賴上限), fill = "purple", alpha = 0.4) +
  
  geom_line(aes(y = 濃度, color = 圖例, linetype = 圖例), linewidth = 1.2) +
  labs(title = sprintf("最終模型預測 vs. 真實 (樣本 #%d, 類型: %s)", sample_idx, sample_type),
       subtitle = plot_subtitle_simple, x = "深度 (mm)", y = "濃度 (%)") +
  scale_color_manual(name = "圖例", values = c("真實曲線" = "black", "模型預測" = "red")) +
  scale_linetype_manual(name = "圖例", values = c("真實曲線" = "dashed", "模型預測" = "solid")) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")

print(plot_simple_comparison)

# 儲存圖片
simple_comparison_filename <- file.path(output_dir, sprintf("D_Simple_Comparison_Sample_%d.png", sample_idx))
ggsave(simple_comparison_filename, plot_simple_comparison, width = 10, height = 6)
cat(sprintf("已儲存簡潔版對比圖至: %s\n", simple_comparison_filename))

cat(sprintf("\n此單一樣本【整條預測曲線】的 RMSE 為: %.5f\n", single_point_rmse))

} # 迴圈結束
# =========================================================================
#           PART 7: 繪製測試集性能指標分佈圖 (RMSE & MAPE)
# =========================================================================
cat("\n--- PART 7: 正在繪製 RMSE 和 MAPE 的分佈直方圖 ---\n")

# --- 繪圖 A: RMSE 直方圖 ---
mean_rmse <- mean(results_test_df$curve_rmse)
plot_rmse_histogram <- ggplot(results_test_df, aes(x = curve_rmse)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.001, fill = "#377EB8", color = "white", alpha = 0.8) +
  geom_density(color = "red", linewidth = 1) +
  geom_vline(xintercept = mean_rmse, color = "darkgreen", linetype = "dashed", linewidth = 1.2) +
  annotate("text", x = mean_rmse, y = Inf, label = paste(" 平均 RMSE\n", round(mean_rmse, 5)), 
           color = "darkgreen", hjust = -0.1, vjust = 1.5, angle = 90) +
  labs(
    title = sprintf("模型在 %d 筆獨立測試集上的 RMSE 分佈", nrow(results_test_df)),
    x = "整條曲線的 RMSE", y = "密度 (Density)"
  ) +
  theme_bw(base_size = 14)
print(plot_rmse_histogram)

histogram_filename_rmse <- file.path(output_dir, "E_RMSE_Histogram.png")
ggsave(histogram_filename_rmse, plot_rmse_histogram, width = 10, height = 6)
cat(sprintf("已儲存 RMSE 直方圖至: %s\n", histogram_filename_rmse))

# --- 繪圖 B: MAPE 直方圖 ---
mean_mape <- mean(results_test_df$curve_mape)
plot_mape_histogram <- ggplot(results_test_df, aes(x = curve_mape)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, fill = "#E41A1C", color = "white", alpha = 0.8) +
  geom_density(color = "blue", linewidth = 1) +
  geom_vline(xintercept = mean_mape, color = "darkgreen", linetype = "dashed", linewidth = 1.2) +
  annotate("text", x = mean_mape, y = Inf, label = paste(" 平均 MAPE\n", round(mean_mape, 2), "%"), 
           color = "darkgreen", hjust = -0.1, vjust = 1.5, angle = 90) +
  labs(
    title = sprintf("模型在 %d 筆獨立測試集上的 MAPE 分佈", nrow(results_test_df)),
    x = "整條曲線的 MAPE (%)", y = "密度 (Density)"
  ) +
  theme_bw(base_size = 14)
print(plot_mape_histogram)

histogram_filename_mape <- file.path(output_dir, "F_MAPE_Histogram.png")
ggsave(histogram_filename_mape, plot_mape_histogram, width = 10, height = 6)
cat(sprintf("已儲存 MAPE 直方圖至: %s\n", histogram_filename_mape))


# ---  找出 RMSE 或 MAPE 最高的幾個案例 ---
cat("\n--- 預測 RMSE 最差的 5 個樣本 ---\n")
worst_rmse_predictions <- results_test_df %>%
  arrange(desc(curve_rmse)) %>%
  head(5)
print(worst_rmse_predictions)

cat("\n--- 預測 MAPE 最差的 5 個樣本 ---\n")
worst_mape_predictions <- results_test_df %>%
  arrange(desc(curve_mape)) %>%
  head(5)
print(worst_mape_predictions)

# 根據 Stage_Type 分組，計算每一種類型的平均 MAPE
results_test_df %>%
  group_by(Stage_Type) %>%
  summarise(
    avg_mape = mean(curve_mape),
    count = n()
  )

# 繪製按製程類型分的 MAPE 小提琴圖
plot_mape_by_type <- ggplot(results_test_df, aes(x = Stage_Type, y = curve_mape, fill = Stage_Type)) +
  geom_violin(trim = FALSE, alpha = 0.8) + # 繪製小提琴圖
  geom_jitter(width = 0.1, alpha = 0.5, height = 0) + # 加上實際的數據點
  labs(
    title = "模型在不同製程類型上的 MAPE 表現",
    subtitle = "證實了不同複雜度的製程對預測精度的影響",
    x = "製程類型 (Stage Type)",
    y = "整條曲線的 MAPE (%)"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") # x軸已有標示，無需圖例

print(plot_mape_by_type)

# 儲存這張圖
ggsave(file.path(output_dir, "G_MAPE_by_StageType.png"), plot_mape_by_type, width = 10, height = 6)