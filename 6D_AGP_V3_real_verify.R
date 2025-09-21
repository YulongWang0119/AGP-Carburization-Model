# # install.packages("caret") 
# library(caret) # 載入 caret 套件
# # K-Fold 交叉驗證
# 
# print("--- 開始執行 5-Fold 交叉驗證 ---")
# 
# # 0. 參數設定
# set.seed(20250910) 
# K <- 5         
# 
# # =========================================================================
# # 在 createFolds 之前先打亂資料順序，避免因資料排序導致切分出標準差為零的訓練集
# # =========================================================================
# cat("--- 正在打亂原始資料順序以增加交叉驗證的隨機性 ---\n")
# shuffle_indices <- sample(1:nrow(X_train_raw))
# X_train_raw_shuffled <- X_train_raw[shuffle_indices, ]
# Y_train_raw_shuffled <- Y_train_raw[shuffle_indices, , drop = FALSE]
# # =========================================================================
# folds <- createFolds(1:nrow(X_train_raw_shuffled), k = K, list = TRUE, returnTrain = FALSE)
# 
# 
# # 用於儲存每一折預測結果的容器
# all_predictions <- data.frame(real_Y = numeric(), predicted_Y = numeric())
# all_rmses <- numeric(K) # 新增一個容器來儲存每一折的 RMSE
# 
# # 1. 開始 K-Fold 迴圈
# for (k in 1:K) {
#   
#   cat(sprintf("\n--- 正在執行第 %d / %d 折 ---\n", k, K))
#   
#   # 2. 分割數據
#   # 根據 folds 列表，找出當前這一折的測試集索引和訓練集索引
#   test_indices <- folds[[k]]
#   train_indices <- setdiff(1:nrow(X_train_raw_shuffled), test_indices)
#   
#   # 從打亂後的數據中分割出訓練集和測試集
#   X_cv_train_raw <- X_train_raw_shuffled[train_indices, ]
#   Y_cv_train_raw <- Y_train_raw_shuffled[train_indices, , drop = FALSE]
#   
#   X_cv_test_raw <- X_train_raw_shuffled[test_indices, ]
#   Y_cv_test_raw <- Y_train_raw_shuffled[test_indices, , drop = FALSE]
#   
#   
#   # 3. 在當前訓練集上進行正規化
#   #    正規化因子必須只從訓練集中計算，然後應用到測試集上
#   
#   # 計算 X 的正規化因子
#   X_cv_min <- apply(X_cv_train_raw, 2, min, na.rm = TRUE) # 增加 na.rm = TRUE 以防萬一
#   X_cv_max <- apply(X_cv_train_raw, 2, max, na.rm = TRUE)
#   scale_cv_factors <- X_cv_max - X_cv_min
#   
#   # =========================================================================
#   # 偵錯程式碼區塊 START，診斷每一折的資料特性
#   # =========================================================================
#   cat(sprintf("--- Fold %d 診斷報告 ---\n", k))
#   cat(sprintf("訓練集大小: %d 個樣本\n", length(train_indices)))
#   cat(sprintf("測試集大小: %d 個樣本\n", length(test_indices)))
#   
#   # --- 偵測：Y 的標準差是否為零 ---
#   Y_cv_sd <- sd(Y_cv_train_raw, na.rm = TRUE)
#   cat(sprintf("訓練集 Y 的標準差 (sd): %f\n", Y_cv_sd))
#   if (is.na(Y_cv_sd) || Y_cv_sd < 1e-9) {
#     cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
#     cat(sprintf("警告！第 %d 折的訓練集 Y 標準差為零或極小。\n", k))
#     cat("這會導致除以零的錯誤。正在跳過此折...\n")
#     cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n")
#     all_rmses[k] <- NA # 記錄此折失敗
#     next # 直接跳到下一折
#   }
#   
#   # --- 偵測：X 的某些維度標準差是否為零 ---
#   zero_var_cols <- which(scale_cv_factors < 1e-9)
#   if (length(zero_var_cols) > 0) {
#     cat("--------------------------------------------------------\n")
#     cat(sprintf("警告！第 %d 折的訓練集 X 在以下維度上標準差為零或極小:\n", k))
#     print(colnames(X_cv_train_raw)[zero_var_cols])
#     cat("這可能導致後續超參數優化不穩定。\n")
#     cat("--------------------------------------------------------\n")
#   }
#   # =========================================================================
#   # 【偵錯程式碼區塊 END】
#   # =========================================================================
#   
#   # 防止除以零
#   scale_cv_factors[scale_cv_factors == 0] <- 1
#   
#   # 計算 Y 的正規化因子 (因為上面已經算過，這裡可以直接用)
#   Y_cv_mean <- mean(Y_cv_train_raw, na.rm = TRUE)
#   
#   # =========================================================================
#   # 數據淨化與類型強制
#   # =========================================================================
#   
#   # 對 X_train 進行正規化並確保是純矩陣
#   X_cv_train_scaled <- as.matrix(scale(X_cv_train_raw, center = X_cv_min, scale = scale_cv_factors))
# 
#   # 對 Y_train 進行正規化並確保是純 n x 1 矩陣
#   # 先把 Y_cv_train_raw 轉成向量，再轉回矩陣，確保維度正確
#   Y_cv_train_scaled <- matrix((as.numeric(Y_cv_train_raw) - Y_cv_mean) / Y_cv_sd, ncol = 1)
#   
#   # 對 X_test 進行正規化並確保是純矩陣
#   X_cv_test_scaled <- as.matrix(scale(X_cv_test_raw, center = X_cv_min, scale = scale_cv_factors))
#   
#   # =========================================================================
#   # ... K-fold 迴圈內部 ...
#   # =========================================================================
#   # 在呼叫模型訓練函式之前，將所有輸入資料強制轉換為最純粹的形式
#   # =========================================================================
#   X_final_train <- as.matrix(X_cv_train_scaled)
#   Y_final_train_scaled <- as.matrix(as.vector(Y_cv_train_scaled))
#   
#   # 移除所有維度名稱 (dimnames)
#   dimnames(X_final_train) <- NULL
#   dimnames(Y_final_train_scaled) <- NULL
#   
#   cat("\n--- 最終資料淨化後，準備送入模型訓練 ---\n")
#   cat("X_final_train 的維度: "); print(dim(X_final_train))
#   cat("Y_final_train_scaled 的維度: "); print(dim(Y_final_train_scaled))
#   cat("---------------------------------------------\n\n")
#   # =========================================================================
#   
#   # 4. 在當前的訓練集上訓練 AGP 模型
#   #    現在傳入的是被徹底淨化過的資料
#   cv_model <- train_AGP_model(
#     X_train = X_final_train,            # <--- 使用淨化後的 X
#     Y_train_scaled = Y_final_train_scaled,  # <--- 使用淨化後的 Y
#     lower_b_log = lower_bounds_log, 
#     upper_b_log = upper_bounds_log,
#     y_mean = Y_cv_mean, 
#     y_sd = Y_cv_sd
#   )
#   
#   # 5. 在當前的測試集上進行預測
#   predictions_k <- data.frame(
#     real_Y = as.vector(Y_cv_test_raw),
#     predicted_Y = rep(NA, nrow(X_cv_test_raw))
#   )
#   
#   for (j in 1:nrow(X_cv_test_scaled)) {
#     # 提取一個測試點 (已正規化)
#     x_new_point <- X_cv_test_scaled[j, ]
#     
#     # 進行預測 (預測函式會自動返回【真實尺度】的值)
#     pred_result <- predict_AGP(as.numeric(x_new_point), cv_model)
#     predictions_k$predicted_Y[j] <- pred_result$mean
#   }
#   
#   # 6. 將這一折的預測結果儲存到總容器中
#   all_predictions <- rbind(all_predictions, predictions_k)
#   
#   # 計算並儲存這一折的 RMSE
#   rmse_k <- sqrt(mean((predictions_k$real_Y - predictions_k$predicted_Y)^2, na.rm = TRUE))
#   all_rmses[k] <- rmse_k
#   cat(sprintf("--- 第 %d 折完成，此折的 RMSE: %f ---\n", k, rmse_k))
# }
# 
# 
# # 第六部分計算與展示交叉驗證結果
# 
# print("\n--- 5-Fold 交叉驗證結果匯總 ---")
# 
# # 1. 計算性能指標
# #    從 all_predictions 計算總體指標，並從 all_rmses 計算平均指標
# #    RMSE (均方根誤差)：單位與 Y 相同
# rmse_overall <- sqrt(mean((all_predictions$real_Y - all_predictions$predicted_Y)^2, na.rm = TRUE))
# print(paste("交叉驗證的【總體】RMSE:", round(rmse_overall, 4)))
# 
# rmse_avg <- mean(all_rmses, na.rm = TRUE)
# print(paste("交叉驗證的【平均每折】RMSE:", round(rmse_avg, 4)))
# 
# 
# #    R² ：表示模型能解釋數據變異的百分比，越接近 1 越好
# r_squared <- cor(all_predictions$real_Y, all_predictions$predicted_Y)^2
# print(paste("交叉驗證的 R-squared (決定係數):", round(r_squared, 4)))
# 
# #    MAE ：比 RMSE 對異常值不敏感
# mae <- mean(abs(all_predictions$real_Y - all_predictions$predicted_Y), na.rm = TRUE)
# print(paste("交叉驗證的 MAE (平均絕對誤差):", round(mae, 4)))
# 
# 
# # 2. 繪製真實值 vs. 預測值散點圖
# library(ggplot2)
# if(nrow(all_predictions) > 0) {
#   plot_cv_results <- ggplot(all_predictions, aes(x = real_Y, y = predicted_Y)) +
#     geom_point(alpha = 0.7, color = "blue") +
#     # 畫一條 y=x 的預測線，用於對比
#     geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 1) +
#     # 讓 x 軸和 y 軸的範圍一致，使 y=x 線呈 45 度角
#     coord_fixed(ratio = 1, xlim = range(all_predictions, na.rm = TRUE), ylim = range(all_predictions, na.rm = TRUE)) +
#     labs(
#       title = "5-Fold 交叉驗證結果",
#       subtitle = paste("R² =", round(r_squared, 4), " | 全局 RMSE =", round(rmse_overall, 4)),
#       x = "真實有效深度",
#       y = "模型預測有效深度"
#     ) +
#     theme_minimal()
#   
#   print(plot_cv_results)
# } else {
#   print("沒有任何一折成功，無法繪圖。請檢查診斷報告。")
# }

#-------------------------------------------------------------------------------------------
#以下是模組化
# =========================================================================
# 將 K-Fold 流程封裝成函式
# =========================================================================

library(caret) # 載入 caret 套件
validate_gp_model <- function(X_raw, Y_raw, 
                              model_train_func, model_predict_func,
                              hyper_params_bounds, k_folds = 5, seed = 123) {
  
  # 檢查輸入
  if (!is.function(model_train_func) || !is.function(model_predict_func)) {
    stop("model_train_func 和 model_predict_func 必須是函式。")
  }
  
  print(paste("--- 開始執行", k_folds, "-Fold 交叉驗證 ---"))
  
  set.seed(seed)
  
  # --- 1. 超參數邊界設定 ---
  lower_b_log <- log(hyper_params_bounds$lower)
  upper_b_log <- log(hyper_params_bounds$upper)
  
  # --- 2. 資料預處理與切分 ---
  shuffle_indices <- sample(1:nrow(X_raw))
  X_raw_shuffled <- X_raw[shuffle_indices, ]
  Y_raw_shuffled <- Y_raw[shuffle_indices, , drop = FALSE]
  
  folds <- createFolds(1:nrow(X_raw_shuffled), k = k_folds, list = TRUE, returnTrain = FALSE)
  
  # --- 3. 初始化容器 ---
  all_predictions <- data.frame(real_Y = numeric(), predicted_Y = numeric())
  
  # --- 4. K-Fold 主迴圈 ---
  for (k in 1:k_folds) {
    cat(sprintf("\n--- 正在執行第 %d / %d 折 ---\n", k, k_folds))
    
    test_indices <- folds[[k]]
    train_indices <- setdiff(1:nrow(X_raw_shuffled), test_indices)
    
    X_cv_train_raw <- X_raw_shuffled[train_indices, ]
    Y_cv_train_raw <- Y_raw_shuffled[train_indices, , drop = FALSE]
    X_cv_test_raw <- X_raw_shuffled[test_indices, ]
    Y_cv_test_raw <- Y_raw_shuffled[test_indices, , drop = FALSE]
    
    # 內部正規化
    X_cv_min <- apply(X_cv_train_raw, 2, min, na.rm = TRUE)
    X_cv_max <- apply(X_cv_train_raw, 2, max, na.rm = TRUE)
    scale_cv_factors <- X_cv_max - X_cv_min
    scale_cv_factors[scale_cv_factors == 0] <- 1
    
    Y_cv_mean <- mean(Y_cv_train_raw, na.rm = TRUE)
    Y_cv_sd <- sd(Y_cv_train_raw, na.rm = TRUE)
    
    if (is.na(Y_cv_sd) || Y_cv_sd < 1e-9) {
      cat(sprintf("警告！第 %d 折訓練集Y標準差為零，跳過此折。\n", k))
      next
    }
    
    X_cv_train_scaled <- as.matrix(scale(X_cv_train_raw, center = X_cv_min, scale = scale_cv_factors))
    Y_cv_train_scaled <- matrix((as.numeric(Y_cv_train_raw) - Y_cv_mean) / Y_cv_sd, ncol = 1)
    X_cv_test_scaled <- as.matrix(scale(X_cv_test_raw, center = X_cv_min, scale = scale_cv_factors))
    
    # 資料淨化
    dimnames(X_cv_train_scaled) <- NULL
    dimnames(Y_cv_train_scaled) <- NULL
    
    # --- 呼叫傳入的模型訓練函式 ---
    cv_model <- model_train_func(
      X_train = X_cv_train_scaled,
      Y_train_scaled = Y_cv_train_scaled,
      lower_b_log = lower_b_log,
      upper_b_log = upper_b_log,
      y_mean = Y_cv_mean, 
      y_sd = Y_cv_sd
    )
    
    # 預測
    predictions_k <- data.frame(
      real_Y = as.vector(Y_cv_test_raw),
      predicted_Y = rep(NA, nrow(X_cv_test_raw))
    )
    
    for (j in 1:nrow(X_cv_test_scaled)) {
      x_new_point <- X_cv_test_scaled[j, ]
      # 呼叫傳入的模型預測函式
      pred_result <- model_predict_func(as.numeric(x_new_point), cv_model)
      predictions_k$predicted_Y[j] <- pred_result$mean
    }
    
    all_predictions <- rbind(all_predictions, predictions_k)
    cat(sprintf("--- 第 %d 折完成 ---\n", k))
  }
  
  # --- 5. 結果匯總 ---
  print("\n--- 交叉驗證結果匯總 ---")
  
  rmse_overall <- sqrt(mean((all_predictions$real_Y - all_predictions$predicted_Y)^2, na.rm = TRUE))
  r_squared <- NA
  if (nrow(all_predictions) > 1 && sd(all_predictions$predicted_Y, na.rm=TRUE) > 0) {
    r_squared <- cor(all_predictions$real_Y, all_predictions$predicted_Y)^2
  }
  mae <- mean(abs(all_predictions$real_Y - all_predictions$predicted_Y), na.rm = TRUE)
  
  # 繪圖
  plot_object <- NULL
  if (nrow(all_predictions) > 1) {
    plot_object <- ggplot(all_predictions, aes(x = real_Y, y = predicted_Y)) +
      geom_point(alpha = 0.7, color = "blue") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      coord_fixed(ratio = 1) +
      labs(
        title = "K-Fold 交叉驗證結果",
        subtitle = paste("R² =", round(r_squared, 4), " | 全局 RMSE =", round(rmse_overall, 4)),
        x = "真實 Y 值", y = "模型預測 Y 值"
      ) +
      theme_bw()
    print(plot_object)
  }
  
  # 回傳結果
  results <- list(
    rmse = rmse_overall,
    r_squared = r_squared,
    mae = mae,
    predictions = all_predictions,
    plot = plot_object
  )
  
  return(results)
}

#-----------------------------------------------------------------------------
#要先執行 6D_AGP_V3_REAL_RUN
#  使用新框架重新驗證 AGP 模型
# =========================================================================

# 1. 確保所有 AGP 相關函式都已載入 (執行 model_main.R)
source("C:/Users/USER/Desktop/PYCProfessor/定期紀錄/6D_AGP_V3_real_run.R") 

# 2. 定義 AGP 的超參數邊界
agp_bounds <- list(
  lower = c(rep(1e-2, 6), rep(1e-2, 3)), 
  upper = c(rep(20, 6),    rep(20, 3))
)

# 3. 呼叫驗證函式，傳入完整混合數據和AGP模型函式
agp_cv_results <- validate_gp_model(
  X_raw = X_train_raw, 
  Y_raw = Y_train_raw,
  model_train_func = train_AGP_model, # <--- 訓練函式的名字
  model_predict_func = predict_AGP,   # <--- 預測函式的名字
  hyper_params_bounds = agp_bounds,
  k_folds = 5
)

# 4. 檢視結果
print("AGP 模型在混合數據上的最終交叉驗證效能：")
print(paste("R-squared:", round(agp_cv_results$r_squared, 4)))
print(paste("RMSE:", round(agp_cv_results$rmse, 4)))