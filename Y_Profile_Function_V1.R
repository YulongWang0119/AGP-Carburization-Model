library(dplyr)
library(ggplot2)
library(globpso)
library(lhs) #載入 LHS 套件 By Foreester Section 1.4.3
library(readxl)
library(readr)
library(caret) # 載入 caret 套件

# --- Step2 :函式定義區 ---
# --- A: 實作高斯核函數 ---
# 這個函數會計算兩點 x1 和 x2 之間的關聯性 (對應教科書 Eq. 2.20)
# theta 是一個超參數，控制關聯性隨距離下降的速度
# theta 越大，關聯性下降得越快
#kernel_gaussian <- function(x1, x2, theta) {
#  return(exp(-theta * (x1 - x2)^2))
#}

# --- A: AGP 核函數 (包含一階交互作用) ---
# 這個核函數需要 9 個超參數
# Forrester Page 70 Equation (2.19) 
# --- 能處理混合維度的 AGP 核函數 ---
kernel_AGP_mixed <- function(x1, x2, params) {
  
  is_x1_full <- !any(is.na(x1[4:6]))
  is_x2_full <- !any(is.na(x2[4:6]))
  
  thetas_f1 <- params[1:3]
  thetas_f2 <- params[4:6]
  sigma_f1_sq <- params[7]
  sigma_f2_sq <- params[8]
  sigma_interaction_sq <- params[9]
  
  # 1. 計算第一階段核函數 (永遠都算)
  dist_sq_f1 <- sum(thetas_f1 * (x1[1:3] - x2[1:3])^2)
  k1 <- sigma_f1_sq * exp(-dist_sq_f1)
  
  # 2. 初始化第二階段與交互作用核函數為 0
  k2 <- 0
  k_interaction <- 0
  
  # 3. 只有當兩個點都是六維時，才計算它們
  if (is_x1_full && is_x2_full) {
    dist_sq_f2 <- sum(thetas_f2 * (x1[4:6] - x2[4:6])^2)
    k2 <- sigma_f2_sq * exp(-dist_sq_f2)
    
    # 交互作用核函數可以是兩個基礎核的乘積
    k_interaction <- sigma_interaction_sq * exp(-dist_sq_f1) * exp(-dist_sq_f2)
  }
  
  # 4. 將所有部分加總
  return(k1 + k2 + k_interaction)
}

# --- B：計算AGP (含交互作用)負對數概似函數 ---
# 這個函數就是要交給 globpso 去優化的目標
# 它會根據傳入的 theta，計算出一個分數（-ln(L)），globpso 會想辦法讓這個分數變到最小
# 這個函數現在接收 9 個參數
# =========================================================================
# 【最終版本】calculate_AGP_neg_log_likelihood_logscale 函式
# =========================================================================
# --- B：計算AGP (含交互作用)負對數概似函數 ---
calculate_AGP_neg_log_likelihood_logscale <- function(params_log, X_train, Y_train) {
  
  # 將傳入的對數尺度參數，用 exp() 轉換回原始尺度
  params <- exp(params_log)
  
  n <- nrow(X_train)
  one_vector <- matrix(1, nrow = n, ncol = 1)
  
  # 根據傳入的 9 個參數和訓練資料，建立 Psi 矩陣
  Psi <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      Psi[i, j] <- kernel_AGP_mixed(X_train[i, ], X_train[j, ], params)
    }
  }
  
  # 加上 nugget 保持穩定
  Psi <- Psi + diag(n) * 1e-6
  
  # Cholesky 分解
  U <- tryCatch(chol(Psi), error = function(e) NULL)
  if (is.null(U)) { 
    return(1e10) # 如果參數組合不好給予懲罰
  }
  
  Psi_inv_y <- backsolve(U, forwardsolve(t(U), Y_train))
  Psi_inv_1 <- backsolve(U, forwardsolve(t(U), one_vector))
  
  # 計算 μ_hat (Eq. 2.30)
  mu_hat <- (t(one_vector) %*% Psi_inv_y)[1, 1] / (t(one_vector) %*% Psi_inv_1)[1, 1]
  
  # 計算 σ_hat² (Eq. 2.31)
  sigma_sq_hat <- as.numeric((t(Y_train - one_vector * mu_hat) %*% Psi_inv_y) / n)
  
  if (sigma_sq_hat <= 0) {
    return(1e10)
  }
  
  # Forrester Equation (2.32) (Page 55)
  ln_det_Psi <- 2 * sum(log(diag(U)))
  log_likelihood <- -n/2 * log(sigma_sq_hat) - 0.5 * ln_det_Psi
  
  return(-log_likelihood)
}


# --- Step3：AGP (含交互作用) 模型訓練函式 ---
train_AGP_model <- function(X_train, Y_train_scaled, lower_b_log, upper_b_log, y_mean, y_sd) {
  
  objective_function_wrapper <- function(p_log) {
    return(calculate_AGP_neg_log_likelihood_logscale(
      params_log = p_log, 
      X_train = X_train, 
      Y_train = Y_train_scaled
    ))
  }
  
  print("--- 正在使用 globpso 訓練交互作用 AGP 模型 ---")
  
  pso_settings <- getPSOInfo(
    nSwarm = 50, 
    maxIter = 100, 
    psoType = "quantum"
  )
  
  globpso_result <- globpso(
    objFunc = objective_function_wrapper,
    lower = lower_b_log,
    upper = upper_b_log,
    PSO_INFO = pso_settings
  )
  
  best_params_log <- globpso_result$par
  best_params <- exp(best_params_log) 
  
  print("globpso 優化器找到的最佳 9 個超參數為:")
  print(best_params)
  
  # 使用 best_params 計算並儲存所有固定不變的模型參數
  n <- nrow(X_train)
  one_vector <- matrix(1, nrow = n, ncol = 1)
  
  # 用新的核函數和最佳參數來建立最終的 Psi 矩陣
  Psi <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) { 
    for (j in 1:n) { 
      Psi[i, j] <- kernel_AGP_mixed(X_train[i, ], X_train[j, ], best_params) 
    } 
  }
  Psi <- Psi + diag(n) * 1e-6
  U <- chol(Psi)
  
  Psi_inv_y <- backsolve(U, forwardsolve(t(U), Y_train_scaled))
  Psi_inv_1 <- backsolve(U, forwardsolve(t(U), one_vector))
  mu_hat_scaled <- (t(one_vector) %*% Psi_inv_y)[1, 1] / (t(one_vector) %*% Psi_inv_1)[1, 1]
  
  sigma_sq_hat_scaled <- as.numeric((t(Y_train_scaled - one_vector * mu_hat_scaled) %*% Psi_inv_y) / n)
  
  model_object <- list(
    X_train = X_train, 
    Y_train_scaled = Y_train_scaled, 
    params = best_params,
    mu_hat_scaled = mu_hat_scaled,
    sigma_sq_hat_scaled = sigma_sq_hat_scaled, 
    U = U,
    Y_mean = y_mean, 
    Y_sd = y_sd
  )
  
  return(model_object)
}


# --- Step5：使用 AGP (含交互作用) 模型進行預測 ---
predict_AGP <- function(x_new, model) {
  X_train <- model$X_train
  Y_train_scaled <- model$Y_train_scaled
  params <- model$params 
  mu_hat_scaled <- model$mu_hat_scaled
  sigma_sq_hat_scaled <- model$sigma_sq_hat_scaled
  U <- model$U
  Y_mean <- model$Y_mean
  Y_sd <- model$Y_sd 
  
  n <- nrow(X_train)
  one_vector <- matrix(1, nrow = n, ncol = 1)
  
  psi_vector <- matrix(0, nrow = n, ncol = 1)
  for (i in 1:n) { 
    psi_vector[i] <- kernel_AGP_mixed(X_train[i, ], x_new, params) 
  }
  
  # --- 預測平均值 y_hat  ---
  diff_for_pred <- as.vector(Y_train_scaled) - as.vector(one_vector * mu_hat_scaled)
  Psi_inv_diff <- backsolve(U, forwardsolve(t(U), matrix(diff_for_pred, ncol = 1)))
  y_hat_scaled <- mu_hat_scaled + t(psi_vector) %*% Psi_inv_diff
  
  # --- 預測變異數 s²(x)  ---
  # Forrester Equation (3.1) (Page 84)
  Psi_inv_psi <- backsolve(U, forwardsolve(t(U), psi_vector))
  
  # 需要知道核函數在 x_new 點自身的共變異數 k(x_new, x_new)
  # 對於 AGP 核，當 x1=x2 時，距離為0，exp()項為1。
  # 需要計算 x_new 是三維還是六維
  is_x_new_full <- !any(is.na(x_new[4:6]))
  
  # 計算 k(x_new, x_new)
  k_x_new_x_new <- params[7] # sigma_f1_sq
  if (is_x_new_full) {
    k_x_new_x_new <- k_x_new_x_new + params[8] + params[9]
  }
  
  # 穩健的變異數公式
  s_sq_scaled <- sigma_sq_hat_scaled * (k_x_new_x_new - t(psi_vector) %*% Psi_inv_psi)
  
  s_sq_scaled <- max(0, s_sq_scaled) 
  
  # --- 反正規化  ---
  y_hat_real <- y_hat_scaled * Y_sd + Y_mean
  sd_real <- sqrt(s_sq_scaled) * Y_sd
  
  return(list(mean = as.numeric(y_hat_real), sd = as.numeric(sd_real)))
}

#----------------------------------------------------------------------------------
#驗證函數
# =========================================================================
#             validate_gp_model_strict (V2 - 內建完整繪圖功能)
# =========================================================================
# 這是一個全新的函式，用於執行最嚴謹的交叉驗證

validate_gp_model <- function(X_raw, Y_matrix_raw, j_weight_to_validate,
                              model_train_func, model_predict_func,
                              hyper_params_bounds, k_folds = 5, seed = 123,
                              strata_labels = NULL, decompose_func) {
  
  cat(sprintf("\n--- 開始對【權重 %d】進行最嚴謹的 %d-Fold 交叉驗證 ---\n", j_weight_to_validate, k_folds))
  
  set.seed(seed)
  lower_b_log <- log(hyper_params_bounds$lower)
  upper_b_log <- log(hyper_params_bounds$upper)
  
  # --- 資料預處理與切分 (同時切分 X 和 Y_matrix) ---
  shuffle_indices <- sample(1:nrow(X_raw))
  X_raw_shuffled <- X_raw[shuffle_indices, ]
  Y_matrix_shuffled <- Y_matrix_raw[, shuffle_indices]
  #打亂數據：在進行 K-摺分割前，先將數據的原始順序完全打亂。可以避免資料原始的排列順序對分割結果造成偏差。
  #sample(1:nrow(X_raw))：產生一個從 1 到總樣本數的隨機排列組合。
  #X_raw[shuffle_indices, ]：按照這個隨機排列重新排序 X 矩陣的列。
  
  if (!is.null(strata_labels)) {
    strata_labels_shuffled <- strata_labels[shuffle_indices]
    folds <- createFolds(strata_labels_shuffled, k = k_folds, list = TRUE, returnTrain = FALSE)
  } else {
    folds <- createFolds(1:nrow(X_raw_shuffled), k = k_folds, list = TRUE, returnTrain = FALSE)
  }
#提供了 strata_labels（分層標籤），createFolds 會進行分層抽樣，確保每一折中單階段和雙階段數據的比例大致相同。
#對於處理不平衡數據集重要
#folds：是一個 list，包含了 k_folds 個元素。folds[[1]] 存的是第一折的樣本索引，folds[[2]] 是第二折的，依此類推。
  
  all_predictions <- data.frame()
  #建立一個空的 data.frame。用來逐一收集並儲存每一折產生的預測結果。
  
  # --- K-Fold 主迴圈 ---
  for (k in 1:k_folds) {
    cat(sprintf("\n--- 正在執行第 %d / %d 折 ---\n", k, k_folds))
    
    test_indices <- folds[[k]]
    train_indices <- setdiff(1:nrow(X_raw_shuffled), test_indices)
    #for (k in 1:k_folds)：開始一個迴圈，從 k=1 執行到 k=5。
    #test_indices <- folds[[k]]：在第 k 次迴圈中，取出第 k 折的索引作為當次的測試集索引。
    #setdiff(...)：從全部索引中移除測試集索引，剩下的就是當次的訓練集索引。
    
    X_cv_train_raw <- X_raw_shuffled[train_indices, ]
    Y_matrix_cv_train <- Y_matrix_shuffled[, train_indices]
    X_cv_test_raw <- X_raw_shuffled[test_indices, ]
    Y_matrix_cv_test <- Y_matrix_shuffled[, test_indices]
    #根據這兩組索引，從打亂後的數據中分割出當次迴圈專用的訓練集 (_cv_train) 和測試集 (_cv_test)
    
    cat("  > 正在對訓練集進行內部 KL 降維...\n")
    decomp_train <- decompose_func(Y_matrix_cv_train, variance_threshold = 0.999)
    #在迴圈內部，只對當次的訓練集 Y_matrix_cv_train (例如 160 筆) 進行 KL 降維。
    #這意味著每一折的「模型基礎」（平均曲線、基函數）都是獨立計算的，完全沒有用到當次測試集的任何資訊。
    # # ==========================================================
    # #           【新增的繪圖程式碼，從這裡開始】
    # # ==========================================================
    # 
    # # 建立繪圖用的資料框
    # depth_vector_cv <- first_profile$depth_mm # 假設所有 profile 深度點都一樣
    # plot_df_base_cv <- data.frame(depth = depth_vector_cv, mean_curve = decomp_train$Y_mean_curve)
    # 
    # # 準備模式數據
    # num_modes_cv <- decomp_train$num_components
    # plot_df_wide_cv <- plot_df_base_cv %>%
    #   bind_cols(as.data.frame(decomp_train$eigenvectors))
    # colnames(plot_df_wide_cv)[-c(1,2)] <- paste0("模式", 1:num_modes_cv)
    # 
    # plot_df_tidy_cv <- plot_df_wide_cv %>%
    #   tidyr::pivot_longer(
    #     cols = starts_with("模式"),
    #     names_to = "模式",
    #     values_to = "變化量"
    #   )
    # 
    # # 繪圖
    # plot_internal_modes <- ggplot() +
    #   geom_line(data = plot_df_base_cv, aes(x = depth, y = mean_curve, color = "平均曲線"), linewidth = 1.2) +
    #   geom_line(data = plot_df_tidy_cv, aes(x = depth, y = 變化量, color = 模式), linetype = "dashed") +
    #   labs(
    #     title = paste("第", k, "折 (Fold) 的內部變化模式"),
    #     subtitle = paste("基於", ncol(Y_matrix_cv_train), "筆訓練數據"),
    #     x = "深度 (mm)", y = "濃度 / 變化量"
    #   ) +
    #   theme_bw(base_size = 12)
    # 
    # # 印出這張圖
    # print(plot_internal_modes)
    # 
    # # ==========================================================
    # #           【新增的繪圖程式碼，到這裡結束】
    # # ==========================================================
    
    if(j_weight_to_validate > decomp_train$num_components){
      cat(sprintf("  > 警告：在第 %d 折中，主成分數量 (%d) 不足 %d，跳過此折。\n", k, decomp_train$num_components, j_weight_to_validate))
      next
    }
    #一個安全檢查。如果某一折的訓練數據剛好比較單純，降維後的主成分數量 J 可能會比較少。
    #如果 J 小於想驗證的權重編號（例如 J=3 但想驗證 權重4），那就沒有對應的資料。
    #next 指令會直接跳過這次迴圈的剩餘部分，進入下一折。
    
    Y_cv_train_raw <- as.matrix(t(decomp_train$weights[j_weight_to_validate, , drop = FALSE]))
    #從當次降維的結果 decomp_train$weights 中，提取出第 j_weight_to_validate 行，這就是本次迴圈要訓練權重係數
    
    # --- 模型訓練 ---
    X_cv_min <- apply(X_cv_train_raw, 2, min, na.rm = TRUE); X_cv_max <- apply(X_cv_train_raw, 2, max, na.rm = TRUE)
    scale_cv_factors <- X_cv_max - X_cv_min; scale_cv_factors[scale_cv_factors[scale_cv_factors == 0 | is.na(scale_cv_factors)]] <- 1
    Y_cv_mean <- mean(Y_cv_train_raw, na.rm = TRUE); Y_cv_sd <- sd(Y_cv_train_raw, na.rm = TRUE)
    if (is.na(Y_cv_sd) || Y_cv_sd < 1e-9) { next }
    X_cv_train_scaled <- as.matrix(scale(X_cv_train_raw, center = X_cv_min, scale = scale_cv_factors))
    Y_cv_train_scaled <- matrix((Y_cv_train_raw - Y_cv_mean) / Y_cv_sd, ncol = 1)
    
    cv_model <- model_train_func(X_train = X_cv_train_scaled, Y_train_scaled = Y_cv_train_scaled,
                                 lower_b_log = lower_b_log, upper_b_log = upper_b_log,
                                 y_mean = Y_cv_mean, y_sd = Y_cv_sd)
    #只用當次的訓練數據 X_cv_train_raw 和 Y_cv_train_raw 來計算正規化參數，並訓練出一個臨時的、只在本次迴圈有效的模型 cv_model
    # --- 標準答案 ---
    Y_centered_test <- Y_matrix_cv_test - decomp_train$Y_mean_curve
    Lambda_inv_sqrt <- diag(decomp_train$eigenvalues^(-1/2), nrow = decomp_train$num_components)
    Phi_t <- t(decomp_train$eigenvectors)
    Y_weights_cv_test_all_real <- Lambda_inv_sqrt %*% Phi_t %*% Y_centered_test
    
    real_Y <- as.vector(Y_weights_cv_test_all_real[j_weight_to_validate, ])
    #計算標準答案。
    #它使用當次訓練集 (_cv_train) 降維得到的基準 (decomp_train$Y_mean_curve 等)，去對當次的測試集剖面 Y_matrix_cv_test 進行投影，得到真實的權重係數 real_Y
    #文獻3的16
    
    X_cv_test_scaled <- as.matrix(scale(X_cv_test_raw, center = X_cv_min, scale = scale_cv_factors))
    predicted_Y <- rep(NA, length(real_Y))
    for (i in 1:length(real_Y)) {
      predicted_Y[i] <- model_predict_func(X_cv_test_scaled[i,], cv_model)$mean
    }
    #計算模型的預測。
    #使用當次訓練集計算出的正規化參數 (X_cv_min等)，去正規化當次的測試集輸入 X_cv_test_raw。
    #然後用當次訓練出的 cv_model，對正規化後的測試集 X_cv_test_scaled 逐一進行預測，得到預測的權重係數 predicted_Y
    
    source_labels_this_fold <- if(!is.null(strata_labels)) strata_labels_shuffled[test_indices] else NA
    all_predictions <- rbind(all_predictions, 
                             data.frame(real_Y = real_Y, 
                                        predicted_Y = predicted_Y, 
                                        source = source_labels_this_fold))
  }
  #將這一折得到的 real_Y 和 predicted_Y 配對起來，存入 all_predictions 這個總表中
  
  # --- 結果匯總與最終繪圖 ---
  rmse <- sqrt(mean((all_predictions$real_Y - all_predictions$predicted_Y)^2, na.rm = TRUE))
  r_squared <- cor(all_predictions$real_Y, all_predictions$predicted_Y, use = "pairwise.complete.obs")^2
  #在所有 k_folds 次迴圈都結束後，all_predictions 中已經匯集了 200 筆來自不同折的預測結果。
  #函式對這全部 200 筆結果，計算總體的 RMSE 和 R²。
  #繪製一張包含全部 200 個點的性能圖，並將最終的性能指標和圖表物件打包成 list 回傳。
  
  plot_object <- ggplot(all_predictions, aes(x = real_Y, y = predicted_Y)) + 
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    coord_fixed(ratio = 1) +
    labs(
      title = paste("最嚴謹交叉驗證結果：權重", j_weight_to_validate),
      subtitle = paste("R² =", round(r_squared, 4), "| RMSE =", round(rmse, 4)),
      x = "真實權重值", y = "模型預測權重值"
    ) +
    theme_bw(base_size = 14)
  
  if (!is.null(strata_labels)) {
    plot_object <- plot_object +
      geom_point(aes(color = source, shape = source), alpha = 0.8, size = 3) + 
      scale_color_manual(name = "數據來源", values = c("2-stage" = "blue", "1-stage" = "green")) +
      scale_shape_manual(name = "數據來源", values = c("2-stage" = 16, "1-stage" = 17)) +
      theme(legend.position = "bottom")
  } else {
    plot_object <- plot_object + geom_point(alpha = 0.7, color = "blue", size = 2.5)
  }
  
  print(plot_object)
  
  return(list(rmse = rmse, r_squared = r_squared, predictions = all_predictions, plot = plot_object))
}

#---------------------------------------------------------------------------------------------
#呼叫這個函數，把測試集的 Y_test 傳進去，它會用跟訓練時完全相同的基準（Y_mean_curve, eigenvectors, eigenvalues 
#都是從訓練集算出來的），將 Y_test 投影到特徵空間，計算出對應的真實權重，稱為 weights_true。
#公式16
project_Y_to_weights <- function(Y_matrix_new, Y_mean_curve, eigenvectors, eigenvalues) {
  # 確保與 decompose_Y 的權重計算方式完全一致
  Y_centered_new <- sweep(Y_matrix_new, 1, Y_mean_curve, "-")
  
  # 提取需要的部分
  J <- ncol(eigenvectors)
  Lambda_J_inv_sqrt <- diag(eigenvalues[1:J]^(-1/2), nrow = J, ncol = J)
  Phi_J <- eigenvectors[, 1:J]
  
  # 進行與 decompose_Y 相同的矩陣運算
  weights_matrix <- Lambda_J_inv_sqrt %*% t(Phi_J) %*% Y_centered_new
  #t(Phi_J)：將 N x J 的基函數矩陣 Φ̂_J 轉置為 J x N。
  #%*% Y_centered_new：用 J x N 的 t(Phi_J) 左乘 N x n_new 的 Y_centered_new。
  #這一步是「投影」。它將每一條新的、中心化後的曲線（Y_centered_new 的每一欄）投影到 J 個基函數（t(Phi_J) 的每一行）上。
  #結果是一個 J x n_new 的矩陣。第 (j, i) 個元素代表第 i 條測試曲線在第 j 個基函數上的投影值。
  #Lambda_J_inv_sqrt %*% ...：用 J x J 的對角矩陣 Λ̂_J^(-1/2) 左乘上一步的投影結果。
  #這一步是「標準化它將投影矩陣的每一行除以對應特徵值的平方根，使得計算出的權重係數在尺度上與訓練時的權重係數保持一致。
  
  return(t(weights_matrix)) # 回傳 n x J 的矩陣
}
  #t(weights_matrix)：將 J x n_new 的矩陣轉置成 n_new x J。

# =========================================================================
#：KL 降維函式

decompose_Y <- function(Y_matrix, variance_threshold = 0.999) {
  #Y_matrix：原始數據一個矩陣
  #variance_threshold：用來決定要保留多少主成分的閾值，預設值為 0.999（即保留能解釋 99.9% 變異數的主成分）
  
  # N 是深度點數，n 是實驗筆數，目前切成101個網格8mm，訓練集200筆混和數據
  N <- nrow(Y_matrix)
  n <- ncol(Y_matrix)
  
  # 步驟 1: 計算平均曲線 (論文中的 μ̂_M)
  #rowMeans 計算出在每個空間點上的平均值， Y_mean_curve 是一個長度為 N 的向量，N 是深度
  Y_mean_curve <- rowMeans(Y_matrix)
  
  # 步驟 2: 中心化數據 (讓每條曲線減去平均曲線)
  #讓 Y_matrix 的每一欄都減去它。實現了數據中心化
  #中心化的目的是移除數據的整體趨勢，專注於分析數據的變化模」。KL 展開理論是建立在中心化的隨機過程上的
  Y_centered <- Y_matrix - Y_mean_curve
  
  # 步驟 3: 計算樣本共變異數矩陣 S (論文公式 14)
  # 注意 R 的 cov() 函式是計算欄位間的共變異數，我要的是列之間的。
  # 所以需要先轉置 t()，讓每一列變成一個觀測。
  # S 矩陣的維度會是 N x N。
  # cov(...)：計算這個 n x N 矩陣的共變異數。它會計算 N 個欄位（空間點）之間的兩兩共變異關係，因此輸出的 S 是一個 N x N 的矩陣
  S <- cov(t(Y_centered))
  #S <- (1 / n) * (Y_centered %*% t(Y_centered))
  #公式的意義：它在計算 N 個空間點之間的樣本共變異數。
  #矩陣 S 的第 (j,k) 個元素，就是 (1/n) * Σ [(Y_i[j] - Ȳ[j]) * (Y_i[k] - Ȳ[k])
  
  # 步驟 4: 對 S 進行特徵分解 (主成分分析) (論文公式 15)
  eigen_decomp <- eigen(S)
  eigenvalues <- eigen_decomp$values      # 特徵值 (λ̂)
  eigenvectors <- eigen_decomp$vectors   # 特徵向量 (Φ̂，即主要變化模式)
  #公式 (15) S = Φ̂Λ̂Φ̂ᵀ。
  #eigenvectors 就是 Φ̂，它的每一欄 Φ̂_i 就是一個基函數 φ_i(s) 的離散版本，代表一種主要的變化模式
  #eigenvalues 就是對角矩陣 Λ̂ 的對角線元素 λ̂_i，代表了每種模式的重要性（方差）
  
  # 步驟 5: 決定要保留多少個主成分 J
  cumulative_variance <- cumsum(eigenvalues) / sum(eigenvalues)
  J <- min(which(cumulative_variance >= variance_threshold))
  cat(sprintf("使用 %d 個主成分來解釋 %.1f%% 的變異。\n", J, variance_threshold * 100))
  
  # 步驟 6: 計算權重係數矩陣 Xi_hat (論文公式 16)
  # 公式: Ξ̂_i = λ̂_i^(-1/2) * Φ̂_i^T * (Y - Ȳ)
  # 一次計算所有 J 個權重
  
  # 提取需要的部分
  Lambda_J_inv_sqrt <- diag(eigenvalues[1:J]^(-1/2))
  Phi_J <- eigenvectors[, 1:J]
  #diag 的意思 在左邊乘以一個對角矩陣
  #Phi_J 是 N x J，所以 t(Phi_J) 是 J x N。Y_centered 是 N x n
  #投影 t(Phi_J) %*% Y_centered
  #這個矩陣的每一列，對應一個基函數；每一欄，對應一次實驗。
  #第 (i, j) 個元素代表第 j 次實驗在 φ_i 這個基函數上的投影值（還未縮放）。
  
    # 進行矩陣運算
  weights_matrix <- Lambda_J_inv_sqrt %*% t(Phi_J) %*% Y_centered
  #weights_matrix 就是 Ξ̂*，它的維度是 J x n。將原始的 N x n 函數型數據，壓縮成了 J x n 的係數數據
  #J個主成分 * n筆訓練集數據
  #接下來對這個 weights_matrix 的每一列建立一個 GP 模型（x -> ξ_j(x)）。
  
  # 步驟 7: 回傳所有重要結果
  results <- list(
    Y_mean_curve = Y_mean_curve,
    eigenvalues = eigenvalues[1:J],
    eigenvectors = Phi_J,
    weights = weights_matrix,
    num_components = J
  )
  
  return(results)
}



# =========================================================================
#  預測與重構 

# --- 4a. 創建一個新的預測與重構函式 ---
predict_and_reconstruct_profile <- function(x_new_raw, model_list, decomposition_results, X_min, X_scale_factors) {
  #model_list: 一個列表，包含了為 J 個係數 ξ 訓練好的 J 個 GP 模型。
  #ecomposition_results: 一個列表，包含了 decompose_Y 函數計算出的所有結果 (μ_M, Φ_J, λ_J, J 等)。
  #X_min, X_scale_factors: 這是在訓練模型前，對輸入 X_train 進行標準化/正規化時所使用的均值和標準差
  
  # 1. 正規化新的輸入點 x_new
  x_new_scaled <- as.numeric(scale(matrix(x_new_raw, nrow=1), center = X_min, scale = X_scale_factors))
  #GP 模型是在正規化後的輸入數據 X_train_scaled 上訓練的，所以任何新的輸入點 x_new 在進行預測前，
  #必須使用與訓練數據完全相同的參數（X_min 和 X_scale_factors）進行正規化
  
  J <- decomposition_results$num_components
  predicted_weights <- numeric(J)
  #從 decomposition_results 中取出主成分的數量 J
  #創建一個長度為 J 的空向量 predicted_weights 來儲存預測結果
  
  # 2. 用 J 個模型分別預測 J 個權重
  for (j in 1:J) {
    prediction <- predict_AGP(x_new_scaled, model_list[[j]])
    predicted_weights[j] <- prediction$mean
  }
  #在預測文獻三中的係數向量 Ξ̂(x_new)
  #predicted_weights 就是 (ξ̂₁(x_new), ξ̂₂(x_new), ..., ξ̂_J(x_new)) 這個向量
  
  # 3. 線性重構回完整的 Profile 曲線
  Y_mean_curve <- decomposition_results$Y_mean_curve
  eigenvectors <- decomposition_results$eigenvectors
  eigenvalues <- decomposition_results$eigenvalues # 論文公式17，重構時可能需要乘回 sqrt(lambda)
  # # Tan (2017) 的重構公式更接近 Y_hat = Y_mean + Phi * Lambda^(1/2) * Xi_hat
  # reconstructed_profile <- Y_mean_curve + eigenvectors %*% (sqrt(eigenvalues) * predicted_weights)
  # return(as.vector(reconstructed_profile))
  
  # 【新增】建立一個矩陣來儲存每個模式的貢獻曲線
  # 維度是 N x J (深度點數 x 模式數)
  mode_contributions <- matrix(0, nrow = length(Y_mean_curve), ncol = J)
  
  # 計算加權後的係數向量
  weighted_coeffs <- sqrt(eigenvalues) * predicted_weights
  
  # 計算每個模式的貢獻
  for (j in 1:J) {
    # 第 j 個模式的貢獻 = Φ_j * (sqrt(λ_j) * ξ_j)
    mode_contributions[, j] <- eigenvectors[, j] * weighted_coeffs[j]
  }
  
  # 重構最終曲線 (現在是平均曲線 + 所有模式貢獻的總和)
  reconstructed_profile <- Y_mean_curve + rowSums(mode_contributions)
  
  # --- 4. 回傳一個包含所有資訊的列表 (【修改處】) ---
  return(
    list(
      final_prediction = as.vector(reconstructed_profile),
      predicted_weights = predicted_weights,
      mean_curve = Y_mean_curve,
      mode_contributions = as.data.frame(mode_contributions)
    )
  )
}