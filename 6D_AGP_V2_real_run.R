library(dplyr)
library(ggplot2)
library(globpso)
library(lhs) #載入 LHS 套件 By Foreester Section 1.4.3
library(readxl)
library(readr)
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
  
  thetas <- params[1:6]
  sigma_f1_sq <- params[7]
  sigma_f2_sq <- params[8]
  sigma_interaction_sq <- params[9]
  
  # 1. Covariance 計算 
  dist_sq_f1 <- sum(thetas[1:3] * (x1[1:3] - x2[1:3])^2)
  k1 <- exp(-dist_sq_f1)
  covariance <- sigma_f1_sq * k1
  
  if (is_x1_full && is_x2_full) {
    dist_sq_f2 <- sum(thetas[4:6] * (x1[4:6] - x2[4:6])^2)
    k2 <- exp(-dist_sq_f2)
    k_interaction <- k1 * k2
    covariance <- covariance + sigma_f2_sq * k2 + sigma_interaction_sq * k_interaction
  }
  
  # 2. 分別計算 x1 和 x2 各自的總信號方差
  total_variance_x1 <- sigma_f1_sq
  if (is_x1_full) {
    total_variance_x1 <- total_variance_x1 + sigma_f2_sq + sigma_interaction_sq
  }
  
  total_variance_x2 <- sigma_f1_sq
  if (is_x2_full) {
    total_variance_x2 <- total_variance_x2 + sigma_f2_sq + sigma_interaction_sq
  }
  
  # 3. 使用兩個總方差的幾何平均值作為正規化分母
  #    sqrt(a*b) 是 a 和 b 的幾何平均值
  normalizing_factor <- sqrt(total_variance_x1 * total_variance_x2)
  
  # 4. 返回相關性
  return(covariance / (normalizing_factor + 1e-9))
}

# --- B：計算AGP (含交互作用)負對數概似函數 ---
# 這個函數就是要交給 globpso 去優化的目標
# 它會根據傳入的 theta，計算出一個分數（-ln(L)），globpso 會想辦法讓這個分數變到最小
# 這個函數現在接收 9 個參數
calculate_AGP_neg_log_likelihood_logscale <- function(params_log, X_train, Y_train) {
  
  # 將傳入的對數尺度參數，用 exp() 轉換回原始尺度
  params <- exp(params_log)
  
  n <- nrow(X_train)
  one_vector <- matrix(1, nrow = n, ncol = 1)
  #為了 Eq2.30 2.31 做準備
  
  # 根據傳入的 9 個參數和 6D 訓練資料，建立 Ψ 矩陣 ByForrester Equation (2.21) (Page 51)
  #  X_train 是矩陣，所以用 X_train[i, ] 來取一整列
  # 計算第 i 個點和第 j 個點之間的關聯性分數，並填入 Psi 矩陣的對應位置
  Psi <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      Psi[i, j] <- kernel_AGP_mixed(X_train[i, ], X_train[j, ], params)
    }
  }
  
  #加上 nugget 保持穩定，防止 Psi 矩陣因為訓練點太接近而變得接近 singular 
  #By Forrester Page 57 (ill conditioning)
  Psi <- Psi + diag(n) * 1e-8 
  
  # Cholesky 分解
  U <- tryCatch(chol(Psi), error = function(e) NULL)
  if (is.null(U)) { return(1e10) } # 如果參數組合不好給予懲罰
  
  Psi_inv_y <- backsolve(U, forwardsolve(t(U), Y_train))
  Psi_inv_1 <- backsolve(U, forwardsolve(t(U), one_vector))
  
  # 計算 μ_hat (Eq. 2.30)
  mu_hat <- (t(one_vector) %*% Psi_inv_y)[1, 1] / (t(one_vector) %*% Psi_inv_1)[1, 1]
  # 計算 σ_hat² (Eq. 2.31)，經過推導所以長這樣
  sigma_sq_hat <- as.numeric((t(Y_train - one_vector * mu_hat) %*% Psi_inv_y) / n)
  if (sigma_sq_hat <= 0) return(1e10)
  
  #Forrester Equation (2.32) (Page 55)
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
  
  # --- 使用 getPSOInfo() 創建詳細的優化器設定 ---
  # nSwarm: 粒子數量
  # maxIter: 最大迭代次
  pso_settings <- getPSOInfo(
    nSwarm = 50, 
    maxIter = 100, 
    psoType = "quantum" 
  )
  
  # --- 呼叫 globpso，傳入設定物件 ---
  globpso_result <- globpso(
    objFunc = objective_function_wrapper,
    lower = lower_b_log,
    upper = upper_b_log,
    PSO_INFO = pso_settings # <--- 透過 PSO_INFO 傳遞設定
  )
  # 將優化結果從對數尺度轉換回原始尺度 
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
  Psi <- Psi + diag(n) * 1e-8
  U <- chol(Psi)
  
  Psi_inv_y <- backsolve(U, forwardsolve(t(U), Y_train_scaled))
  Psi_inv_1 <- backsolve(U, forwardsolve(t(U), one_vector))
  mu_hat_scaled <- (t(one_vector) %*% Psi_inv_y)[1, 1] / (t(one_vector) %*% Psi_inv_1)[1, 1]
  sigma_sq_hat_scaled <- as.numeric((t(Y_train_scaled - one_vector * mu_hat_scaled) %*% Psi_inv_y) / n)
  
  # 打包成新的模型物件
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

#  ----- Step4 主流程 -----
# 1. 定義 9 維超參數的搜尋邊界 
# Forrester Page 55
#  ----- Step 4: 主流程 - 參數設定與模型訓練 -----
#    前 6 個是 thetas (或 length-scales), 後 3 個是 sigmas (權重)
lower_bounds <- c(rep(1e-2, 6), rep(1e-2, 3)) 
upper_bounds <- c(rep(100, 6),    rep(20, 3))  

# 2. 將邊界轉換到對數尺度
#    因為目標函數和優化器在對數尺度上運作
lower_bounds_log <- log(lower_bounds)
upper_bounds_log <- log(upper_bounds)

# 3. 訓練 AGP 模型
#    使用抽樣並正規化後的真實數據 (X_train, Y_train_scaled)
#    傳入對數尺度的邊界
#    將全局的 Y 正規化因子 (Y_mean, Y_sd) 傳入進模型
real_data_model <- train_AGP_model(
  X_train = as.matrix(X_train), # 確保傳入的是矩陣格式
  Y_train_scaled = Y_train_scaled,
  lower_b_log = lower_bounds_log, 
  upper_b_log = upper_bounds_log,
  y_mean = Y_mean, 
  y_sd = Y_sd
)

# 4. 查看訓練結果
print("--- 真實數據 AGP 模型訓練完成 ---")
print("找到的最佳超參數 (原始尺度):")

# 從模型物件中取出最終的最佳參數
best_params_real <- real_data_model$params 

# 為參數向量命名，方便閱讀
names(best_params_real) <- c("theta1", "theta2", "theta3", "theta4", "theta5", "theta6", 
                             "sigma_f1_sq", "sigma_f2_sq", "sigma_interact_sq")

print(round(best_params_real, 4))

# --- Step5：使用 AGP (含交互作用) 模型進行預測 ---
predict_AGP <- function(x_new, model) {
  # 從模型物件中 unpack 所有已經算好的參數
  X_train <- model$X_train
  Y_train <- model$Y_train_scaled
  params <- model$params 
  mu_hat_scaled <- model$mu_hat_scaled       # 估計出的全局均值 (在正規化尺度上)
  sigma_sq_hat_scaled <- model$sigma_sq_hat_scaled
  U <- model$U
  Y_mean <- model$Y_mean
  Y_sd <- model$Y_sd 
  
  n <- nrow(X_train)
  one_vector <- matrix(1, nrow = n, ncol = 1)
  
  # --- 計算 psi_vector (新點 x_new 與所有訓練點的關聯性向量) ---
  # Forrester Equation (2.33) (Page 59)
  psi_vector <- matrix(0, nrow = n, ncol = 1)
  for (i in 1:n) { 
    psi_vector[i] <- kernel_AGP_mixed(X_train[i, ], x_new, params) 
  }
  
  # --- 預測平均值 y_hat 和變異數 s²(x)
  Psi_inv_diff <- backsolve(U, forwardsolve(t(U), (Y_train_scaled - one_vector * mu_hat_scaled)))
  y_hat_scaled <- mu_hat_scaled + t(psi_vector) %*% Psi_inv_diff
  
  # Forrester Equation (3.1) (Page 84)
  Psi_inv_psi <- backsolve(U, forwardsolve(t(U), psi_vector))
  s_sq_scaled <- sigma_sq_hat_scaled * (1 - t(psi_vector) %*% Psi_inv_psi)
  s_sq_scaled <- max(0, s_sq_scaled) 
  # 將結果反正規化
  y_hat_real <- y_hat_scaled * Y_sd + Y_mean
  sd_real <- sqrt(s_sq_scaled) * Y_sd
  
  return(list(mean = as.numeric(y_hat_real), sd = as.numeric(sd_real)))
}

# --- Step 6 & 7: 最終結果的視覺化 ---

# 1. 建立繪圖網格
grid_size <- 50
x1_grid <- seq(0, 1, length.out = grid_size)
x4_grid <- seq(0, 1, length.out = grid_size)
test_grid <- expand.grid(X1 = x1_grid, X4 = x4_grid)

# 2. 設定固定變數的值 
#    這些值是在正規化的 [0, 1] 空間中的
fixed_vars <- c(0.5, 0.5, 0.5, 0.5)

# 3. 計算網格點的預測值 (平均值和不確定性)
#    對網格上的每個點，呼叫在 Step 5 中定義好的預測函式
grid_predictions <- apply(test_grid, 1, function(row) {
  x_new <- c(row[1], fixed_vars[1:2], row[2], fixed_vars[3:4])
  # 呼叫正確的函式名稱，並傳入正確的模型物件
  predict_AGP(x_new, real_data_model)$mean 
})
test_grid$Z <- grid_predictions

grid_uncertainty <- apply(test_grid, 1, function(row) {
  x_new <- c(row[1], fixed_vars[1:2], row[2], fixed_vars[3:4])
  predict_AGP(x_new, real_data_model)$sd
})
test_grid$SD <- grid_uncertainty


# 4. 繪製預測曲面圖
# 從模型中提取權重
s1_sq <- round(real_data_model$params[7], 2)
s2_sq <- round(real_data_model$params[8], 2)
s_int_sq <- round(real_data_model$params[9], 2)

# 確保繪圖時參考的訓練數據是正規化後的 X_train
plot_mean <- ggplot(test_grid, aes(x = X1, y = X4, z = Z)) +
  geom_contour_filled(alpha = 0.8) +
  geom_point(data = as.data.frame(X_train), aes(x = Stage1_Temp_C, y = Stage2_Temp_C, z = NULL), 
             color = "white", size = 2, shape = 21, fill = "black") +
  labs(
    title = "【真實數據】AGP 模型預測曲面",
    subtitle = paste("權重 σ²: f1 =", s1_sq, ", f2 =", s2_sq, ", f1*f2 =", s_int_sq),
    x = "維度 1 (第一階段溫度，正規化)", 
    y = "維度 4 (第二階段溫度，正規化)",
    fill = "有效深度"
  ) +
  theme_minimal()

print(plot_mean)

# 5. 繪製不確定性圖
plot_sd <- ggplot(test_grid, aes(x = X1, y = X4, z = SD)) +
  geom_contour_filled(alpha = 0.8) +
  geom_point(data = as.data.frame(X_train), aes(x = Stage1_Temp_C, y = Stage2_Temp_C, z = NULL), 
             color = "white", size = 2, shape = 21, fill = "black") +
  labs(
    title = "【真實數據】AGP 模型預測不確定性 (標準差)",
    subtitle = "顏色越黃代表模型在該區域越不確定",
    x = "維度 1 (第一階段溫度，正規化)", 
    y = "維度 4 (第二階段溫度，正規化)",
    fill = "標準差 (sd)"
  ) +
  theme_minimal()

print(plot_sd)


# 6. 檢查輸出
print(paste("學到的全局過程變異數 (正規化尺度):", round(real_data_model$sigma_sq_hat_scaled, 4))) 
print(paste("圖上的最大標準差 (真實尺度):", round(max(grid_uncertainty), 4)))

#--------------------------------------------------------------------
# --- Step 8: 驗證實驗繪製穿過某個特定真實訓練點的二維切片圖 ---

print("--- 執行高維切片驗證實驗 ---")

# 1. 從正規化後的訓練集中，挑選一個點作為錨點
#    可以選擇任意一個點，例如第 50 個
anchor_point_index <- 50 
anchor_point <- X_train[anchor_point_index, ]

print("選定的錨點座標為 (在正規化 [0,1] 空間中)：")
print(anchor_point)

# 2. 根據錨點，設定新的固定維度值
#    依然讓 x1 和 x4 變化，但將 x2, x3, x5, x6 固定在錨點的值上
fixed_vars_anchor <- c(
  anchor_point$Stage1_Time_hr,  # x2
  anchor_point$Stage1_Cs_pct,   # x3
  anchor_point$Stage2_Time_hr,  # x5
  anchor_point$Stage2_Cs_pct    # x6
)

# 3. 創建繪圖網格
grid_size_anchor <- 50
x1_grid_anchor <- seq(0, 1, length.out = grid_size_anchor)
x4_grid_anchor <- seq(0, 1, length.out = grid_size_anchor)
test_grid_anchor <- expand.grid(X1 = x1_grid_anchor, X4 = x4_grid_anchor)

# 4. 對這個新的切片平面上的每個點，重新計算其預測標準差
grid_uncertainty_anchor <- apply(test_grid_anchor, 1, function(row) {
  x_new <- c(row[1], fixed_vars_anchor[1:2], row[2], fixed_vars_anchor[3:4])
  
  predict_AGP(x_new, real_data_model)$sd 
})

# 5. 將不確定性結果添加到新的網格數據中
test_grid_anchor$SD <- grid_uncertainty_anchor

# 6. 繪製新的不確定性圖
plot_sd_anchor <- ggplot(test_grid_anchor, aes(x = X1, y = X4, z = SD)) +
  geom_contour_filled(alpha = 0.8) +
  # 疊加上所有訓練點的投影以便比較
  geom_point(data = as.data.frame(X_train), aes(x = Stage1_Temp_C, y = Stage2_Temp_C, z = NULL), 
             color="white", size=2, shape = 21, fill = "black") +
  # 特別標示錨點
  geom_point(
    data = data.frame(X1 = anchor_point$Stage1_Temp_C, X4 = anchor_point$Stage2_Temp_C),
    aes(x = X1, y = X4, z = NULL),
    color = "red", size = 5, shape = 3, stroke = 1.5
  ) +
  labs(
    title = paste("【驗證實驗】穿過訓練點 #", subset_indices[anchor_point_index], " 的不確定性切片圖", sep=""),
    subtitle = paste("固定維度 (x2,x3,x5,x6) =", paste(round(fixed_vars_anchor, 2), collapse=", ")),
    x = "維度 1 (第一階段溫度，正規化)", 
    y = "維度 4 (第二階段溫度，正規化)",
    fill = "標準差 (sd)"
  ) +
  theme_minimal()

print(plot_sd_anchor)