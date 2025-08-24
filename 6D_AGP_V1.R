library(dplyr)
library(ggplot2)
library(globpso)
library(lhs) #載入 LHS 套件 By Foreester Section 1.4.3
library(readxl)
library(readr)

# --- Step1： 準備包含交互作用的模擬 6D 訓練資料 ---
true_function_agp_interaction <- function(x) {
  # x 是一個長度為 6 的向量
  # f1 部分的貢獻 (模擬第一階段)
  f1 <- 2*x[1] + 1.5*x[2] + 1*x[3]
  
  # f2 部分的貢獻 (模擬第二階段)
  f2 <- 1.2*x[4] + 0.8*x[5] + 0.5*x[6]
  
  # 交互作用部分的貢獻
  interaction_effect <- 1.5 * f1 * f2
  return(f1 + f2 + interaction_effect)
}
#係數本身沒有物理意義。目的是建立一個已知結構
#用這個函式來生成一組模擬數據，包含了第一階段、第二階段和交互作用這三種效應。
#希望看到模型能夠成功地從數據中學習到這三種效應的權重程度


# 使用 lhs 套件產生 40 個點的 6D 實驗設計
set.seed(456) 
n_points <- 40
X_train <- randomLHS(n_points, 6) # 點的範圍在 [0,1]  
#By Foreester 2.1.1, Page 34

# 計算 Y_train，並加上一點雜訊
Y_train <- apply(X_train, 1, true_function_agp_interaction) + rnorm(n_points, mean = 0, sd = 0.05)
Y_train <- as.matrix(Y_train)

print("--- 6D AGP (含交互作用) 模擬訓練資料已生成 ---")
head(X_train)
head(Y_train)

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
kernel_AGP_2stage_with_interaction <- function(x, x_prime, params) {
  
  # 1. 解包參數
  thetas <- params[1:6]             # 6 個 length-scale 相關參數
  sigma_f1_sq <- params[7]          # 權重 for 階段1 (獨立效應)
  sigma_f2_sq <- params[8]          # 權重 for 階段2 (獨立效應)
  sigma_interaction_sq <- params[9] # 權重 for 交互作用
  
  # 2. 計算基礎核 k1 和 k2
  dist_sq_f1 <- sum(thetas[1:3] * (x[1:3] - x_prime[1:3])^2)
  k1 <- exp(-dist_sq_f1)
  
  dist_sq_f2 <- sum(thetas[4:6] * (x[4:6] - x_prime[4:6])^2)
  k2 <- exp(-dist_sq_f2)
  
  # 3. 交互作用核
  k_interaction <- k1 * k2
  
  # 4. 最終核函數是三個加權部分的總和  By Duvenaud 
  #    k_total = σ₁²*k₁ + σ₂²*k₂ + σ_int²*(k₁*k₂)
  return(sigma_f1_sq * k1 + sigma_f2_sq * k2 + sigma_interaction_sq * k_interaction)
}

# --- B：計算AGP (含交互作用)負對數概似函數 ---
# 這個函數就是要交給 globpso 去優化的目標
# 它會根據傳入的 theta，計算出一個分數（-ln(L)），globpso 會想辦法讓這個分數變到最小
# 這個函數現在接收 9 個參數
calculate_AGP_interaction_neg_log_likelihood <- function(params, X_train, Y_train) {
  n <- nrow(X_train)
  one_vector <- matrix(1, nrow = n, ncol = 1)
  #為了 Eq2.30 2.31 做準備
  
  # 根據傳入的 9 個參數和 6D 訓練資料，建立 Ψ 矩陣 ByForrester Equation (2.21) (Page 51)
  #  X_train 是矩陣，所以用 X_train[i, ] 來取一整列
  # 計算第 i 個點和第 j 個點之間的關聯性分數，並填入 Psi 矩陣的對應位置
  Psi <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      Psi[i, j] <- kernel_AGP_2stage_with_interaction(X_train[i, ], X_train[j, ], params)
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
  mu_hat <- as.numeric((t(one_vector) %*% Psi_inv_y) / (t(one_vector) %*% Psi_inv_1))
  # 計算 σ_hat² (Eq. 2.31)，經過推導所以長這樣
  sigma_sq_hat <- as.numeric((t(Y_train - one_vector * mu_hat) %*% Psi_inv_y) / n)
  if (sigma_sq_hat <= 0) return(1e10)
  
  #Forrester Equation (2.32) (Page 55)
  ln_det_Psi <- 2 * sum(log(diag(U)))
  log_likelihood <- -n/2 * log(sigma_sq_hat) - 0.5 * ln_det_Psi
  
  return(-log_likelihood)
}


# --- Step3：AGP (含交互作用) 模型訓練函式 ---
train_AGP_interaction_model <- function(X_train, Y_train, lower_b, upper_b) {
  
  # --- 使用 globpso 尋找最佳的 9 個超參數 ---
  print("--- 正在使用 globpso 訓練交互作用 AGP 模型 ---")
  globpso_result <- globpso(
    objFunc = calculate_AGP_interaction_neg_log_likelihood, 
    lower = lower_b, #  9 維向量
    upper = upper_b, #  9 維向量
    X_train = X_train,
    Y_train = Y_train
  )
  
  # 變數改名
  best_params <- globpso_result$par
  print("globpso 優化器找到的最佳 9 個超參數為:")
  print(best_params)
  
  # 使用 best_params 計算並儲存所有固定不變的模型參數
  n <- nrow(X_train)
  one_vector <- matrix(1, nrow = n, ncol = 1)
  
  # 用新的核函數和最佳參數來建立最終的 Psi 矩陣
  Psi <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) { 
    for (j in 1:n) { 
      Psi[i, j] <- kernel_AGP_2stage_with_interaction(X_train[i, ], X_train[j, ], best_params) 
    } 
  }
  Psi <- Psi + diag(n) * 1e-8
  U <- chol(Psi)
  
  Psi_inv_y <- backsolve(U, forwardsolve(t(U), Y_train))
  Psi_inv_1 <- backsolve(U, forwardsolve(t(U), one_vector))
  mu_hat <- as.numeric((t(one_vector) %*% Psi_inv_y) / (t(one_vector) %*% Psi_inv_1))
  sigma_sq_hat <- as.numeric((t(Y_train - one_vector * mu_hat) %*% Psi_inv_y) / n)
  
  # 打包成新的模型物件
  model_object <- list(
    X_train = X_train, 
    Y_train = Y_train, 
    params = best_params, # <-- 改名
    mu_hat = mu_hat,
    sigma_sq_hat = sigma_sq_hat, 
    U = U
  )
  
  return(model_object)
}

#  ----- Step4 主流程 -----
# 1. 定義 9 維超參數的搜尋邊界 
# Forrester Page 55
# 前 6 個 a_j, 後 3 個 sigma_f^2
lower_bounds <- c(rep(1e-3, 6), rep(1e-2, 3)) 
upper_bounds <- c(rep(100, 6),   rep(20, 3))   

# 2. 訓練交互作用 AGP 模型
agp_model_interaction <- train_AGP_interaction_model(
  X_train, 
  Y_train, 
  lower_b = lower_bounds, 
  upper_b = upper_bounds
)

# 3. 查看訓練結果 Duvenaud Section 3
print("--- 交互作用 AGP 模型訓練完成 ---")
print("找到的最佳超參數 (6個 a_j + 3個 sigma_f^2):")
best_params <- agp_model_interaction$params # 從模型物件中取出
names(best_params) <- c("a1", "a2", "a3", "a4", "a5", "a6", 
                        "sigma_f1_sq", "sigma_f2_sq", "sigma_interact_sq")
print(round(best_params, 4))

# --- Step5：使用 AGP (含交互作用) 模型進行預測 ---
predict_with_AGP_interaction <- function(x_new, model) {
  # 從模型物件中 unpack 所有已經算好的參數
  X_train <- model$X_train
  Y_train <- model$Y_train
  params <- model$params # <-- 解包 9 個參數
  mu_hat <- model$mu_hat
  sigma_sq_hat <- model$sigma_sq_hat
  U <- model$U
  
  n <- nrow(X_train)
  one_vector <- matrix(1, nrow = n, ncol = 1)
  
  # --- 計算 psi_vector (新點 x_new 與所有訓練點的關聯性向量) ---
  # Forrester Equation (2.33) (Page 59)
  psi_vector <- matrix(0, nrow = n, ncol = 1)
  for (i in 1:n) { 
    psi_vector[i] <- kernel_AGP_2stage_with_interaction(X_train[i, ], x_new, params) 
  }
  
  # --- 預測平均值 y_hat 和變異數 s²(x)
  Psi_inv_diff <- backsolve(U, forwardsolve(t(U), (Y_train - one_vector * mu_hat)))
  y_hat <- mu_hat + t(psi_vector) %*% Psi_inv_diff
  
  # Forrester Equation (3.1) (Page 84)
  Psi_inv_psi <- backsolve(U, forwardsolve(t(U), psi_vector))
  s_sq <- sigma_sq_hat * (1 - t(psi_vector) %*% Psi_inv_psi)
  s_sq <- max(0, s_sq) 
  
  return(list(mean = as.numeric(y_hat), sd = sqrt(as.numeric(s_sq))))
}

# --- Step6 : AGP 模型預測與視覺化 ---
library(ggplot2)

# 建立一個 2D 的繪圖網格 (例如，觀察維度1 和 維度4 的影響)
grid_size <- 30
x1_grid <- seq(0, 1, length.out = grid_size)
x4_grid <- seq(0, 1, length.out = grid_size)
test_grid <- expand.grid(X1 = x1_grid, X4 = x4_grid)

# 固定其他變數為平均值 0.5 (x2, x3, x5, x6)
fixed_vars <- c(0.5, 0.5, 0.5, 0.5)

# 對網格上每個點進行預測
grid_predictions <- apply(test_grid, 1, function(row) {
  # 將變動的維度(x1, x4)和固定的維度(x2,x3,x5,x6)組合成 6D 向量
  x_new <- c(row[1], fixed_vars[1:2], row[2], fixed_vars[3:4])
  predict_with_AGP_interaction(x_new, agp_model_interaction)$mean # 只取預測平均值
})

# 將預測結果加到網格資料中
test_grid$Z <- grid_predictions

# 從模型中提取權重，用於標題
s1_sq <- round(agp_model_interaction$params[7], 2)
s2_sq <- round(agp_model_interaction$params[8], 2)
s_int_sq <- round(agp_model_interaction$params[9], 2)

ggplot(test_grid, aes(x = X1, y = X4, z = Z)) +
  geom_contour_filled(alpha = 0.8) +
  # 加上訓練點的位置 (只顯示對應的維度)
  geom_point(data = as.data.frame(X_train), aes(x=V1, y=V4, z=NULL), 
             color="white", size=2, shape = 21, fill = "black") +
  labs(
    title = "AGP 模型預測曲面 (含交互作用)",
    subtitle = paste("權重 σ²: f1 =", s1_sq, ", f2 =", s2_sq, ", f1*f2 =", s_int_sq),
    x = "維度 1 (模擬溫度1)", 
    y = "維度 4 (模擬溫度2)"
  ) +
  theme_minimal()

# --- Step 7 : 視覺化模型的不確定性 ---
# 1. 對網格上每個點，計算其預測的標準差
grid_uncertainty <- apply(test_grid, 1, function(row) {
  x_new <- c(row[1], fixed_vars[1:2], row[2], fixed_vars[3:4])

  predict_with_AGP_interaction(x_new, agp_model_interaction)$sd 
})

# 2. 將不確定性結果，作為一個新的欄位 SD 加到網格資料中
test_grid$SD <- grid_uncertainty

plot_sd <- ggplot(test_grid, aes(x = X1, y = X4, z = SD)) +
  geom_contour_filled(alpha = 0.8) +
  # 加上訓練點的位置
  geom_point(data = as.data.frame(X_train), aes(x=V1, y=V4, z=NULL), 
             color="white", size=2, shape = 21, fill = "black") +
  labs(
    title = "AGP 模型預測不確定性 (標準差)",
    subtitle = "顏色越黃代表模型在該區域越不確定",
    x = "維度 1 (模擬溫度1)", 
    y = "維度 4 (模擬溫度2)",
    fill = "標準差 (sd)"
  ) +
  theme_minimal()

print(plot_sd)

#-----------------------------------------------------------------------
#檢查用但 sd 結果應該不合理
print(paste("學到的全局過程變異數 (sigma_sq_hat):", agp_model_interaction$sigma_sq_hat)) 
print(paste("不確定性圖上的最大 sd 值為:", max(grid_uncertainty)))