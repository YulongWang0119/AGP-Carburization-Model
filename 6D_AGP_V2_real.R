#install.packages("FNN")
#install.packages("randomForest") 
library(randomForest)

recipes_file_path <- "C:/Users/USER/Desktop/PYCProfessor/CarburizationDB_Full/recipes.csv" 
all_recipes <- read.csv(recipes_file_path)

X_real_raw <- all_recipes[, c("Stage1_Temp_C", "Stage1_Time_hr", "Stage1_Cs_pct", 
                              "Stage2_Temp_C", "Stage2_Time_hr", "Stage2_Cs_pct")]

# 將其轉換為數值矩陣
X_real_raw <- as.matrix(X_real_raw)

# 查看數據維度是 3600 x 6
print(dim(X_real_raw))

#-------------------------------------------------------------------------

# 函式：從一個 profile 檔案路徑計算有效滲碳深度
calculate_effective_depth <- function(profile_file_path, target_concentration = 0.2) {
  profile_data <- tryCatch({
    read.csv(profile_file_path)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(profile_data)) {
    return(NA) # 如果檔案讀取失敗，返回 NA
  }
  
  # 確保欄位名稱正確
  colnames(profile_data) <- c("depth_mm", "concentration_pct")
  
  # 找到第一個濃度 <= target_concentration 的位置
  # which() 會返回所有滿足條件的索引
  indices_below_target <- which(profile_data$concentration_pct <= target_concentration)
  
  if (length(indices_below_target) == 0) {
    # 如果所有點的濃度都高於目標值，表示滲碳非常深，可能返回最大深度或 NA
    return(max(profile_data$depth_mm)) 
  }
  
  # 第一個低於目標的點的索引
  first_index_below <- min(indices_below_target)
  
  # 如果第一個點就低於目標，深度為 0
  if (first_index_below == 1) {
    return(0)
  }
  
  # 進行線性內插，找到精確的深度
  # p1 是目標點的前一個點，p2 是目標點本身
  p1 <- profile_data[first_index_below - 1, ]
  p2 <- profile_data[first_index_below, ]
  
  # 避免除以零
  if (p1$concentration_pct == p2$concentration_pct) {
    return(p1$depth_mm)
  }
  
  # 線性內插公式: x = x1 + (y - y1) * (x2 - x1) / (y2 - y1)
  effective_depth <- p1$depth_mm + 
    (target_concentration - p1$concentration_pct) * 
    (p2$depth_mm - p1$depth_mm) / 
    (p2$concentration_pct - p1$concentration_pct)
  
  return(effective_depth)
}

# 對所有 3600 筆數據計算 Y 值
profile_paths <- all_recipes$Profile_File_Path 
base_path <- "C:/Users/USER/Desktop/PYCProfessor/CarburizationDB_Full/"

# 使用 sapply 迴圈對每個路徑計算 Y 值
Y_real_raw <- sapply(profile_paths, function(p_path) {
  full_path <- file.path(base_path, p_path) # 組合出完整路徑
  calculate_effective_depth(full_path)
})

Y_real_raw <- as.matrix(Y_real_raw)

# 檢查一下是否有計算失敗的 NA 值，並處理它們
print(paste("計算出 Y 值的數量:", length(Y_real_raw)))
print(paste("其中失敗(NA)的數量:", sum(is.na(Y_real_raw))))

# 移除那些 X 和 Y 中包含 NA 的行 (如果有的話)
valid_indices <- !is.na(Y_real_raw)
X_real_raw <- X_real_raw[valid_indices, ]
Y_real_raw <- Y_real_raw[valid_indices, , drop = FALSE]

print(paste("清理後剩餘的數據筆數:", nrow(X_real_raw)))

#-----------------------------------------------
#  X_real_raw 和 Y_real_raw 是一一對應的乾淨數據
data_full_raw <- cbind(X_real_raw, Y_real_raw)
colnames(data_full_raw)[7] <- "Effective_Depth"

# =========================================================================
# 訓練模型

print("--- 用全部 3600 筆數據訓練一個隨機森林模型 ---")
# randomForest 需要 y 是向量，所以用 as.vector() 轉換
teacher_model <- randomForest(x = X_real_raw, y = as.vector(Y_real_raw), ntree = 100)
print("教師模型訓練完成。")

# ---------------------------------------------------
# Maximin Distance 子集抽樣
# ---------------------------------------------------
# N 是想要抽出的樣本數，例如 N=100
N <- 60

print(paste("開始執行 Maximin Distance 抽樣，目標樣本數:", N))

# 1. 對 X 數據進行正規化 (Min-Max Scaling)，距離計算應在正規化空間中進行
X_min <- apply(X_real_raw, 2, min)
X_max <- apply(X_real_raw, 2, max)

# 避免分母為零的情況
scale_factors <- X_max - X_min
scale_factors[scale_factors == 0] <- 1 
X_real_scaled <- scale(X_real_raw, center = X_min, scale = scale_factors)

# 2. 初始化子集
set.seed(2025) 
remaining_indices <- 1:nrow(X_real_scaled)

# 隨機選取第一個點
start_index <- sample(remaining_indices, 1)
subset_indices <- c(start_index)
remaining_indices <- remaining_indices[remaining_indices != start_index]

# 記錄迴圈開始時間，用於計算預估剩餘時間
start_time <- Sys.time()

# 3. 迭代地選取剩下的 N-1 個點
for (i in 2:N) {
  
  # 初始化尋找下一個最佳點的變數
  best_next_index <- -1
  max_min_dist <- -Inf # 初始化為負無窮大
  
  # 從尚未被選中的點中，遍歷尋找最佳點
  for (candidate_index in remaining_indices) {
    
    # 提取候選點的座標 (正規化後)
    candidate_point <- X_real_scaled[candidate_index, ]
    
    # 提取當前子集中所有點的座標 (正規化後)
    subset_points <- X_real_scaled[subset_indices, , drop = FALSE]
    
    # 計算該候選點到子集中【所有點】的距離
    # dist() 函式可以高效計算兩個矩陣之間的距離
    # 使用歐幾里得距離 (method = "euclidean")，對應 Forrester 公式 (1.5) p=2
    distances_to_subset <- as.matrix(dist(rbind(candidate_point, subset_points)))[1, -1]
    
    # 找到這個候選點與當前子集的最小距離
    min_dist <- min(distances_to_subset)
    
    # 核心邏輯：如果這個點的最小距離，比我們目前記錄的最佳點的最小距離還要大
    # 那麼它就是新的最佳候選點
    if (min_dist > max_min_dist) {
      max_min_dist <- min_dist
      best_next_index <- candidate_index
    }
  }
  
  # 將找到的最佳點加入子集，並從候選池中移除
  subset_indices <- c(subset_indices, best_next_index)
  remaining_indices <- remaining_indices[remaining_indices != best_next_index]
  # 計算已經過的時間
  time_now <- Sys.time()
  time_elapsed_secs <- as.numeric(difftime(time_now, start_time, units = "secs"))
  
  # 迴圈從 i=2 開始，所以已經完成了 i-1 次迭代
  avg_time_per_iter <- time_elapsed_secs / (i - 1)
  
  # 預估剩餘時間
  time_remaining_secs <- avg_time_per_iter * (N - i)
  mins_remaining <- floor(time_remaining_secs / 60)
  secs_remaining <- round(time_remaining_secs %% 60, 0)
  
  # 計算進度百分比
  progress_percent <- (i / N) * 100
  
  # 組合要顯示的訊息字串
  progress_message <- sprintf(
    "進度: %5.1f%% (%d/%d) | 已耗時: %.0fs | 預計剩餘: %d分 %02d秒",
    progress_percent, i, N, time_elapsed_secs, mins_remaining, secs_remaining
  )
  
  # 使用 cat 和 \r 在同一行刷新訊息
  cat("\r", progress_message)
  # --- 進度顯示結束 ---
}

# 在迴圈結束後換行，讓後續的輸出在新的一行開始
cat("\n")

# 4. 根據最終選出的索引，從原始未正規化的數據中提取 X 和 Y
# 將變數名改為 _real，以區分後續的增強數據
X_train_raw_real <- X_real_raw[subset_indices, ]
Y_train_raw_real <- Y_real_raw[subset_indices, , drop = FALSE]

print("---------------------------------------------------")
print(paste("已成功抽出", nrow(X_train_raw_real), "筆 Maximin Distance 樣本")) 

#--------------------------------------------------------------------------------------
# --- Step3：數據增強模組

N_augment <- 40 # 要生成多少個增強點，總數 60 + 40 = 100

print(paste("數據增強生成樣本數:", N_augment))

# 初始化用於存放增強數據的矩陣
X_train_raw_augmented <- matrix(NA, nrow = N_augment, ncol = ncol(X_real_raw))
Y_train_raw_augmented <- matrix(NA, nrow = N_augment, ncol = 1)

set.seed(202509) 

for (i in 1:N_augment) {
  # 從整個 3600 點數據集中隨機挑選兩個不同的點 A 和 B
  pair_indices <- sample(1:nrow(X_real_raw), 2, replace = FALSE)
  
  point_A_X <- X_real_raw[pair_indices[1], ]
  point_B_X <- X_real_raw[pair_indices[2], ]
  
  point_A_Y <- Y_real_raw[pair_indices[1], ]
  point_B_Y <- Y_real_raw[pair_indices[2], ]
  
  # 2. 隨機生成一個 0.25 到 0.75 之間的權重 w
  #    避免生成太靠近原始點的數據
  w <- runif(1, min = 0.25, max = 0.75)
  
  # 3. 進行線性內插，生成新的中間點 C
  #    X_C = w * X_A + (1-w) * X_B
  #    Y_C = w * Y_A + (1-w) * Y_B
  new_X <- w * point_A_X + (1 - w) * point_B_X
  
  # =========================================================================
  # 對生成的 new_X 進行四捨五入修飾
  # =========================================================================
  # 溫度欄位 (第 1 和 第 4 欄) 四捨五入到整數位
  new_X[1] <- round(new_X[1], digits = 0)
  new_X[4] <- round(new_X[4], digits = 0)
  
  # 時間和碳勢欄位 (第 2, 3, 5, 6 欄) 四捨五入到小數點後一位
  new_X[2] <- round(new_X[2], digits = 1)
  new_X[3] <- round(new_X[3], digits = 1)
  new_X[5] <- round(new_X[5], digits = 1)
  new_X[6] <- round(new_X[6], digits = 1)
  # =========================================================================
  
  # 用修飾過的 new_X 去預測 Y
  new_X_df <- as.data.frame(t(new_X))
  colnames(new_X_df) <- colnames(X_real_raw)
  new_Y <- predict(teacher_model, newdata = new_X_df)
  new_Y <- round(new_Y, digits = 2)
  # 4. 將修飾後的新點存儲起來
  X_train_raw_augmented[i, ] <- new_X
  Y_train_raw_augmented[i, ] <- new_Y
}

print(paste("已成功生成", N_augment, "筆增強樣本"))

# --- Step 4: 創建最終的混合訓練集 (60+40+100) ---

# 先將真實雙階段和增強雙階合併成一個 100 筆的雙階段數據集
X_2stage_combined_raw <- rbind(X_train_raw_real, X_train_raw_augmented) # 60 + 40 = 100 筆
Y_2stage_combined_raw <- rbind(Y_train_raw_real, Y_train_raw_augmented)

print("---------------------------------------------------")
print(paste("已合併", nrow(X_2stage_combined_raw), "筆雙階段數據 (真實 + 增強)"))


# B. 創建 100 筆真實單階段數據
print("---創建 100 筆單階段數據 ---")
set.seed(2024) 

num_1stage <- 100
indices_1stage <- sample(1:nrow(X_real_raw), num_1stage)

# 從原始數據庫中提取 X 和 Y
X_1stage_raw_full <- X_real_raw[indices_1stage, ]
Y_1stage_raw <- Y_real_raw[indices_1stage, , drop = FALSE]

# 將 X 的第二階段資訊抹除
X_1stage_raw <- X_1stage_raw_full
X_1stage_raw[, 4:6] <- NA


# C. 最終合併將雙階段數據集和單階段數據集合併
X_train_raw <- rbind(X_2stage_combined_raw, X_1stage_raw) # 100 + 100 = 200 筆
Y_train_raw <- rbind(Y_2stage_combined_raw, Y_1stage_raw)
Y_train_raw <- round(Y_train_raw, digits = 2)

print(paste("最終混合數據集創建完成，總共包含", nrow(X_train_raw), "筆數據。"))


# --- Step 5: 最終正規化 ---
print("--- 對最終的 200 筆混合數據進行正規化 ---")

# 正規化輸入 X (使用全局邊界 X_min, X_max)
X_train <- scale(X_train_raw, center = X_min, scale = scale_factors)
X_train <- as.data.frame(X_train)
colnames(X_train) <- colnames(X_real_raw)

# 正規化輸出 Y (使用全局 Y_mean, Y_sd)
Y_mean <- mean(Y_real_raw)
Y_sd <- sd(Y_real_raw)
Y_train_scaled <- as.matrix((Y_train_raw - Y_mean) / Y_sd)

print("最終訓練集 X_train (已正規化) 和 Y_train_scaled (已正規化) 已準備就緒。")
print("預覽最終訓練集的前幾行 (雙階段):")
print(head(X_train))
print("預覽最終訓練集的後幾行 (應包含 NA):")
print(tail(X_train))


# --- Step 6: 儲存最終生成的訓練集 ---
print("--- Step 6: 正在儲存最終生成的訓練數據集 ---")

output_directory <- "C:/Users/USER/Desktop/PYCProfessor/CarburizationDB_Full"
x_output_file <- file.path(output_directory, "final_X_train_raw_200.csv")
y_output_file <- file.path(output_directory, "final_Y_train_raw_200.csv")

# 為了寫入 CSV 將 X 和 Y 轉換為 data.frame 並設定好欄位名
X_to_save <- as.data.frame(X_train_raw)
colnames(X_to_save) <- colnames(X_real_raw)

Y_to_save <- as.data.frame(Y_train_raw)
colnames(Y_to_save) <- "Effective_Depth"

# row.names = FALSE 不要將 R 預設的行號也寫入檔案
write.csv(X_to_save, file = x_output_file, row.names = FALSE, na = "") # na="" 讓 NA 儲存為空格
write.csv(Y_to_save, file = y_output_file, row.names = FALSE)

print(paste("成功將 X 數據儲存至:", x_output_file))
print(paste("成功將 Y 數據儲存至:", y_output_file))