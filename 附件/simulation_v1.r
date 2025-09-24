# --- Step 0: 加载所有需要的库和函数 ---
# 确保你已经将 NullModel_Ordinal 和 Ordinal_ScoreTest 函数加载到你的R环境中

# 确保加载所有需要的库
library(ordinal)
library(Matrix)
library(data.table)

# ==============================================================================
# 这里粘贴你提供的 NullModel_Ordinal 函数
# ==============================================================================
#' 拟合有序多分类零模型 (修正版 - 使用标准残差和权重)
#'
#' @description
#' 这个修正版的函数使用标准的广义残差 (y - E[y]) 和权重 (Var[y])
#' 来替代原来复杂的手动计算，以解决假阳性膨胀的问题。
#'
#' 拟合有序多分类零模型 (真正最终版 - 手动计算概率)
#'
#' @description
#' 放弃使用 fitted.values 和 predict，改为根据模型系数手动计算概率。
#' 这是最透明、最稳健的方法，可以彻底解决所有与概率矩阵结构相关的问题。
#'
NullModel_Ordinal <- function(phenofile, outcomeCol, sampleCol,
                              covCol = NULL, PRSCol = NULL, verbose = FALSE) {
  # --- Part 1, 1.5, 2 (保持不变) ---
  if (is.character(phenofile)) {
    if (!file.exists(phenofile)) stop("Phenotype file does not exist!")
    use_data <- as.data.frame(data.table::fread(phenofile, data.table = FALSE))
  } else {
    use_data <- as.data.frame(phenofile)
  }
  all_cols_to_check <- c(sampleCol, outcomeCol, covCol, PRSCol)
  cols_exist <- all_cols_to_check[all_cols_to_check %in% colnames(use_data)]
  use_data <- use_data[, cols_exist]
  use_data <- use_data[complete.cases(use_data), ]
  use_data[[sampleCol]] <- as.character(use_data[[sampleCol]])
  if(any(duplicated(use_data[[sampleCol]]))) stop("Sample IDs are not unique!")
  rownames(use_data) <- use_data[[sampleCol]]
  all_covars_original <- c(covCol, PRSCol)
  if(length(all_covars_original) > 0){
    is_factor_or_char <- sapply(use_data[, all_covars_original, drop = FALSE], function(x) is.factor(x) || is.character(x))
    factor_covar_names <- all_covars_original[is_factor_or_char]
    numeric_covar_names <- all_covars_original[!is_factor_or_char]
    data_for_model <- use_data
    all_covars_final <- numeric_covar_names
    if(length(factor_covar_names) > 0) {
      dummy_formula <- as.formula(paste("~", paste(factor_covar_names, collapse=" + ")))
      dummy_vars <- model.matrix(dummy_formula, data = use_data)[, -1, drop=FALSE]
      colnames(dummy_vars) <- gsub("\\.", "_", make.names(colnames(dummy_vars)))
      data_for_model <- cbind(data_for_model, dummy_vars)
      all_covars_final <- c(all_covars_final, colnames(dummy_vars))
      data_for_model <- data_for_model[, !colnames(data_for_model) %in% factor_covar_names, drop = FALSE]
    }
  } else {
    data_for_model <- use_data
    all_covars_final <- NULL
  }
  if (!is.ordered(data_for_model[[outcomeCol]])) {
    data_for_model[[outcomeCol]] <- as.ordered(data_for_model[[outcomeCol]])
  }
  if (is.null(all_covars_final) || length(all_covars_final) == 0) {
    formula_string <- paste0("`", outcomeCol, "` ~ 1")
  } else {
    safe_covars <- paste0("`", all_covars_final, "`")
    formula_string <- paste0("`", outcomeCol, "` ~ ", paste(safe_covars, collapse = " + "))
  }
  formula_null <- as.formula(formula_string)
  message("Fitting the ordinal null model using ordinal::clm with logit link...")
  clm_obj  <- ordinal::clm(formula = formula_null, data = data_for_model, model = TRUE, link="logit")
  if(verbose) print(summary(clm_obj))

  # --- Part 3: 预计算 (手动计算概率版) ---
  message("Pre-computing components for the score test (manually calculating probabilities)...")

  model_data <- clm_obj$model
  sample_ids <- rownames(model_data)
  n <- clm_obj$n
  J <- length(clm_obj$y.levels)

  y_numeric <- as.numeric(clm_obj$y)
  y_levels_numeric <- 1:J

  # --- 手动计算概率 ---
  # 提取系数
  beta_hat <- clm_obj$beta      # 协变量系数
  theta_hat <- clm_obj$alpha   # 切点 (intercepts)

  # 构造不含截距的协变量设计矩阵
  X_covars_only <- model.matrix(delete.response(terms(clm_obj)), model_data)
  if ("(Intercept)" %in% colnames(X_covars_only)) {
    X_covars_only <- X_covars_only[, -1, drop = FALSE]
  }

  # 计算线性预测值 eta = X * beta
  eta <- rep(0, n)
  if (!is.null(beta_hat) && length(beta_hat) > 0) {
    # 确保矩阵和向量的名称/顺序对齐
    if (!identical(colnames(X_covars_only), names(beta_hat))) {
      X_covars_only <- X_covars_only[, names(beta_hat), drop = FALSE]
    }
    eta <- as.vector(X_covars_only %*% beta_hat)
  }

  # 计算累积概率 P(Y <= j) = logistic(theta_j - eta)
  cum_probs <- matrix(NA, nrow = n, ncol = J - 1)
  for (j in 1:(J - 1)) {
    cum_probs[, j] <- plogis(theta_hat[j] - eta)
  }

  # 计算类别概率 P(Y = j)
  fitted_probs <- matrix(NA, nrow = n, ncol = J)
  fitted_probs[, 1] <- cum_probs[, 1]
  for (j in 2:(J - 1)) {
    fitted_probs[, j] <- cum_probs[, j] - cum_probs[, j - 1]
  }
  fitted_probs[, J] <- 1 - cum_probs[, J - 1]

  # --- 后续计算现在可以安全进行 ---
  expected_y <- fitted_probs %*% y_levels_numeric
  residuals_simple <- y_numeric - expected_y
  expected_y_sq <- fitted_probs %*% (y_levels_numeric^2)
  weights_simple <- expected_y_sq - (expected_y^2)

  eps <- 1e-8
  weights_simple[weights_simple < eps] <- eps

  X_mat <- model.matrix(formula(clm_obj), data = model_data)
  if (!identical(rownames(X_mat), sample_ids)) {
    rownames(X_mat) <- sample_ids
  }

  W_mat_diag <- Diagonal(x = as.numeric(weights_simple))
  X_t_W <- crossprod(X_mat, W_mat_diag)
  XWX_mat <- X_t_W %*% X_mat
  XWX_inv <- solve(XWX_mat + diag(eps, ncol(XWX_mat)))

  # --- Part 4: 组装最终返回对象 ---
  fit_null <- list(
    clm_fit = clm_obj,
    residuals = as.numeric(residuals_simple),
    W_diag = as.numeric(weights_simple),
    sample_ids = sample_ids,
    X_mat = X_mat,
    X_t_W_X_inv = XWX_inv,
    num_unstable_samples = 0
  )

  message("Ordinal null model fitting complete.")
  return(fit_null)
}


# --- Step 1: 生成模拟数据 ---

cat("--- Step 1: Generating Simulation Data ---\n\n")

set.seed(2023)
n_samples <- 2000
J_categories <- 4 # 4个有序类别

# 生成协变量
cov1 <- rnorm(n_samples)
cov2_factor <- factor(sample(c("A", "B", "C"), n_samples, replace = TRUE))

# 生成遗传变异
# SNP_causal: 常见变异 (MAF=20%), 与表型相关
# SNP_null:   稀有变异 (MAF=1%), 与表型无关
maf_causal <- 0.20
maf_null <- 0.01
SNP_causal <- rbinom(n_samples, 2, maf_causal)
SNP_null <- rbinom(n_samples, 2, maf_null)

# 真实参数 (用于数据生成)
beta_cov1 <- 0.5
beta_cov2B <- -0.3 # 因子B相对于A的效应
beta_cov2C <- 0.2  # 因子C相对于A的效应
gamma_causal <- 0.4 # SNP_causal 的真实效应 (OR ≈ 1.5)
gamma_null <- 0.0   # SNP_null 的真实效应为零

# 生成线性预测变量
# 手动创建因子协变量的 design matrix for data generation
cov2_B_dummy <- as.numeric(cov2_factor == "B")
cov2_C_dummy <- as.numeric(cov2_factor == "C")

linear_pred <- (cov1 * beta_cov1 +
                cov2_B_dummy * beta_cov2B +
                cov2_C_dummy * beta_cov2C +
                SNP_causal * gamma_causal +
                SNP_null * gamma_null)

# 生成有序表型，使其分布不均衡
true_intercepts <- c(-1.5, 0, 1.5) # J-1 个切点

prob <- matrix(0, nrow = n_samples, ncol = J_categories)
cum_prob <- matrix(0, nrow = n_samples, ncol = J_categories)

for (j in 1:(J_categories - 1)) {
  cum_prob[, j] <- plogis(true_intercepts[j] - linear_pred)
}
cum_prob[, J_categories] <- 1

prob[, 1] <- cum_prob[, 1]
for (j in 2:J_categories) {
  prob[, j] <- cum_prob[, j] - cum_prob[, j - 1]
}

phenotype_ord <- apply(prob, 1, function(p) sample(1:J_categories, 1, prob = p))

# 组装成数据框
# 在 Step 1 生成模拟数据部分
sim_data <- data.frame(
  sampleID = paste0("ID_", 1:n_samples),
  phenotype = factor(phenotype_ord, ordered = TRUE),
  cov1 = cov1,
  cov2_factor = cov2_factor,
  SNP_causal = SNP_causal,
  SNP_null = SNP_null,
  row.names = paste0("ID_", 1:n_samples) # <--- 增加这一行
)

cat("Phenotype distribution:\n")
print(table(sim_data$phenotype))
cat("\nCausal SNP (MAF=", maf_causal, ") distribution:\n", sep="")
print(table(sim_data$SNP_causal))
cat("\nNull SNP (MAF=", maf_null, ") distribution:\n", sep="")
print(table(sim_data$SNP_null))
cat("\n")

# --- Step 2: 拟合零模型 ---

cat("--- Step 2: Fitting the Null Model using NullModel_Ordinal ---\n\n")

null_model_fit <- NullModel_Ordinal(
  phenofile = sim_data,
  outcomeCol = "phenotype",
  sampleCol = "sampleID",
  covCol = c("cov1", "cov2_factor"),
  verbose = TRUE
)

# --- Step 3: 执行 Score 检验 ---

cat("\n--- Step 3: Performing Score Tests ---\n\n")

# 确保基因型向量与模型样本对齐
# NullModel_Ordinal 内部处理了 complete.cases, 所以我们需要从它的输出中获取对齐的数据
model_samples <- null_model_fit$sample_ids
aligned_G_causal <- sim_data[model_samples, "SNP_causal"]
aligned_G_null <- sim_data[model_samples, "SNP_null"]

# 检验 Causal SNP
cat("Testing Causal SNP (SNP_causal, MAF=", maf_causal, "):\n", sep="")
score_result_causal <- Ordinal_ScoreTest(null_model_fit, aligned_G_causal)
print(score_result_causal)

# 检验 Null SNP
cat("\nTesting Null SNP (SNP_null, MAF=", maf_null, "):\n", sep="")
score_result_null <- Ordinal_ScoreTest(null_model_fit, aligned_G_null)
print(score_result_null)


# --- Step 4: 结果验证 (与 Wald 检验比较) ---

cat("\n--- Step 4: Verifying with Wald Tests from Full Models ---\n\n")

# 验证 Causal SNP
cat("Fitting full model for Causal SNP (SNP_causal)...\n")
full_model_causal <- ordinal::clm(
  phenotype ~ cov1 + cov2_factor + SNP_causal,
  data = sim_data,
  Hess = TRUE
)
summary_causal <- summary(full_model_causal)
p_wald_causal <- summary_causal$coefficients["SNP_causal", "Pr(>|z|)"]
cat("Wald test p-value for SNP_causal:", p_wald_causal, "\n")
cat("Score test p-value for SNP_causal:", score_result_causal$p_value_norm, "\n")


# 验证 Null SNP
cat("\nFitting full model for Null SNP (SNP_null)...\n")
full_model_null <- ordinal::clm(
  phenotype ~ cov1 + cov2_factor + SNP_null,
  data = sim_data,
  Hess = TRUE
)
summary_null <- summary(full_model_null)
p_wald_null <- summary_null$coefficients["SNP_null", "Pr(>|z|)"]
cat("Wald test p-value for SNP_null:", p_wald_null, "\n")
cat("Score test p-value for SNP_null:", score_result_null$p_value_norm, "\n")

