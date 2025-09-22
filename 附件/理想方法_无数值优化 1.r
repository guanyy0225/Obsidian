library(ordinal)
library(data.table)

NullModel_Likelihood <- function(phenofile, outcomeCol, sampleCol, covCol = NULL) {
  # --- Part 1: 加载和准备表型数据 ---
  if (is.character(phenofile)) {
    if (!file.exists(phenofile)) stop("Phenotype file does not exist!")
    use_data <- data.table::fread(phenofile, data.table = FALSE)
  } else {
    use_data <- as.data.frame(phenofile)
  }
  
  # 移除含有NA的行，确保数据完整
  use_data <- na.omit(use_data[, c(sampleCol, outcomeCol, covCol)])
  use_data[[sampleCol]] <- as.character(use_data[[sampleCol]])
  
  # --- Part 2: 构建零模型公式 ---
  if (!is.ordered(use_data[[outcomeCol]])) {
    use_data[[outcomeCol]] <- as.ordered(use_data[[outcomeCol]])
  }
  
  if (is.null(covCol)) {
    formula_null <- as.formula(paste(outcomeCol, "~ 1"))
  } else {
    formula_null <- as.formula(paste(outcomeCol, "~", paste(covCol, collapse = " + ")))
  }
  
  # --- Part 3: 拟合有序零模型 ---
  message("Fitting the ordinal null model using ordinal::clm...")
  clm_obj <- ordinal::clm(formula = formula_null, data = use_data, link = "probit", model = TRUE)
  
  message("Ordinal null model fitting complete.")
  
  # --- Part 4: 准备返回对象 ---
  # clm返回的对象已经包含了大部分我们需要的信息
  # 我们返回模型对象本身，以及清理后的数据框
  # model.frame(clm_obj) 包含了模型实际使用的数据
  fit_null <- list(
    clm_fit = clm_obj,
    model_data = model.frame(clm_obj),
    sample_ids = use_data[rownames(model.frame(clm_obj)), sampleCol],
    outcome = clm_obj$y,
    formula_null = formula(clm_obj)
  )
  
  return(fit_null)
}






library(numDeriv)

#' 对单个遗传变异执行得分检验
#'
#' @param null_model_fit `NullModel_Likelihood` 函数的输出
#' @param G 单个遗传变异的基因型向量 (与 null_model_fit$sample_ids 对应)
#' @return 包含得分统计量和p值的列表
OrdinalScoreTest <- function(null_model_fit, G) {
  
  # --- Step 1: 提取零模型信息 ---
  clm_fit <- null_model_fit$clm_fit
  y <- null_model_fit$outcome
  
  # 零模型的设计矩阵 (不包含截距)
  X_null_no_intercept <- model.matrix(null_model_fit$formula_null, data = null_model_fit$model_data)[, -1, drop = FALSE]
  
  # 零模型下的参数估计值
  alpha_tilde <- clm_fit$alpha
  beta_tilde <- clm_fit$beta # 这是协变量的系数
  
  # --- Step 2: 构建完整模型的设计矩阵和参数向量 (在H0下) ---
  # 完整模型的设计矩阵 (协变量 + 基因型G)，不含截距
  X_full_no_intercept <- cbind(X_null_no_intercept, G)
  
  # 完整模型在H0下的参数向量 (alpha, beta_cov, beta_G=0)
  params_h0 <- c(alpha_tilde, beta_tilde, 0)
  names(params_h0) <- c(names(alpha_tilde), names(beta_tilde), "G")
  
  # --- Step 3: 计算得分向量 (Gradient) ---
  # 我们需要一个包装函数，因为 `grad` 需要一个只接受参数的函数
  # 注意：loglike_probit 返回的是 *负* 对数似然，所以梯度也是负的。
  # 我们需要的是对数似然的梯度，所以要取反。
  grad_func <- function(params) {
    # 重新定义 loglike_probit 以符合 numDeriv 的要求
    # 这里的 logL 返回正值
    n_alpha <- nlevels(y) - 1
    alpha <- params[1:n_alpha]
    beta <- params[-(1:n_alpha)]
    if (is.unsorted(alpha)) return(-Inf) # 返回-Inf表示无效参数
    
    thresholds <- c(-Inf, alpha, Inf)
    eta <- as.vector(X_full_no_intercept %*% beta)
    y_idx <- as.numeric(y)
    
    upper_p <- pnorm(thresholds[y_idx + 1] - eta)
    lower_p <- pnorm(thresholds[y_idx] - eta)
    
    prob <- upper_p - lower_p
    prob[prob <= 0] <- 1e-100
    
    return(sum(log(prob)))
  }
  
  message("Calculating Score vector (Gradient)...")
  U <- grad(func = grad_func, x = params_h0)
  
  # --- Step 4: 计算信息矩阵 (Hessian) ---
  message("Calculating Information matrix (Hessian)...")
  # Hessian of the log-likelihood
  hessian_mat <- hessian(func = grad_func, x = params_h0)
  # Information matrix is the negative of the Hessian
  I <- -hessian_mat
  
  # 检查信息矩阵是否可逆
  # 使用 tryCatch 来处理奇异矩阵 (singular matrix)
  I_inv <- tryCatch({
    solve(I)
  }, error = function(e) {
    # 如果不可逆，可能存在共线性问题
    # 使用广义逆 Moore-Penrose generalized inverse
    warning("Information matrix is singular, using generalized inverse.")
    MASS::ginv(I)
  })
  
  # --- Step 5: 计算得分检验统计量 ---
  # S = U' * I_inv * U
  score_stat <- t(U) %*% I_inv %*% U
  score_stat <- as.numeric(score_stat)
  
  # 在H0下，得分统计量服从自由度为1的卡方分布
  # (因为我们只检验了一个参数 beta_G)
  p_value <- pchisq(score_stat, df = 1, lower.tail = FALSE)
  
  return(list(score_statistic = score_stat, p_value = p_value, score_vector = U))
}


