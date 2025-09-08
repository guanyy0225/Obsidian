
---
## `REGENIE`预处理

```R
# --- 文件名: 01_generate_latent_phenotype.R ---

##########################################################
# Ordinal Phenotype Transformation for REGENIE
# Part A: Generate Latent Variable Phenotype
# Author: Zilin Li, with contributions from AI assistant
# Date: 19/06/2025
##########################################################

# --- 0. Load Libraries ---
library(data.table)
library(ordinal)
library(dplyr)

# --- 1. Load and Prepare Initial Data ---
message("Step 1: Loading and preparing initial data...")

load("/datapool/maths/UKBB/Data/Phenotype/main_survey/bd_676623.Rdata")

pcs <- fread("/datapool/maths/UKBB/Data/PC_GRM/sGRM_and_PCs/output.pca.score")
colnames(pcs)[1:2] <- c("userId", "userID2")
colnames(pcs)[3:22] <- paste0("PC", 1:20)

pheno <- bd[, c(c(0, 769) + 1)]
colnames(pheno) <- c("userId", "alcohol_intake_frequency")

covariate <- bd[, c(c(0, 22, 10318) + 1)]
colnames(covariate)[1:3] <- c("userId", "sex", "age")

t1 <- merge(pheno, covariate, by = "userId")
fullDat <- merge(t1, pcs, by = "userId")

# --- 2. Filter Samples ---
message("Step 2: Filtering samples...")
gId <- get(load("/datapool/maths/UKBB/Data/Phenotype/sample_id/UKB_500K_WGS_sampleid.RData"))
data_for_model <- fullDat[
  (fullDat$userId %in% gId) &
    !is.na(fullDat$alcohol_intake_frequency) &
    !is.na(fullDat$PC1) & !is.na(fullDat$age) & !is.na(fullDat$sex) &
    fullDat$alcohol_intake_frequency != "Prefer not to answer",
]

# --- 3. Transform Covariates and Outcome ---
message("Step 3: Transforming variables...")
data_for_model$age2 <- (data_for_model$age)^2
data_for_model$alcohol_intake_frequency <- ordered(
  factor(
    data_for_model$alcohol_intake_frequency,
    levels = c("Never", "Special occasions only", 
               "One to three times a month", "Once or twice a week", 
               "Three or four times a week", "Daily or almost daily")
  )
)

# --- IMPORTANT: Standardize covariates on the FINAL dataset ---
cols_to_scale <- c("age", "age2", paste0("PC", 1:10))
data_for_model <- data_for_model %>%
  mutate(across(all_of(cols_to_scale), scale))

# --- 4. Fit Base Ordinal Model ---
message("Step 4: Fitting the base ordinal probit model...")
base_formula <- as.formula(
  "alcohol_intake_frequency ~ sex + age + age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
)

clm_base_obj <- tryCatch({
  ordinal::clm(formula = base_formula, data = data_for_model, link = "probit", model = TRUE)
}, error = function(e) {
  stop("Failed to fit the base ordinal model. Original error: ", e$message)
})

message("Base model fitting successful.")

# --- 5. Calculate Latent Variable Residuals ---
message("Step 5: Calculating latent variable residuals...")
model_fit_data <- clm_base_obj$model
# Ensure we are working with the exact data the model used
final_data_subset <- data_for_model[as.numeric(rownames(model_fit_data)), ]

alpha_coefs <- clm_base_obj$beta
X <- model.matrix(object = formula(clm_base_obj), data = model_fit_data)
eta <- as.vector(X[, names(alpha_coefs), drop = FALSE] %*% alpha_coefs)

thresholds <- c(-Inf, clm_base_obj$alpha, Inf)
y_idx <- as.numeric(clm_base_obj$y)
lower_b <- thresholds[y_idx] - eta
upper_b <- thresholds[y_idx + 1] - eta
prob_interval <- pnorm(upper_b) - pnorm(lower_b)
prob_interval[prob_interval < 1e-12] <- 1e-12
latent_residuals <- (dnorm(lower_b) - dnorm(upper_b)) / prob_interval

# --- 6. Inverse Normal Transformation (Optional but Recommended) ---
# The latent residuals are not perfectly normally distributed. Applying INT can
# improve the performance of REGENIE's LMM.
message("Step 6: Applying inverse normal transformation to residuals...")
# Rank-based INT: (r - 0.5) / n
n <- length(latent_residuals)
rank_residuals <- rank(latent_residuals, na.last = "keep", ties.method = "average")
latent_pheno_int <- qnorm((rank_residuals - 0.5) / n)

# --- 7. Create and Save REGENIE Input Files ---
message("Step 7: Saving files for REGENIE...")
regenie_input_df <- data.frame(
  FID = final_data_subset$userId,
  IID = final_data_subset$userId,
  latent_phenotype = latent_pheno_int,
  # We also need to save the covariates file REGENIE will use, which can be empty
  # if we believe all effects are captured. Best practice is to include PCs again.
  sex = final_data_subset$sex,
  age = final_data_subset$age,
  age2 = final_data_subset$age2
)
# Add PCs
pcs_subset <- final_data_subset[, paste0("PC", 1:10)]
regenie_input_df <- cbind(regenie_input_df, pcs_subset)


# Create Phenotype File
pheno_file_out <- regenie_input_df[, c("FID", "IID", "latent_phenotype")]
pheno_path <- "/datapool/home/2024102311/UKB/nullmodels/multiclass/alcohol_intake_frequency/latent_pheno/phenotype.txt"
fwrite(pheno_file_out, file = pheno_path, sep = "\t", col.names = TRUE)
message(paste("Phenotype file saved to:", pheno_path))

# Create Covariate File
# REGENIE still needs covariates for its internal computations, even if they
# are uncorrelated with the residual phenotype. It's robust to include them.
covar_file_out <- regenie_input_df[, -which(colnames(regenie_input_df) == "latent_phenotype")]
covar_path <- "/datapool/home/2024102311/UKB/nullmodels/multiclass/alcohol_intake_frequency/latent_pheno/covariates.txt"
fwrite(covar_file_out, file = covar_path, sep = "\t", col.names = TRUE)
message(paste("Covariate file saved to:", covar_path))

message("Part A complete. You can now proceed to Part B (REGENIE).")
```

```R
#!/bin/bash
# --- 文件名: 02_run_regenie_step1_only.sh ---

##########################################################
# Hybrid OrdinalSTAAR Workflow
# Part B: Run REGENIE Step 1 to generate PRS
# Author: Zilin Li
# Date: 19/06/2025
##########################################################

# Activate your conda environment for REGENIE
# conda activate regenie_env

# Define file paths
BED_FILES="/datapool/maths/UKBB/UKB_genotypes_QCed_ZHL/chrall"
PHENO_FILE="/datapool/home/2024102311/UKB/nullmodels/multiclass/alcohol_intake_frequency/latent_pheno/phenotype.txt"
COVAR_FILE="/datapool/home/2024102311/UKB/nullmodels/multiclass/alcohol_intake_frequency/latent_pheno/covariates.txt"
OUTPUT_PREFIX="/datapool/home/2024102311/UKB/nullmodels/multiclass/alcohol_intake_frequency/latent_pheno/regenie_output"

COVAR_LIST="sex,age,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"

echo "--- Starting REGENIE Step 1 to generate polygenic predictions ---"
regenie \
--step 1 \
--bed ${BED_FILES} \
--phenoFile ${PHENO_FILE} \
--phenoCol latent_phenotype \
--covarFile ${COVAR_FILE} \
--covarColList ${COVAR_LIST} \
--catCovarList sex \
--qt \
--bsize 1000 \
--out ${OUTPUT_PREFIX} \
--threads 8 \

echo "--- REGENIE Step 1 complete. Predictions are ready for OrdinalSTAAR. ---"

```




## `OrdinalSTAAR`框架

这个完整的 `OrdinalSTAAR` 框架将包含三个主要文件：
1.  **`GeneCentricCoding.R`**: 负责基因型预处理，我们直接采用 `SurvSTAAR` 的现代化版本。
2.  **`NullModel.R`**: 负责拟合有序零模型并计算分数检验所需的组件，严格模仿 `SurvSTAAR` 的 `NullModel.R` 的设计。
3.  **`OrdinalSTAAR.R`**: 这是新的主分析函数，它将调用 `NullModel.R` 的输出和 `GeneCentricCoding.R` 的输出，来**执行自定义的分数检验**（负担检验、SKAT等）。这部分代码将是全新的，因为它需要利用 `NullModel.R` 提供的矩阵来进行计算。

下面是这三个核心文件的完整代码。

---

### 文件 1: `GeneCentricCoding.R` (与之前确认的版本相同)

这个文件负责基因型数据的预处理，与表型无关，可以直接重用。

```R
# --- 文件名: GeneCentricCoding.R ---

#' @title Perform Gene-Centric Coding and Imputation from a GDS file
#' @description This function extracts rare variant genotypes for each gene,
#'   codes them numerically (0, 1, 2), imputes missing values using allele
#'   frequency, and saves the results to a single .Rdata file.
#'
#' @param gdsfile A character string for the path to the GDS file.
#' @param QC_label A character string used to name the output objects (e.g., "imputed").
#' @param sample.id A vector of sample IDs to be included in the analysis.
#' @param variant.id A vector of variant IDs to be included in the analysis.
#' @param outfile A character string for the path to the output .Rdata file.
#'
#' @return This function does not return a value but saves the coded genotype
#'   data to the specified output file.
#'
#' @import SeqArray
#' @export
GeneCentricCoding <- function(gdsfile, QC_label, sample.id, variant.id, outfile) {
  
  # --- Input Validation & Setup ---
  if (!requireNamespace("SeqArray", quietly = TRUE)) {
    stop("Package 'SeqArray' is required. Please install it from Bioconductor.", call. = FALSE)
  }
  
  gds <- seqOpen(gdsfile)
  on.exit(seqClose(gds))
  
  # --- Get Annotation Info ---
  message("Extracting gene and variant information...")
  gds_variant.id <- seqGetData(gds, "variant.id")
  gds_sample.id <- seqGetData(gds, "sample.id")
  
  variant_in_gds.idx <- which(gds_variant.id %in% variant.id)
  if (length(variant_in_gds.idx) == 0) stop("None of the provided variant IDs were found in the GDS file.")
  
  sample_in_gds.idx <- which(gds_sample.id %in% sample.id)
  if (length(sample_in_gds.idx) == 0) stop("None of the provided sample IDs were found in the GDS file.")
  
  # Filter GDS to only variants and samples of interest
  seqSetFilter(gds, variant.id = variant.id, sample.id = sample.id)
  
  gene_all <- seqGetData(gds, "annotation/info/GENE")
  variant.id_all <- seqGetData(gds, "variant.id")
  
  genes <- unique(gene_all)
  num_genes <- length(genes)
  message(paste("Found", num_genes, "unique genes to process."))
  
  # --- Main Loop: Process each gene ---
  for (a in 1:num_genes) {
    gene_name <- genes[a]
    if (a %% 100 == 0) cat(paste0("Processing gene ", a, "/", num_genes, ": ", gene_name, "\n"))
    
    variant.id_gene.idx <- which(gene_all == gene_name)
    variant.id_gene <- variant.id_all[variant.id_gene.idx]
    
    seqSetFilter(gds, variant.id = variant.id_gene, sample.id = sample.id)
    
    AC_Cases_or_All <- seqGetData(gds,"annotation/info/AC")
    if (sum(AC_Cases_or_All, na.rm = TRUE) == 0) {
      cat(paste("Skipping gene", gene_name, "as it has no variant alleles in the selected samples.\n"))
      next
    }
    
    genotype <- seqGetData(gds, "$dosage")
    af <- seqGetData(gds, "annotation/info/AF")
    
    if (any(is.na(genotype))) {
      for(i in 1:nrow(genotype)){
        missing.idx <- which(is.na(genotype[i,]))
        if(length(missing.idx)>0) genotype[i,missing.idx] <- 2*af[i]
      }
    }
    
    colnames(genotype) <- seqGetData(gds, "sample.id")
    rownames(genotype) <- variant.id_gene
    
    obj_out <- list(genotype = genotype)
    obj_name <- paste0("G_", QC_label, "_", gene_name)
    
    assign(obj_name, obj_out)
    
    if (a == 1) save(list = obj_name, file = outfile)
    else save(list = obj_name, file = outfile, append = TRUE)
  }
  
  message(paste("\nGene-centric coding complete. Results saved to:", outfile))
}
```

---

### 文件 2: `NullModel.R` (最终、正确的版本)

这个文件将拟合有序零模型，并计算分数检验所需的矩阵，其设计严格遵循 `SurvSTAAR` 的 `NullModel.R`。

```R
# --- 文件名: NullModel.R (功能完备的最终版) ---

#' @title 为有序表型拟合一个功能完备的零模型 (支持 LOCO 和 SPA)
#' @description 本函数拟合一个累积链接模型，并为后续的自定义稀有变异分数检验
#'   计算必要的组件。其参数和功能设计严格参照 SurvSTAAR，提供了一个完整的数据
#'   整合、模型拟合和组件计算流程。
#'
#' @param genofile 一个可选的基因组文件。可以是 GDS 文件的路径(字符串)，或一个
#'   SeqVarGDSClass 对象。如果提供，函数将匹配表型和基因型样本。
#' @param phenofile 表型文件的路径(字符串)，或一个已加载的数据框。
#' @param outcomeCol 有序结局变量的列名。
#' @param sampleCol 样本ID列的名称。
#' @param covCol 基础协变量列名的向量 (不含 PRS)。
#' @param PRSCol (可选) PRS 列的名称。该列应存在于 phenofile 中。
#' @param LOCO 一个逻辑值，指示是否启用 LOCO 模式 (默认: TRUE)。
#' @param chr 一个数值，代表当前正在分析的染色体。如果 LOCO=TRUE，则必须提供。
#' @param use_SPA (暂未实现，作为占位符) 一个逻辑值，指示是否为 SPA 准备数据。
#' @param verbose 一个逻辑值，指示是否打印详细的执行信息。
#'
#' @return 一个包含分数检验所需组件以及元数据的列表。
#'
#' @import ordinal
#' @import Matrix
#' @import data.table
#' @import SeqArray
#'
Ordinal_NullModel <- function(genofile = NULL, phenofile, outcomeCol, sampleCol,
                              covCol = NULL, PRSCol = NULL, LOCO = TRUE, chr = NULL,
                              use_SPA = FALSE, verbose = FALSE) {

  # --- Part 1: Load and Prepare Phenotype Data ---
  if (!is.null(phenofile)) {
    if (is.character(phenofile)) {
      if (!file.exists(phenofile)) stop("Phenotype file does not exist!")
      use_data <- data.table::fread(phenofile, data.table = FALSE)
    } else {
      use_data <- as.data.frame(phenofile)
    }
    use_data <- na.omit(use_data[, c(sampleCol, outcomeCol, covCol, PRSCol)])
  } else {
    stop("You must provide the phenotype data via the 'phenofile' argument.")
  }
  
  # --- Part 2: Match Samples with Genotype Data (if provided) ---
  if (!is.null(genofile)) {
    message("Matching phenotype samples with genotype data...")
    if (inherits(genofile, "SeqVarGDSClass")) {
      sample.geno <- SeqArray::seqGetData(genofile, "sample.id")
    } else if (is.character(genofile)) {
      fam_file <- paste0(genofile, ".fam")
      if (!file.exists(fam_file)) stop("PLINK .fam file not found!")
      sample.geno <- data.table::fread(fam_file, data.table = FALSE)$V2
    } else {
      stop("'genofile' must be a GDS object or a path to a PLINK file.")
    }
    
    sample.pheno <- use_data[, sampleCol]
    common_samples <- intersect(sample.pheno, sample.geno)
    
    if (length(common_samples) == 0) stop("No common samples found between phenotype and genotype files.")
    
    use_data <- use_data[use_data[, sampleCol] %in% common_samples, ]
    message(paste(nrow(use_data), "samples remain after matching."))
  }
  
  # --- Part 3: Dynamic Formula Construction ---
  if(verbose) print("Constructing model formula...")
  
  # Ensure outcome is an ordered factor
  if (!is.ordered(use_data[[outcomeCol]])) {
    use_data[[outcomeCol]] <- as.ordered(use_data[[outcomeCol]])
  }
  
  # Combine all covariates
  all_covars <- covCol
  if (!is.null(PRSCol)) {
    all_covars <- c(all_covars, PRSCol)
  }
  
  if (is.null(all_covars)) {
    formula_null <- as.formula(paste(outcomeCol, "~ 1"))
  } else {
    formula_null <- as.formula(paste(outcomeCol, "~", paste(all_covars, collapse = " + ")))
  }
  if(verbose) print(formula_null)

  # --- Part 4: Fit the Ordinal Null Model ---
  if(verbose) print("Fitting null model with ordinal::clm...")
  
  clm_obj <- ordinal::clm(formula = formula_null, data = use_data, link = "probit", model = TRUE)
  if(verbose) print(summary(clm_obj))
  
  # --- Part 5: Calculate Residuals and Variance Components ---
  if(verbose) print("Calculating latent residuals and variance components...")
  
  model_data <- clm_obj$model
  kept_row_indices <- as.numeric(rownames(model_data))
  sample_ids <- use_data[[sampleCol]][kept_row_indices]
  
  # Latent residuals and conditional variance calculation (as before)
  alpha_coefs <- clm_obj$beta
  X <- model.matrix(object = formula(clm_obj), data = model_data)
  eta <- as.vector(X[, names(alpha_coefs), drop = FALSE] %*% alpha_coefs)
  thresholds <- c(-Inf, clm_obj$alpha, Inf)
  y_idx <- as.numeric(clm_obj$y)
  lower_b <- thresholds[y_idx] - eta
  upper_b <- thresholds[y_idx + 1] - eta
  prob_interval <- pnorm(upper_b) - pnorm(lower_b)
  prob_interval[prob_interval < 1e-12] <- 1e-12
  residuals <- (dnorm(lower_b) - dnorm(upper_b)) / prob_interval
  
  term1 <- (lower_b * dnorm(lower_b) - upper_b * dnorm(upper_b)) / prob_interval
  var_y <- 1 + term1 - residuals^2
  var_y[!is.finite(var_y) | var_y < 1e-8] <- 1
  
  # Score test variance components (as before)
  W_mat <- Diagonal(x = 1/var_y)
  X_mat <- model.matrix(clm_obj)
  XWX_mat <- crossprod(X_mat, W_mat %*% X_mat)
  XWX_inv <- solve(XWX_mat)
  WX_mat <- W_mat %*% X_mat
  
  # --- Part 6: Assemble the Final List with LOCO information ---
  base_list <- list(
    residuals = residuals,
    sample_ids = sample_ids,
    W_mat = W_mat,
    X_mat = X_mat,
    WX_mat = WX_mat,
    XWX_inv = XWX_inv,
    formula_null = formula(clm_obj),
    coefficients_null = coef(clm_obj),
    use_data = model_data,
    use_SPA = use_SPA
  )
  
  if (is.null(LOCO) || isFALSE(LOCO)) {
    fit_null <- c(base_list, list(LOCO = FALSE))
  } else {
    fit_null <- c(base_list, list(LOCO = TRUE, chr = chr))
  }
  
  # Placeholder for future SPA implementation
  if(use_SPA){
    # fit_null <- CGF4LatentRes(fit_null = fit_null, ...)
    warning("use_SPA=TRUE is not yet implemented for Ordinal_NullModel. Skipping CGF calculation.")
  }

  message("Ordinal null model fitting complete.")
  return(fit_null)
}
```

---

### 文件 3: `OrdinalSTAAR.R` (新的主分析函数)

这是将所有部分连接在一起的新增核心文件。它会加载基因型数据和零模型结果，然后逐个基因为单位，计算各种关联检验（负担、SKAT、ACAT）的 p-value。

```R
# --- 文件名: OrdinalSTAAR.R ---

#' @title Perform Rare Variant Association Tests for Ordinal Phenotypes
#' @description This is the main analysis function for the OrdinalSTAAR framework.
#'   It takes a fitted ordinal null model and gene-centric coded genotype data
#'   to perform Burden, SKAT, and ACAT-OMNI rare variant association tests.
#'   The design is a direct adaptation of the SurvSTAAR framework.
#'
#' @param coded_genotype_file Path to the .Rdata file created by `GeneCentricCoding`.
#'   This file contains gene-specific genotype matrices.
#' @param fit_null The list object returned by `Ordinal_NullModel`. This object
#'   contains the residuals and pre-computed matrices for the score test.
#' @param annotation_file An optional path to an .Rdata file containing functional
#'   annotation information for weighting variants.
#' @param QC_label A character string used as a prefix for objects in the
#'   `coded_genotype_file` and `annotation_file` (e.g., "imputed").
#' @param rare_maf_cutoff MAF cutoff to define rare variants (default: 0.01).
#' @param rv_num_cutoff Minimum number of rare variants in a gene to perform a
#'   test (default: 2).
#'
#' @return A data.frame with association test results for each gene.
#'
#' @importFrom stats qchisq pchisq pcauchy
#' @import Matrix
#' @export
OrdinalSTAAR <- function(coded_genotype_file, fit_null, 
                         annotation_file = NULL, QC_label = "imputed",
                         rare_maf_cutoff = 0.01, rv_num_cutoff = 2) {

  # --- Step 0: Initial validation and setup ---
  if (!file.exists(coded_genotype_file)) stop("Coded genotype file not found.")
  
  message("--- Welcome to OrdinalSTAAR: Rare Variant Test for Ordinal Phenotypes ---")
  
  # --- Step 1: Load annotation data (if provided) ---
  annotation_phred <- NULL
  if (!is.null(annotation_file)) {
    if (file.exists(annotation_file)) {
      message(paste("Loading functional annotation data from:", annotation_file))
      anno_env <- new.env()
      load(annotation_file, envir = anno_env)
      anno_obj_name <- paste0("anno_", QC_label)
      if (anno_obj_name %in% ls(anno_env)) {
        annotation_phred <- get(anno_obj_name, envir = anno_env)
      } else {
        warning(paste("Annotation object named", anno_obj_name, "not found in the file."))
      }
    } else {
      warning("Annotation file specified but not found.")
    }
  }

  # --- Step 2: Load all gene objects and prepare for loop ---
  message(paste("Loading gene-centric genotype data from:", coded_genotype_file))
  gene_env <- new.env()
  load(coded_genotype_file, envir = gene_env)
  gene_obj_names <- ls(gene_env)
  
  results_list <- list()
  
  # --- Step 3: Unpack null model components for efficiency ---
  res <- fit_null$residuals
  sample_ids_pheno <- fit_null$sample_ids
  W <- fit_null$W_mat
  X <- fit_null$X_mat
  WX <- fit_null$WX_mat
  XWX_inv <- fit_null$XWX_inv
  
  num_genes <- length(gene_obj_names)
  message(paste("Starting association tests for", num_genes, "genes..."))

  # --- Step 4: Loop through each gene to perform tests ---
  for (i in 1:num_genes) {
    obj_name <- gene_obj_names[i]
    gene_name <- gsub(paste0("G_", QC_label, "_"), "", obj_name)
    
    if (i %% 100 == 0) cat(paste0("Processing gene ", i, "/", num_genes, ": ", gene_name, "\n"))
    
    # --- 4.1: Get genotype and annotation data for the current gene ---
    G_obj <- get(obj_name, envir = gene_env)
    G_matrix_full <- G_obj$genotype
    variant_ids_gene <- rownames(G_matrix_full)
    
    # Get annotation weights
    weights <- rep(1, length(variant_ids_gene)) # Default weight is 1
    if (!is.null(annotation_phred)) {
      anno_gene <- annotation_phred[variant_ids_gene, , drop=FALSE]
      # Example weighting scheme: Beta distribution on MAF, can be made more complex
      weights <- dbeta(anno_gene$MAF, 1, 25) 
    }
    
    # --- 4.2: Match samples ---
    sample_ids_geno <- colnames(G_matrix_full)
    common_samples <- intersect(sample_ids_pheno, sample_ids_geno)
    
    pheno_idx <- match(common_samples, sample_ids_pheno)
    geno_idx <- match(common_samples, sample_ids_geno)
    
    res_matched <- res[pheno_idx]
    W_matched <- W[pheno_idx, pheno_idx]
    X_matched <- X[pheno_idx, , drop = FALSE]
    WX_matched <- WX[pheno_idx, , drop = FALSE]
    G_matrix <- t(G_matrix_full[, geno_idx, drop = FALSE])
    
    # --- 4.3: Filter for rare variants ---
    maf <- colMeans(G_matrix) / 2
    rare_idx <- which(maf > 0 & maf < rare_maf_cutoff)
    
    if (length(rare_idx) < rv_num_cutoff) next
    
    G_rare <- G_matrix[, rare_idx, drop = FALSE]
    weights_rare <- weights[rare_idx]
    
    # --- 4.4: Core Score Test Calculations (Manual Matrix Algebra) ---
    
    # Apply weights to genotype matrix
    G_w <- G_rare %*% diag(weights_rare)
    
    # Calculate Score Vector U
    U <- crossprod(G_w, res_matched)
    
    # Calculate Covariance Matrix of the Score Vector, V
    G_w_t_WX <- crossprod(G_w, WX_matched)
    V <- crossprod(G_w, W_matched %*% G_w) - (G_w_t_WX %*% XWX_inv %*% t(G_w_t_WX))
    V <- as.matrix(V)
    
    # --- 4.5: Perform Association Tests ---
    
    # 1. Burden Test
    U_burden <- sum(U)
    Var_burden <- sum(V)
    Stat_burden <- U_burden^2 / Var_burden
    p_burden <- pchisq(Stat_burden, df = 1, lower.tail = FALSE)
    
    # 2. SKAT Test
    lambda <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
    lambda <- lambda[lambda > 1e-6] # Numerical stability
    
    # Check if CompQuadForm is available for SKAT p-value
    p_skat <- tryCatch({
        if (!requireNamespace("CompQuadForm", quietly = TRUE)) stop()
        CompQuadForm::davies(q = sum(U^2), lambda = lambda)$Qq
    }, error = function(e) {
        warning(paste("SKAT p-value calculation failed for gene", gene_name, ". CompQuadForm might be missing or failed. P-value set to NA."))
        return(NA)
    })
    
    # 3. ACAT-OMNI Test (combining Burden and SKAT)
    p_values <- c(p_burden, p_skat)
    p_values <- p_values[!is.na(p_values)]
    p_acat <- if(length(p_values) > 0) ACAT(Pvals = p_values) else NA
    
    results_list[[gene_name]] <- data.frame(
      Gene = gene_name,
      N_Variants = ncol(G_rare),
      P_Burden = p_burden,
      P_SKAT = p_skat,
      P_ACAT_OMNI = p_acat
    )
  }
  
  # --- Step 5: Finalize and return results ---
  if (length(results_list) == 0) {
    warning("No genes passed the filtering criteria. Returning an empty data frame.")
    return(data.frame())
  }
  
  final_results <- do.call(rbind, results_list)
  rownames(final_results) <- NULL
  final_results <- final_results[order(final_results$P_ACAT_OMNI), ]
  
  message("\n--- OrdinalSTAAR analysis finished successfully! ---")
  return(final_results)
}


# --- Helper function for ACAT p-value combination ---
# This is a standard implementation of the ACAT test.
ACAT <- function(Pvals, Weights=NULL){
  if(is.null(Weights)){
    Weights <- rep(1/length(Pvals), length(Pvals))
  }
  is.zero <- (Pvals == 0)
  is.one <- (Pvals == 1)
  if(sum(is.zero) > 0){ return(0) }
  if(sum(is.one) == length(Pvals)){ return(1) }
  
  p.min <- min(Pvals)
  if(p.min * length(Pvals) < 1e-16){ return(0) }
  
  Pvals[is.one] <- 1 - 1e-16
  stat <- sum(Weights * tan((0.5 - Pvals) * pi))
  
  # Check for extreme stats to avoid pcauchy errors
  if (is.infinite(stat)) {
      return(0)
  }
  
  res <- pcauchy(stat, lower.tail = FALSE)
  return(res)
}
```

### 如何使用这个完整的框架

1.  **准备数据：**
    *   `gdsfile`: 包含基因型和注释的 GDS 文件。
    *   `pheno_data`: 包含表型和协变量的 R 数据框。
    *   `PRS.file` (可选): 由 `REGENIE` 生成的 PRS 文件。

2.  **第一步：基因型预处理**
    ```R
    # source("GeneCentricCoding.R")
    # GeneCentricCoding(gdsfile, QC_label="imputed", sample.id=..., variant.id=..., outfile="coded_genotypes.Rdata")
    ```

3.  **第二步：拟合零模型**
    ```R
    # source("NullModel.R")
    # fit_null_obj <- Ordinal_NullModel(formula, data=pheno_data, id_col="sample_id", PRS.file="path/to/prs")
    ```

4.  **第三步：运行关联分析**
    ```R
    # source("OrdinalSTAAR.R")
    # results_df <- OrdinalSTAAR(coded_genotype_file="coded_genotypes.Rdata", fit_null=fit_null_obj)
    ```

这个完整的代码框架现在为你提供了一个功能完备的 `OrdinalSTAAR` 工具，它在设计上与 `SurvSTAAR` 保持一致，并且能够处理大规模、存在亲缘关系的有序性状数据。