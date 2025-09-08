```R
SurvSTAAR_O = function(Geno, objNull, annotation_rank = NULL, MAC = NULL,
                       use_SPA = NULL, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                       weight_A, weight_B, weight_S,
                       combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 20, verbose = FALSE) {

  if (!inherits(Geno, "matrix") && !inherits(Geno, "Matrix")) stop("Genotype is not a matrix")

  if (ncol(Geno) != nrow(weight_A) | ncol(Geno) != nrow(weight_B) | ncol(Geno) != nrow(weight_S)) stop("Dimensions don't match for genotype and weights")

  if (is.null(use_SPA)) use_SPA = objNull$use_SPA

  if (is.null(MAC)) MAC = colSums(Geno)


  ### individual test
  single_test = exactScore(objNull = objNull, G_mat = Geno, use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff, verbose = verbose)
  Score = single_test$Score
  Covariance = single_test$Covariance
  if(use_SPA) {
    Pvalue = single_test$result$Pvalue_SPA
  } else {
    Pvalue = single_test$result$Pvalue
  }


  ### set based test
  if (verbose) print(paste0("SKAT test:   begin at ", Sys.time()))
  Pvalue_S = SKAT(Geno, Score, Covariance, Pvalue, MAC, weight_S, weight_B, objNull, use_SPA, SPA_filter, SPA_filter_cutoff, combine_ultra_rare, ultra_rare_mac_cutoff, verbose)

  if (verbose) print(paste0("ACAT test:   begin at ", Sys.time()))
  Pvalue_A = ACAT(Geno, Score, Covariance, Pvalue, MAC, weight_A, weight_B, objNull, use_SPA, SPA_filter, SPA_filter_cutoff, combine_ultra_rare, ultra_rare_mac_cutoff, verbose)

  if (verbose) print(paste0("Burden test: begin at ", Sys.time()))
  Pvalue_B = Burden(Geno, Score, Covariance, weight_B, objNull, use_SPA, SPA_filter, SPA_filter_cutoff, verbose)


  ### combine results
  results_pvalue = rbind(Pvalue_S, Pvalue_A, Pvalue_B)
  rownames(results_pvalue) = c("SKAT-(1,25)", "SKAT-(1,1)", "ACAT-(1,25)", "ACAT-(1,1)", "Burden-(1,25)", "Burden-(1,1)")

  if (!is.null(annotation_rank)) {
    colnames(results_pvalue) = c("Beta", colnames(annotation_rank))
  } else {
    colnames(results_pvalue) = "results"
  }


  results_STAAR_O = CCT(c(results_pvalue))
  results_ACAT_O = CCT(c(results_pvalue[, 1]))
  results_STAAR_S = CCT(results_pvalue[1:2, ])
  results_STAAR_A = CCT(results_pvalue[3:4, ])
  results_STAAR_B = CCT(results_pvalue[5:6, ])

  STAAR_S_1_25 = CCT(Pvalue_S[1, ])
  STAAR_S_1_1  = CCT(Pvalue_S[2, ])
  STAAR_A_1_25 = CCT(Pvalue_A[1, ])
  STAAR_A_1_1  = CCT(Pvalue_A[2, ])
  STAAR_B_1_25 = CCT(Pvalue_B[1, ])
  STAAR_B_1_1  = CCT(Pvalue_B[2, ])

  results_SurvSTAAR = c(STAAR_S_1_25, STAAR_S_1_1, STAAR_A_1_25, STAAR_A_1_1, STAAR_B_1_25, STAAR_B_1_1)
  results_pvalue = cbind(results_pvalue, results_SurvSTAAR)


  results = list("SurvSTAAR_O" = results_STAAR_O,
                 "ACAT_O" = results_ACAT_O,
                 "SurvSTAAR_SKAT" = results_STAAR_S,
                 "SurvSTAAR_ACAT" = results_STAAR_A,
                 "SurvSTAAR_Burden" = results_STAAR_B,
                 "SurvSTAAR_pvalue" = results_pvalue)


  return(results)
}

```