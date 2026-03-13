rm(list = ls()) ; cat("\014")
source("R/segmentation_performance_measures.R")

dir_load <- "simulations_limitation/middling/results"
dir_save <- "simulations_limitation/middling/analyzed"

suffix_type <- "j05sim-limit-middlingonly"
suffix_location <- "mixed"
R <- 100L 
n <- 500L
num_items <- 300L
plots <- FALSE

# specify proportion of careless respondents
prop_crls_arr <- c(0.02, 0.04, seq(from = 0.1, to = 0.6, by = 0.1), 0.8, 1.0)

types <- c("middle")
methods <- c("CP_LS_SC", "CP_LS", "CP_SC")


for(prop_crls in prop_crls_arr)
{
  
  load(paste0(dir_load, "/", suffix_type, "_", suffix_location, "_prop-crls=", prop_crls, "_n=", n, "_R=", R, ".Rdata"))
  #load(paste0(dir_load, "/", suffix_type, "_", suffix_location, "_prop-crls=", prop_crls, "_R=", R, ".Rdata"))
  
  ari_perf <- matrix(NA_real_, length(types), length(methods))
  rownames(ari_perf) <- types
  colnames(ari_perf) <- methods
  
  # true onset
  onset_true <- careless_start
  onset_true[idx_careful] <- NA_integer_ # set onset location of truly careful respondents to NA (expected by fun below)
  
  
  for(a in seq_along(alpha))
  {
    
    classif_perf <- matrix(NA_real_, length(methods), 4L)
    rownames(classif_perf) <- methods
    colnames(classif_perf) <- c("TPR", "FNR", "TNR", "FPR")
    div_perf <- rep(list(NA), length(methods))
    names(div_perf) <- methods
    
    for(method in methods){
      
      perf <- matrix(NA_real_, length(types), 3)
      rownames(perf) <- types
      colnames(perf) <- c(paste0("m=", 0:1), "ARI")
      
      out_mod <- get(paste0("out_", method))
      
      perf_random <- perf_straight <- perf_extreme <- perf_middle <- 
        perf_pattern <- perf_careless <- perf_careful <-  matrix(NA_real_, R, 3L)
      cp_arr <- array(NA_real_, dim = c(n, num_items, R))
      
      classif_mat <- matrix(NA_real_, R, 4L)
      colnames(classif_mat) <- c("TPR", "FNR", "TNR", "FPR")
      
      # divergences
      div_careless <- matrix(NA_real_, R, 5L)
      colnames(div_careless) <- c("absdiff_mean", "absdiff_median", "reldiff_mean", "reldiff_median", "fnr")
      div <- matrix(NA_real_, length(types), 5L)
      colnames(div) <- colnames(div_careless)
      rownames(div) <- "careless"
      
      for(r in 1:R)
      {
        flagged <- out_mod[[r]][[a]]
        # cp_mat      <- get_cp_mat(flagged = flagged, num_items = num_items)
        # cp_arr[,,r] <- cp_mat
        # 
        # perf_middle[r,] <- colMeans(
        #   evaluate_endsurvey(cp_mat[idx_middle,,drop = FALSE], 
        #                      start_careless = careless_start[idx_middle]))
        # 
        # if(prop_crls < 1.0){
        #   perf_careful[r,] <- colMeans(
        #     evaluate_endsurvey(cp_mat[idx_careful,,drop = FALSE], 
        #                        start_careless = NA, 
        #                        careful = TRUE))
        # } # IF
        
        ## divergence performances
        div_careless[r,] <- onset_divergence(onset_estimated = flagged,
                                            onset_true = onset_true,
                                            num_items = num_items)
        
        
        ## classification performance
        classif_mat[r,] <- classification(idx_flagged = which(!is.na(flagged)), 
                                          idx_careless = idx_careless, n = n)
        
      } # FOR r
      
      # average across the repetitions
      #perf["middle",] <- colMeans(perf_middle)
      classif_perf[method, ] <- colMeans(classif_mat)
      
      # save as csv
      # write.csv(round(perf,3), 
      #           file = paste0(dir_save, "/RES_sim_base_singleCP_",
      #                         method, "_prop-crls", prop_crls, "_", alpha[a], ".csv"))
      
      # collect ARI results
      #ari_perf[, method] <- perf[, "ARI"]
      
      ## collect divergence results
      # NB: if no changepoint is ever flagged (rarely happens), then the divergences will be NaN, so we omit them in the averaging...
      # ... and the corresponding poor detection performance will be captured by the FNR performance measures
      div["careless",] <- colMeans(div_careless, na.rm = TRUE)
      div_perf[[method]] <- div
      
      
    } # FOR method
    
    save(classif_perf,
         file = paste0(dir_save, "/classification_", suffix_type, "_", suffix_location, 
                       "_prop-crls=", prop_crls, "_alpha=", alpha[a], "_n=", n,
                       "_R=", R, ".Rdata"))
    
    save(div_perf,
         file = paste0(dir_save, "/divergence_", suffix_type, "_", suffix_location, 
                       "_prop-crls=", prop_crls, "_alpha=", alpha[a], "_n=", n,
                       "_R=", R, ".Rdata"))
    
  } # FOR a in seq_along(alpha)
} # FOR prop_crlss
