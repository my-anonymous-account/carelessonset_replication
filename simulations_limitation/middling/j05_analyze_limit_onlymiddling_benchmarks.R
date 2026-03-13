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

# the benchmark methods
methods <- c("reliability", "antonym", "longstring")
cutoffs <- c(reliability = 0.3, antonym = -0.03, longstring = 7)

for(prop_crls in prop_crls_arr)
{
  
  load(paste0(dir_load, "/", suffix_type, "_", suffix_location, "_prop-crls=", prop_crls, "_n=", n, "_R=", R, ".Rdata"))

  classif_perf <- matrix(NA_real_, length(methods), 4L)
  rownames(classif_perf) <- methods
  colnames(classif_perf) <- c("TPR", "FNR", "TNR", "FPR")
  
  for(method in methods)
  {
    ## initialize
    classif_mat <- matrix(NA_real_, R, 4L)
    colnames(classif_mat) <- c("TPR", "FNR", "TNR", "FPR")
    
    ## get scores and cutoff of the method in question
    scores_mat <- get(paste0(method, "_mat"))
    cutoff <- cutoffs[method]
    
    
    for(r in seq_len(R))
    {
      flagged <- which(scores_mat[,r] > cutoff)
      classif_mat[r,] <- classification(idx_flagged = flagged, idx_careless = idx_careless, n = n)
    } # FOR r
    
    ## collect results
    classif_perf[method, ] <- colMeans(classif_mat)
    
  } # FOR method
  
  
  ## save the matrix object
  save(classif_perf,
       file = paste0(dir_save, "/classification_", suffix_type, "_", suffix_location, 
                     "_prop-crls=", prop_crls, "_benchmarks", "_n=", n,
                     "_R=", R, ".Rdata"))
} # FOR prop_crls
