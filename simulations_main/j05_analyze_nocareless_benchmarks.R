rm(list = ls()) ; cat("\014")
source("R/segmentation_performance_measures.R")

dir_load <- "simulations_main/results"
dir_save <- "simulations_main/analyzed"

suffix_type <- "j05sim-main"
R <- 100L 
n <- 500L
num_items <- 300L

# the benchmark methods
methods <- c("reliability", "antonym", "longstring")
cutoffs <- c(reliability = 0.3, antonym = -0.03, longstring = 7) # johnson's cutoff (except LS): careless if exceeded

load(paste0(dir_load, "/", suffix_type, "_prop-crls=", 0.0, "_n=", n, "_R=", R, ".Rdata"))

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
    flagged <- which(scores_mat[,r] > cutoff) # NOTE: some scores might be NA. We ignore them here, which is equivalent to labeling them as attentive
    classif_mat[r,] <- classification(idx_flagged = flagged, idx_careless = integer(), n = n)
  } # FOR r
  
  ## collect results
  classif_perf[method, ] <- colMeans(classif_mat)
  
} # FOR method


## save the matrix object
save(classif_perf,
     file = paste0(dir_save, "/classification_", suffix_type, 
                   "_prop-crls=", 0.0, "_benchmarks", "_n=", n,
                   "_R=", R, ".Rdata"))
