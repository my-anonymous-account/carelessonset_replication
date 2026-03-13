rm(list = ls()) ; cat("\014")
source("R/segmentation_performance_measures.R")

dir_load <- "simulations_additional/multiple_changepoints/results"
dir_save <- "simulations_additional/multiple_changepoints/analyzed"

suffix_type <- "j05sim-multipleCP"
suffix_location <- "early"
R <- 100L 
n <- 500L
num_items <- 300L
plots <- FALSE

# specify proportion of careless respondents
prop_crls_arr <- c(0.02, 0.04, seq(from = 0.1, to = 0.6, by = 0.1), 0.8, 1.0)

# the benchmark methods
methods <- c("reliability", "antonym", "longstring")
cutoffs <- c(reliability = 0.3, antonym = -0.03, longstring = 7) # johnson's cutoff (except LS): careless if exceeded

# schroeders paper: "At least six identical responses in a row (before reverse coding negatively phrased items) were considered conspicuous for the longstring index (see also Niessen et al., 2016)."
# Niessen says: ". The scree plot did not show a clear cutoff, so we were conservative (a cutoff value between five and seven seemed reasonable based on the plot) and chose a cutoff of more than seven consecutive responses"  They also investigat the 6 cutoff, which they say is based on data. So I'd just go with 6 and not look back (shouldn't matter much for us anyway)
# The original johnson paper differentiated cutoffs for different response options (and devised cutoff via scree plot). Meade and craig don't do that and explicitly say so: 
# "While Johnson computed this index for each individual response option (e.g., the maximum number of consec- utive responses of 1, 2, etc.), we simply computed the maximum number of items with consecutive response, regardless of the value of that response"



# the subgroups
subgroups <- c("careless", "straightlining", "extreme", "pattern", "random")


for(prop_crls in prop_crls_arr)
{
  
  load(paste0(dir_load, "/", suffix_type, "_", 
              suffix_location, 
              "_prop-crls=", prop_crls, "_n=", n, "_R=", R, ".Rdata"))

  ## classification 
  classif_perf <- matrix(NA_real_, length(methods), 4L)
  rownames(classif_perf) <- methods
  colnames(classif_perf) <- c("TPR", "FNR", "TNR", "FPR")

  
  ## FNR by subgroup
  fnr_perf <- matrix(NA_real_, length(methods), length(subgroups))
  rownames(fnr_perf) <- methods
  colnames(fnr_perf) <- subgroups
  
  for(method in methods)
  {
    ## initialize
    classif_mat <- matrix(NA_real_, R, 4L)
    colnames(classif_mat) <- c("TPR", "FNR", "TNR", "FPR")
    subgroups_fnr_mat <- matrix(NA_real_, R, 5L)
    colnames(subgroups_fnr_mat) <- subgroups
    
    ## get scores and cutoff of the method in question
    scores_mat <- get(paste0(method, "_mat"))
    cutoff <- cutoffs[method]
    
    
    for(r in seq_len(R))
    {
      idx_flagged <- which(scores_mat[,r] > cutoff) # NOTE: some scores might be NA. We ignore them here, which is equivalent to labeling them as attentive
      
      ## classification matrix
      classif_mat[r,] <- classification(idx_flagged = idx_flagged, idx_careless = idx_careless, n = n)
      
      subgroups_fnr_mat[r, "careless"] <- 
        benchmarks_fnr(idx_flagged = idx_flagged, 
                       idx_positives = idx_careless,
                       n = n)
      subgroups_fnr_mat[r, "straightlining"] <- 
        benchmarks_fnr(idx_flagged = intersect(idx_flagged, idx_straight), 
                       idx_positives = intersect(idx_careless, idx_straight),
                       n = n)
      subgroups_fnr_mat[r, "extreme"] <- 
        benchmarks_fnr(idx_flagged = intersect(idx_flagged, idx_extreme), 
                       idx_positives = intersect(idx_careless, idx_extreme),
                       n = n)
      subgroups_fnr_mat[r, "pattern"] <- 
        benchmarks_fnr(idx_flagged = intersect(idx_flagged, idx_pattern), 
                       idx_positives = intersect(idx_careless, idx_pattern),
                       n = n)
      subgroups_fnr_mat[r, "random"] <- 
        benchmarks_fnr(idx_flagged = intersect(idx_flagged, idx_random), 
                       idx_positives = intersect(idx_careless, idx_random),
                       n = n)
    } # FOR r
    
    ## collect results
    classif_perf[method, ] <- colMeans(classif_mat)
    fnr_perf[method, ] <- colMeans(subgroups_fnr_mat)
    
  } # FOR method
  
  
  ## save the matrix object
  save(classif_perf, fnr_perf,
       file = paste0(dir_save, "/classification_", suffix_type, "_", suffix_location, 
                     "_prop-crls=", prop_crls, "_benchmarks", "_n=", n,
                     "_R=", R, ".Rdata"))
} # FOR prop_crls
