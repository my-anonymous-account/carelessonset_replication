rm(list = ls()) ; cat("\014")
source("R/segmentation_performance_measures.R")

dir_load <- "simulations_limitation/shortsurvey/results"
dir_save <- "simulations_limitation/shortsurvey/analyzed"

suffix_type <- "j05sim-limit-shortsurvey"
nums_dropped <- 1:4 # number of traits dropped
suffix_location <- c("mixed")
R <- 100L 
n <- 500L

# specify proportion of careless respondents
prop_crls_arr <- c(0.02, 0.04, seq(from = 0.1, to = 0.6, by = 0.1), 0.8, 1.0)

types <- c("random", "straightlining", "extreme", "pattern", "careless", "careful")
methods <- c("CP_LS_SC", "CP_LS", "CP_SC")

# initialize
RESULTS <- NULL


for(num_dropped in nums_dropped)
{
  
  ### LOOP OVER SCENARIOS WITH NONZEO CARELESSNESS
  for(prop_crls in prop_crls_arr)
  {
    
    load(paste0(dir_load, "/", suffix_type, 
                "_num-traits-dropped=", num_dropped, "_", suffix_location, 
                "_prop-crls=", prop_crls, "_n=", n, "_R=", R, ".Rdata"))
    
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
        div_random <- div_straight <- div_extreme <- div_pattern <- div_careless
        div <- matrix(NA_real_, length(types) - 1L, 5L)
        rownames(div) <- types[-which(types == "careful")]
        colnames(div) <- colnames(div_careless)
        
        
        for(r in 1:R)
        {
          flagged <- out_mod[[r]][[a]]
          
          ## divergence performances
          div_careless[r,] <- onset_divergence(onset_estimated = flagged,
                                               onset_true = onset_true,
                                               num_items = num_items)
          
          div_extreme[r,] <- onset_divergence(onset_estimated = flagged[idx_extreme],
                                              onset_true = onset_true[idx_extreme],
                                              num_items = num_items)
          
          div_random[r,] <- onset_divergence(onset_estimated = flagged[idx_random],
                                             onset_true = onset_true[idx_random],
                                             num_items = num_items)
          
          div_pattern[r,] <- onset_divergence(onset_estimated = flagged[idx_pattern],
                                              onset_true = onset_true[idx_pattern],
                                              num_items = num_items)
          
          div_straight[r,] <- onset_divergence(onset_estimated = flagged[idx_straight],
                                               onset_true = onset_true[idx_straight],
                                               num_items = num_items)
          
          
          ## classification performance
          classif_mat[r,] <- classification(idx_flagged = which(!is.na(flagged)), 
                                            idx_careless = idx_careless, n = n)
          
        } # FOR r
        
        ## average across the repetitions
        classif_perf_method <- colMeans(classif_mat)
        
        
        ## collect divergence results
        # NB: if no changepoint is ever flagged (rarely happens), then the divergences will be NaN, so we omit them in the averaging...
        # ... and the corresponding poor detection performance will be captured by the FNR performance measures
        div["random",] <- colMeans(div_random, na.rm = TRUE)
        div["straightlining",] <- colMeans(div_straight, na.rm = TRUE)
        div["extreme",] <- colMeans(div_extreme, na.rm = TRUE)
        div["pattern",] <- colMeans(div_pattern, na.rm = TRUE)
        div["careless",] <- colMeans(div_careless, na.rm = TRUE)
        div_perf_method <- div
        
        ## gather relevant results
        df <- as.data.frame(div_perf_method[,c("absdiff_mean", "fnr")])
        df$type <- rownames(div_perf_method)
        df$fpr <- NA_real_ # initialize
        
        ## append with 'attentive' type (no FNR or MAE can be computed here)
        attntve <- data.frame(absdiff_mean = NA_real_, 
                              fnr = NA_real_,
                              type = "attentive", 
                              fpr = classif_perf_method["FPR"])
        df <- rbind(df, attntve)
        rownames(df) <- NULL
        
        ## append columns for methods, alpha, contamination, numitems
        df$method <- method
        df$prop_crls <- prop_crls
        df$alpha <- alpha[a]
        df$numitems <- 300L - num_dropped * 60L
        
        ## append to RESULTS
        RESULTS <- rbind(RESULTS, df)
        
      } # FOR method
    } # FOR a in seq_along(alpha)
  } # FOR prop_crlss
} # FOR num_dropped



### NOW LOOP OVER ZERO-CARELESSNESS SETTINGS
for(num_dropped in nums_dropped)
{
  # load results without contamination
  load(paste0(dir_load, "/", suffix_type, "_num-traits-dropped=", num_dropped, 
              "_prop-crls=", 0.0, "_n=", n,  "_R=", R, ".Rdata"))
  
  
  for(a in seq_along(alpha))
  {
    for(method in methods){
      
      perf <- matrix(NA_real_, 1, 3)
      rownames(perf) <- "careful"
      colnames(perf) <- c(paste0("m=", 0:1), "ARI")
      
      out_mod <- get(paste0("out_", method))
      
      perf_careful <-  matrix(NA_real_, R, 3L)
      cp_arr <- array(NA_real_, dim = c(n, num_items, R))
      
      classif_mat <- matrix(NA_real_, R, 4L)
      colnames(classif_mat) <- c("TPR", "FNR", "TNR", "FPR")
      
      for(r in 1:R)
      {
        flagged <- out_mod[[r]][[a]]
        classif_mat[r,] <- classification(idx_flagged = which(!is.na(flagged)), 
                                          idx_careless = integer(), n = n)
        
      } # FOR r
      
      classif_perf_method <- colMeans(classif_mat)
      
      ## append with 'attentive' type (no FNR or MAE can be computed here)
      attntve <- data.frame(absdiff_mean = NA_real_, 
                            fnr = NA_real_,
                            type = "attentive", 
                            fpr = classif_perf_method["FPR"])
      rownames(attntve) <- NULL
      
      ## append columns for methods, alpha, contamination, numitems
      attntve$method <- method
      attntve$prop_crls <- 0.0
      attntve$alpha <- alpha[a]
      attntve$numitems <- 300L - num_dropped * 60L
      
      ## append to RESULTS
      RESULTS <- rbind(RESULTS, attntve)
      
    } # FOR method
  } # FOR a in seq_along(alpha)
} # FOR num-_dropped



## update column names to be more meaningful
colnames(RESULTS) <- c("MAE", "FNR", "type", "FPR", "method", "prop_crls", "sig_level", "numitems")

# convert to long format
RESULTS_long <- reshape2::melt(data = RESULTS, 
                              id.vars = c("method", "type", "prop_crls", "numitems", "sig_level"), 
                              variable.name = "measure", 
                              value.name = "value")

## change the factor levels
library("dplyr")
RESULTS_out <- RESULTS_long

method_labels <- c(CP_LS = "LSP only", 
                   CP_SC = "RE only",
                   CP_LS_SC = "Both\ndimensions")
                   
measure_labels <- c(FNR = "False~negative~rate",
                    FPR = "False~positive~rate",
                    MAE = "Mean~absolute~error")
type_levels <- c(attentive =  "Attentive",
                 careless = "Careless", 
                 random = "Random", 
                 extreme = "Extreme", 
                 straightlining = "Straightlining", 
                 pattern = "Pattern")

RESULTS_out$method <- factor(recode(RESULTS_long$method,
                    CP_LS_SC = method_labels["CP_LS_SC"],
                    CP_SC = method_labels["CP_SC"],
                    CP_LS = method_labels["CP_LS"],
                    ), levels = method_labels)

RESULTS_out$measure <- factor(
  recode(RESULTS_long$measure,
         MAE = measure_labels["MAE"],
         FNR = measure_labels["FNR"],
         FPR = measure_labels["FPR"]),
  levels = measure_labels)

RESULTS_out$type <- factor(
  recode(RESULTS_long$type,
         attentive = type_levels["attentive"],
         careless = type_levels["careless"],
         random = type_levels["random"],
         extreme = type_levels["extreme"],
         straightlining = type_levels["straightlining"],
         pattern = type_levels["pattern"]),
  levels = type_levels)


RESULTS_out$numitems <-
  paste0('p*" = "*', as.character(RESULTS_long$numitems)) # to-be-parsed labels
RESULTS_out$numitems <- factor(RESULTS_out$numitems, levels = unique(RESULTS_out$numitems))

RESULTS_out$sig_level <- 
  paste0(as.character(100 * RESULTS_long$sig_level), "%") # to-be-parsed labels

## update naming
RESULTS <- RESULTS_out

## save results
save(RESULTS, 
     file = paste0(dir_save, "/gatheredresults_CODERS_", 
                   suffix_type, "_n=", n, "_R=", R, ".Rdata"))
