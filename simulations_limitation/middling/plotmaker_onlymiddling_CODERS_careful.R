rm(list = ls()) ; cat("\014")
source("R/segmentation_performance_measures.R")

dir_load <- "simulations_limitation/middling/results"
dir_save <- "simulations_limitation/middling/plots"
suffix_type <- "j05sim-limit-middlingonly"
suffix_locations <- "mixed"
types <- "middling"
R <- 100L 
n <- 500L
num_items <- 300L

# specify proportion of careless respondents
prop_crls_arr <- c(0.02, 0.04, seq(from = 0.1, to = 0.6, by = 0.1), 0.8, 1.0)
methods <- c("CP_LS_SC", "CP_LS", "CP_SC")

# initialize
RESULTS <- NULL


for(suffix_location in suffix_locations)
{
  
  ### LOOP OVER SCENARIOS WITH NONZEO CARELESSNESS
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
        
        
        for(r in 1:R)
        {
          flagged <- out_mod[[r]][[a]]
          
        
          ## classification performance
          classif_mat[r,] <- classification(idx_flagged = which(!is.na(flagged)), 
                                            idx_careless = idx_careless, n = n)
          
        } # FOR r
        
        ## average across the repetitions
        classif_perf_method <- colMeans(classif_mat)
        
        
        ## append with 'attentive' type (no FNR or MAE can be computed here)
        attntve <- data.frame(absdiff_mean = NA_real_, 
                              fnr = NA_real_,
                              type = "attentive", 
                              fpr = classif_perf_method["FPR"])
        df <- attntve
        rownames(df) <- NULL
        
        ## append columns for methods, alpha, contamination, location
        df$method <- method
        df$prop_crls <- prop_crls
        df$alpha <- alpha[a]
        df$onset <- suffix_location
        
        ## append to RESULTS
        RESULTS <- rbind(RESULTS, df)
        
      } # FOR method
    } # FOR a in seq_along(alpha)
  } # FOR prop_crlss
} # FOR suffix_location



### NOW LOOP OVER ZERO-CARELESSNESS SETTINGS

## load data
load(paste0("simulations_main/results/j05sim-main", "_prop-crls=", 0.0, "_n=", n, "_R=", R, ".Rdata"))

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
    
    ## append columns for methods, alpha, contamination, location
    attntve$method <- method
    attntve$prop_crls <- 0.0
    attntve$alpha <- alpha[a]
    
    ## NB: there is no carelessness onset point here (no carelessness at all),
    # but it is easier for later to have the no-careless results thrice (once for
    # each regime)
    df <- NULL
    
    for(suffix_location in suffix_locations)
    {
      tmp <- attntve
      tmp$onset <- suffix_location
      df <- rbind(tmp, df)
    }
    
    ## append to RESULTS
    RESULTS <- rbind(RESULTS, df)
    
  } # FOR method
} # FOR a in seq_along(alpha)


## update column names to be more meaningful
colnames(RESULTS) <- c("MAE", "FNR", "type", "FPR", "method", "prop_crls", "sig_level", "onset")

## reduce data frame to the relevant columns
RESULTS <- RESULTS[, c("FPR", "method", "prop_crls", "sig_level")]

# convert to long format
RESULTS_long <- reshape2::melt(data = RESULTS, 
                              id.vars = c("method", "prop_crls", "sig_level"), 
                              variable.name = "measure", 
                              value.name = "value")

## change the factor levels
library("dplyr")
RESULTS_out <- RESULTS_long

method_labels <- c(CP_LS = "LSP only", 
                   CP_SC = "RE only",
                   CP_LS_SC = "Both\ndimensions")
                   
measure_labels <- c(FPR = "False~positive~rate")


RESULTS_out$method <- factor(recode(RESULTS_long$method,
                    CP_LS_SC = method_labels["CP_LS_SC"],
                    CP_SC = method_labels["CP_SC"],
                    CP_LS = method_labels["CP_LS"],
                    ), levels = method_labels)


RESULTS_out$sig_level <- 
  paste0(as.character(100 * RESULTS_long$sig_level), "%") # to-be-parsed labels

## update naming
DF <- RESULTS_out

## drop uninformative cases with 100% carelessness (no FPR possible here)
DF <- DF %>% filter(prop_crls < 1)

## plot FPR
P <- ggplot(data = DF, mapping = aes(x = prop_crls, y = value, color = method, linetype = sig_level)) +
  geom_line() +
  ylim(0,1) +
  scale_linetype_manual("Significance\nlevel", values = c("dashed", "solid", "twodash")) +
  scale_alpha_manual("Significance\nlevel", values = c(0.4, 1, 0.4)) + 
  scale_color_manual("\n\nCODERS", values = c("#619CFF", "#F8766D", "black")) + # ugly hack to avoid top cutoff of legend
  xlab("Prevalence of partially careless respondents") + 
  ylab("False positive rate") +
  theme_bw() +
  theme(legend.position = "right", legend.key.width = unit(19, "pt"), 
        # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1)
  ) +
  guides(color = guide_legend(nrow = 3, reverse = TRUE, order = 1))

# save the plot
ggsave(filename = paste0("plot-attentive_", suffix_type, "_", suffix_location, "_n=", n, "_R=", R, ".pdf"),
       plot = P,
       path = dir_save, 
       device = "pdf",
       width = 7.5,
       height = 2.5)
