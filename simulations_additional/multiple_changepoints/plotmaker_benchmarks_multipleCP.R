rm(list = ls()) ; cat("\014")

library("dplyr")
library("ggplot2")
library("tidyr")

dir_load <- "simulations_additional/multiple_changepoints/analyzed"
dir_save <- "simulations_additional/multiple_changepoints/plots"
suffix_type <- "j05sim-multipleCP"
suffix_locations <- "early"

## contamination levels = proportion of respondents who become careless
prop_crlss_arr <- c(0.02, 0.04, seq(from = 0.1, to = 0.6, by = 0.1), 0.8, 1.0)
R <- 100L
n <- 500L

## initialize
df_all <- NULL

for(suffix_location in suffix_locations)
{
  
  df <- NULL 
  
  for(prop_crls in c(0.0, prop_crlss_arr))
  {
    # load classification results
    if(prop_crls > 0.0)
    {
      load(paste0(dir_load, "/classification_", suffix_type, "_", suffix_location, 
                  "_prop-crls=", prop_crls, "_benchmarks", "_n=", n,
                  "_R=", R, ".Rdata"))
    } else
    {
      # load classification performance without contamination
      load(paste0("simulations_main/analyzed", "/classification_", "j05sim-main", 
                  "_prop-crls=", 0.0, "_benchmarks", "_n=", n,
                  "_R=", R, ".Rdata"))
    }
    
    
    # prepare the data frame that is to be appended to 'df'
    tmp <- as.data.frame(classif_perf)
    tmp <- tmp[,c("FPR", "FNR")] # drop uninformative FNR and FNR 
    tmp$method <- rownames(tmp)
    rownames(tmp) <- NULL
    
    # convert to long format
    tmp <- reshape2::melt(data = tmp, 
                          id.vars = "method", 
                          variable.name = "measure", 
                          value.name = "value")
    
    # add a column for contamination level
    tmp$contam <- prop_crls
    
    # iteratively append to the main data frame
    df <- rbind(df, tmp)
    
  } # FOR prop_crlss
  
  # limit data frame to informative cases (NA's only occur when there are either no inattentive or no careful respondents)
  df <- df %>% 
    filter(!is.na(value))
  
  ## change the factor levels
  method_labels <- c(reliability = "Reliability",
                     antonym = "Antonym",
                     longstring = "LongString")
  measure_labels <- c(FPR = "False positive rate",
                      FNR = "False negative rate")
  df$method <- factor(method_labels[as.character(df$method)], levels = method_labels)
  df$measure <- factor(measure_labels[as.character(df$measure)], levels = measure_labels)
  
  
  ## make the plot
  df$`Method` <- df$method
  # p <- 
  #   ggplot(data = df, mapping = aes(x = contam, y = value, color = `Method`)) +
  #   geom_line() +
  #   ylim(0,1) +
  #   facet_grid(cols = vars(measure), scales = "free_x") +
  #   #scale_color_manual("", values = c("black", "#F8766D", "#619CFF")) +
  #   xlab("Prevalence of partially careless respondents") + 
  #   theme_bw() +
  #   theme(legend.position = "right", axis.title.y = element_blank()) +
  #   guides(color = guide_legend(ncol = 1))
  # 
  # # save the plot
  # ggsave(filename = paste0("plot-simresults_", suffix_type, "_", suffix_location, "_benchmarks", "_n=", n, "_classification_R=", R, ".pdf"),
  #        plot = p,
  #        path = dir_save,
  #        device = "pdf",
  #        width = 8,
  #        height = 2)
  
  ## append df
  df$setting <- suffix_location
  
  df_all <- rbind(df, df_all)
  
} # FOR suffix_location


## prepare joint plot
df_all$Setting <- factor(
  recode(df_all$setting,
         mixed = "Baseline",
         early = "Early onset",
         late = "Late onset"),
  levels = c("Baseline", "Late onset", "Early onset"))


p_all <- 
  ggplot(data = df_all, mapping = aes(x = contam, y = value, color = `Method`)) +
  geom_line() +
  theme_bw() +
  ylim(0,1) +
  ggh4x::facet_grid2(rows = vars(measure),  # separate x-axes
                     scales = "free_x", independent = "x", switch = "y") +
  labs(y = NULL) + 
  # scale_color_brewer(palette = "Dark2") +
  scale_color_manual(values = c("#A3A500", "#00BF7D", "#E76BF3")) + 
  xlab("Prevalence of partially careless respondents") + 
  theme(legend.position = "right", legend.key.width = unit(19, "pt"), 
        # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1),
        strip.background.y = element_blank(), 
        strip.placement = "outside", 
        strip.text.y = element_text(size = 11), 
        strip.switch.pad.grid = unit(0, "pt"),
        axis.title.y = element_blank()) +
  guides(color = guide_legend(ncol = 1))

# save the plot
ggsave(filename = paste0("plot-simresults_", suffix_type, "_benchmarks", "_n=", n, "_R=", R, ".pdf"),
       plot = p_all,
       path = dir_save,
       device = "pdf",
       width = 7.5,
       height = 4.25)
