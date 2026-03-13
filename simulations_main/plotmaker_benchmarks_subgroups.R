rm(list = ls()) ; cat("\014")

library("dplyr")
library("ggplot2")
library("tidyr")

dir_load <- "simulations_main/analyzed"
dir_save <- "simulations_main/plots"
suffix_type <- "j05sim-main"
suffix_locations <- c("mixed", "late", "early")

## contamination levels = proportion of respondents who become careless
prop_crlss_arr <- c(0.02, 0.04, seq(from = 0.1, to = 0.6, by = 0.1), 0.8, 1.0)
R <- 100L
n <- 500L

## initialize
df <- NULL


for(suffix_location in suffix_locations)
{
  for(prop_crls in c(prop_crlss_arr))
  {
    # load classification results
    load(paste0(dir_load, "/classification_", suffix_type, "_", suffix_location, 
                "_prop-crls=", prop_crls, "_benchmarks", "_n=", n,
                "_R=", R, ".Rdata"))
    
    
    # prepare the data frame that is to be appended to 'df'
    tmp <- as.data.frame(fnr_perf)
    tmp$method <- rownames(tmp)
    rownames(tmp) <- NULL
    
    # convert to long format
    tmp <- reshape2::melt(data = tmp, 
                          id.vars = "method", 
                          variable.name = "type", 
                          value.name = "fnr")
    
    # add a column for contamination level and setting
    tmp$contam <- prop_crls
    tmp$setting <- suffix_location
    
    # iteratively append to the main data frame
    df <- rbind(df, tmp)
    
  } # FOR prop_crlss
} # FOR suffix_location



## change the factor levels
method_labels <- c(reliability = "Reliability",
                   antonym = "Antonym",
                   longstring = "LongString")
df$method <- factor(method_labels[as.character(df$method)], levels = method_labels)
df$Type <- factor(
  recode(df$type,
         careless = "Careless",
         random = "Random",
         extreme = "Extreme",
         straightlining = "Straightlining",
         pattern = "Pattern"),
  levels = c("Careless", "Random", "Extreme", "Straightlining", "Pattern"))

df$Setting <- factor(
  recode(df$setting,
         mixed = "Baseline",
         early = "Early",
         late = "Late"),
  levels = c("Early", "Baseline", "Late"))

## make the plot for the 4 subgroups
df_subgroups <- df %>% filter(Type != "Careless")
df_subgroups$Method <- df_subgroups$method

p_subgroups <- 
  ggplot(data = df_subgroups, mapping = aes(x = contam, y = fnr, color = Method)) +
  geom_line() +
  ylim(c(0,1)) +
  ylab("False negative rate") +
  facet_grid(rows = vars(Type), cols = vars(Setting), scales = "free") +
  # scale_color_brewer(palette = "Dark2") +
  scale_color_manual(values = c("#A3A500", "#00BF7D", "#E76BF3")) +
  xlab("Prevalence of partially careless respondents") + 
  theme_bw() +
  theme(legend.position = "right", legend.key.width = unit(19, "pt"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1),
        #axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))


# save the plot
ggsave(filename = paste0("plot-simresults_", suffix_type, "_benchmarks", "_n=", n, "_subgroups_R=", R, ".pdf"), 
       plot = p_subgroups,
       path = dir_save, 
       device = "pdf",
       width = 7.5,
       height = 4.5) 
