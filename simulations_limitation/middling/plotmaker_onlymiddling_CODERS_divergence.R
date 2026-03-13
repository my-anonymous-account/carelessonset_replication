library("dplyr")
library("ggplot2")
library("tidyr")

dir_save <- "simulations_limitation/middling/plots"
dir_load <- "simulations_limitation/middling/analyzed"
suffix_type <- "j05sim-limit-middlingonly"
suffix_location <- "mixed"
target <- "absdiff_mean" # don't change: code below is specialized for 'absdiff_mean'

# specify the sizes
size_arr <- c(0.001, 0.005, 0.01)

# for each size, specify the contamination levels that we have simulated so far
# contamination levels = proportion of respondents who become careless
# don't include the no-careless case here (covered by classification plot)
prop_crlss_arr <- c(0.02, 0.04, seq(from = 0.1, to = 0.6, by = 0.1), 0.8, 1.0)
R <- 100L
n <- 500L

DF <- NULL

for(i in seq_along(size_arr))
{
  size <- size_arr[i]
  df <- NULL
  
  for(prop_crls in prop_crlss_arr)
  {
    load(paste0(dir_load, "/divergence_", suffix_type, "_", suffix_location, 
                "_prop-crls=", prop_crls, "_alpha=", size, "_n=", n,
                "_R=", R, ".Rdata"))
    
    div_perf_df <- NULL
    methods <- names(div_perf)
    
    for(m in seq_along(div_perf))
    {
      tmp <- as.data.frame(div_perf[[m]])
      tmp$type <- rownames(tmp)
      rownames(tmp) <- NULL
      tmp$method <- methods[m]
      div_perf_df <- rbind(tmp, div_perf_df)
    } # FOR methods
    
    
    # prepare the data frame that is to be appended to 'df'
    div_perf_df <- div_perf_df[, c(target, "fnr", "type", "method")] 
    
    # convert to long format
    div_perf_df <- reshape2::melt(data = div_perf_df, 
                          id.vars = c("method", "type"), 
                          variable.name = "measure", 
                          value.name = "value")
    
    # add a column for contamination level
    div_perf_df$contam <- prop_crls
    
    # iteratively append to the main data frame
    df <- rbind(df, div_perf_df)
    
  } # prop_crlss
  
  
  ## change the factor levels
  method_labels <- c(CP_LS = "LSP only", 
                     CP_SC = "RE only",
                     CP_LS_SC = "Both\ndimensions")
  measure_labels <- c(fnr = "False~negative~rate",
                      absdiff_mean = "Mean~absolute~error")
  df$method <- factor(method_labels[as.character(df$method)], levels = method_labels)
  df$measure <- factor(measure_labels[as.character(df$measure)], levels = measure_labels)
  df$type2 <- factor(
    recode(df$type, 
           careless = "Careless"),
    levels = c("Careless"))
           
  
  ## make the plot for the careless respondents only
  df_crls <- df %>% filter(type == "careless")
  df_scales <- split(df_crls, df_crls$measure)
  scales <- lapply(seq_along(df_scales), function(i) {
      x <- df_scales[[i]]
      if(i == 1)
      {
        scale_y_continuous(limits = c(0,1)) # FNR
      } else{
        scale_y_continuous(limits = c(x$ymin, x$ymax))
      }
  })

  # collect results
  DF <- rbind(DF, data.frame(df_crls, size = size))

} # size


DF$Size <- paste0(as.character(100 * DF$size), "%") # to-be-parsed labels

P <-  ggplot(data = DF, mapping = aes(x = contam, y = value, color = method, linetype = Size, alpha = Size)) +
  geom_line() +
  facet_grid(rows = vars(measure), scales = "free", labeller = label_parsed,  switch = "y") +
  labs(y = NULL) + 
  ggh4x::facetted_pos_scales(y = scales) +
  scale_linetype_manual("Significance\nlevel", values = c("dashed", "solid", "twodash")) +
  scale_alpha_manual("Significance\nlevel", values = c(0.4, 1, 0.4)) + 
  scale_color_manual("\n\nCODERS", values = c("#619CFF", "#F8766D", "black")) + # ugly hack to avoid top cutoff of legend
  xlab("Prevalence of partially careless respondents") + 
  theme_bw() +
  theme(legend.position = "right", 
        legend.key.width = unit(19, "pt"), 
        # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1), 
        strip.background.y = element_blank(), 
        strip.placement = "outside", 
        strip.text.y = element_text(size = 11), 
        strip.switch.pad.grid = unit(0, "pt")) +
  guides(color = guide_legend(nrow = 3, reverse = TRUE, order = 1))

# save the plots
ggsave(filename = paste0("plot-careless_", suffix_type, "_n=", n, "_R=", R, ".pdf"),
       plot = P,
       path = dir_save,
       device = "pdf",
       width = 7.5,
       height = 4.25)
