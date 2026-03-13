rm(list = ls()) ; cat("\014")
library("ggplot2")
library("dplyr")

dir_load <- "simulations_limitation/shortsurvey/analyzed"
dir_save <- "simulations_limitation/shortsurvey/plots"
suffix_type <- "j05sim-limit-shortsurvey"

n <- 500L
R <- 100L

load(paste0(dir_load, "/gatheredresults_CODERS_", 
       suffix_type, "_n=", n, "_R=", R, ".Rdata"))

###### subgroups: FNR and MAE ##############
df_subgroups <- RESULTS %>% 
  filter(type %in% c("Random", "Extreme", "Straightlining", "Pattern"))
df_fnr <- df_subgroups %>% filter(measure == "False~negative~rate")
df_mae <- df_subgroups %>% filter(measure == "Mean~absolute~error")

## FNR
P_crls_fnr <- 
  ggplot(data = df_fnr, 
         mapping = aes(x = prop_crls, y = value, color = method, linetype = sig_level, alpha = sig_level)) +
  geom_line() +
  ylim(0,1) +
  facet_grid(rows = vars(type), cols = vars(numitems), scales = "fixed", labeller = label_parsed) +
  scale_linetype_manual("Significance\nlevel", values = c("dashed", "solid", "twodash")) +
  scale_alpha_manual("Significance\nlevel", values = c(0.4, 1, 0.4)) + 
  scale_color_manual("CODERS", values = c("#619CFF", "#F8766D", "black")) +
  xlab("Prevalence of partially careless respondents") + 
  theme_bw() +
  theme(legend.position = "top", legend.key.width = unit(19, "pt"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1)) +
  ylab("False negative rate") +
  guides(color = guide_legend(nrow = 1, reverse = TRUE, order = 1))

ggsave(filename = paste0("plot-subgroups_", suffix_type, "_n=", n, "_R=", R, "_fnr.pdf"),
       plot = P_crls_fnr,
       path = dir_save,
       device = "pdf",
       width = 7.5,
       height = 6.5)


## MAE
P_crls_mae <- 
  ggplot(data = df_mae, 
         mapping = aes(x = prop_crls, y = value, color = method, linetype = sig_level, alpha = sig_level)) +
  geom_line() +
  facet_grid(rows = vars(type), cols = vars(numitems), scales = "fixed", labeller = label_parsed) +
  scale_linetype_manual("Significance\nlevel", values = c("dashed", "solid", "twodash")) +
  scale_alpha_manual("Significance\nlevel", values = c(0.4, 1, 0.4)) + 
  scale_color_manual("CODERS", values = c("#619CFF", "#F8766D", "black")) +
  xlab("Prevalence of partially careless respondents") + 
  theme_bw() +
  theme(legend.position = "top", legend.key.width = unit(19, "pt"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1)) +
  ylab("Mean absolute error") +
  guides(color = guide_legend(nrow = 1, reverse = TRUE, order = 1))

ggsave(filename = paste0("plot-subgroups_", suffix_type, "_n=", n, "_R=", R, "_mae.pdf"),
       plot = P_crls_mae,
       path = dir_save,
       device = "pdf",
       width = 7.5,
       height = 6.5)



####### careless respondents aggregated ########
df_careless <- RESULTS %>% 
  filter(type == "Careless")

# suppress the empty FPR facet by removing NAs (FPR in careless respondents which cannot be computed)
df_careless <- na.omit(df_careless)

# fix axes limits for the two facets
df_scales <- split(df_careless, df_careless$measure)
scales <- lapply(seq_along(df_scales), function(i) {
  x <- df_scales[[i]]
  if(i == 1)
  {
    scale_y_continuous(limits = c(0,1)) # FNR
  } else{
    scale_y_continuous(limits = c(x$ymin, x$ymax))
  }
})


P_crls <- 
  ggplot(data = df_careless, 
         mapping = aes(x = prop_crls, y = value, color = method, linetype = sig_level, alpha = sig_level)) +
  geom_line() +
  facet_grid(rows = vars(measure), cols = vars(numitems), scales = "free_y", labeller = label_parsed, switch = "y") +
  labs(y = NULL) + 
  scale_linetype_manual("Significance\nlevel", values = c("dashed", "solid", "twodash")) +
  scale_alpha_manual("Significance\nlevel", values = c(0.4, 1, 0.4)) + 
  scale_color_manual("CODERS", values = c("#619CFF", "#F8766D", "black")) +
  ggh4x::facetted_pos_scales(y = scales) +
  xlab("Prevalence of partially careless respondents") + 
  theme_bw() +
  theme(legend.position = "top", 
        legend.key.width = unit(19, "pt"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1), 
        strip.background.y = element_blank(), 
        strip.placement = "outside", 
        strip.text.y = element_text(size = 11), 
        strip.switch.pad.grid = unit(0, "pt")) +
  guides(color = guide_legend(nrow = 1, reverse = TRUE, order = 1))

ggsave(filename = paste0("plot-careless_", suffix_type, "_n=", n, "_R=", R, ".pdf"),
       plot = P_crls,
       path = dir_save,
       device = "pdf",
       width = 7.5,
       height = 4.33)

####### attentive respondents ########
df_attentive <- RESULTS %>% 
  filter(type == "Attentive")

# suppress empty facets (MAE and FNR, which cannot be computed for attentive respondents)
df_attentive <- na.omit(df_attentive)

P_att <- 
  ggplot(data = df_attentive, 
         mapping = aes(x = prop_crls, y = value, color = method, linetype = sig_level, alpha = sig_level)) +
  geom_line() +
  ylim(0,1) + 
  facet_grid(cols = vars(numitems), labeller = label_parsed) +
  scale_linetype_manual("Significance\nlevel", values = c("dashed", "solid", "twodash")) +
  scale_alpha_manual("Significance\nlevel", values = c(0.4, 1, 0.4)) + 
  scale_color_manual("CODERS", values = c("#619CFF", "#F8766D", "black")) + # ugly hack to avoid top cutoff of legend
  xlab("Prevalence of partially careless respondents") + 
  ylab("False positive rate") +
  theme_bw() +
  theme(legend.position = "top", 
        legend.key.width = unit(19, "pt"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1)) +
  guides(color = guide_legend(nrow = 1, reverse = TRUE, order = 1))

ggsave(filename = paste0("plot-attentive_", suffix_type, "_n=", n, "_R=", R, ".pdf"),
       plot = P_att,
       path = dir_save,
       device = "pdf",
       width = 7.5,
       height = 2.8)
