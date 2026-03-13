rm(list = ls()) ; cat("\014")
library("ggplot2")


NA_type <- "all"
dir_save <- "application/plots"


# load results
load(paste0("application/results/", "johnson2005",
            NA_type, "_results.Rdata"))

alpha_arr <- as.character(c(0.001, 0.005, 0.01))

df_location <- df_mean <- NULL

for(alpha. in alpha_arr)
{
  onset <- cp_LS_RE[[alpha.]]
  flagged <- which(!is.na(onset))
  onset_location <- sort(onset[flagged], decreasing = FALSE)
  mn <- mean(onset_location) 
  df_tmp <- data.frame(Respondent = seq_along(onset_location), Onset = onset_location, level = alpha.)
  df_location <- rbind(df_location, df_tmp)
  df_mean <- rbind(df_mean, data.frame(mean = mn, level = alpha.))
}

df_location$level <- paste0('alpha*" = "*', as.character(df_location$level)) 
df_mean$level <- paste0('alpha*" = "*', as.character(df_mean$level)) 


lp <- ggplot(data = df_location, aes(y = Onset, x = Respondent)) +
  facet_grid(cols = vars(level), scales = "free_x", labeller = label_parsed) + 
  geom_line(color = "blue") +
  geom_hline(data = df_mean, aes(yintercept = mean), linetype = "dashed") +
  # theme(axis.text=element_text(size=24),
  #       axis.title=element_text(size=24),
  #       legend.text=element_text(size=20)) +
  xlab("Index of flagged respondent") +
  ylab("Carelessness onset item") +
  theme_bw()


ggsave(filename = "onset-location-line.pdf", plot = lp, device = "pdf", 
       path = dir_save, width = 7.5, height = 2.5)
