rm(list = ls()) ; cat("\014")
NA_type <- "all"
alpha <- "0.001"

# load onset results
load(paste0("application/results/", "johnson2005",
            NA_type, "_results.Rdata"))

# load LSP and RE scores
load(paste0("application/results/johnson2005", 
            NA_type, "_results_LSP.Rdata"))
load(paste0("application/results/johnson2005", 
            NA_type, "_results_scores.Rdata"))

## indices of selected participants
idx <- c(12995, 14956, 6715, 19660, 17193, 11465) 
num_items <- ncol(scores)

## get onset locations
locations <- cp_LS_RE[[alpha]][idx] 
names(locations) <- idx
locations
# > locations
# 12995 14956  6715 19660 17193 11465 
#   148   246    63    58   142    NA 


## prepare list
out <- list()

for(i in seq_along(idx)) 
{
  
  # break loop for the respondent without flagged changepoint
  if(is.na(locations[i])) next
  
  id <- idx[i]
  cp <- locations[i]
  range_pre  <- seq_len(cp - 1)
  range_post <- seq(cp, num_items, 1L)
  
  ## prepare output
  outtab <- matrix(NA_real_, 2, 4)
  rownames(outtab) <- c("RE", "LSP")
  colnames(outtab) <- c("pre-mean", "pre-sd", "post-mean", "post-sd")
  
  ## reconstruction error stats
  outtab["RE", "pre-mean"]  <- mean(scores[id, range_pre])
  outtab["RE", "pre-sd"]    <- sd(scores[id, range_pre])
  outtab["RE", "post-mean"] <- mean(scores[id, range_post])
  outtab["RE", "post-sd"]   <- sd(scores[id, range_post])
  
  ## LSP stats
  outtab["LSP", "pre-mean"]  <- mean(lsp[id, range_pre])
  outtab["LSP", "pre-sd"]    <- sd(lsp[id, range_pre])
  outtab["LSP", "post-mean"] <- mean(lsp[id, range_post])
  outtab["LSP", "post-sd"]   <- sd(lsp[id, range_post])
  
  out[[i]] <- outtab
  
}


## stats for  respondent without flagged changepoint
i <- which(is.na(locations))
id <- idx[i]
outtab <- matrix(NA_real_, 2, 2)
rownames(outtab) <- c("RE", "LSP")
colnames(outtab) <- c("mean", "sd")
outtab["RE", "mean"]  <- mean(scores[id,])
outtab["RE", "sd"]    <- sd(scores[id,])
outtab["LSP", "mean"]  <- mean(lsp[id,])
outtab["LSP", "sd"]    <- sd(lsp[id,])
out[[i]] <- outtab

names(out) <- idx


lapply(seq_along(out), function(i) round(out[[i]], 3))
# [[1]]
#     pre-mean pre-sd post-mean post-sd
# RE     0.019  0.035     0.069   0.082
# LSP    4.122  1.428     3.373   1.464
# 
# [[2]]
#     pre-mean pre-sd post-mean post-sd
# RE     0.017  0.026     0.070   0.096
# LSP    2.886  1.189     3.236   1.815
# 
# [[3]]
#     pre-mean pre-sd post-mean post-sd
# RE     0.095  0.119     0.024   0.036
# LSP    2.565  0.917     2.668   0.925
# 
# [[4]]
#     pre-mean pre-sd post-mean post-sd
# RE     0.078  0.103     0.025   0.035
# LSP    2.158  0.621    19.708  11.027
# 
# [[5]]
#     pre-mean pre-sd post-mean post-sd
# RE     0.018  0.029     0.024   0.046
# LSP  146.000  0.000    23.597  25.878
# 
# [[6]]
#      mean    sd
# RE  0.079 0.105
# LSP 2.813 1.115


# load visualizing functions
source("application/plotfuns.R")

# create plots
dir_save <- "application/plots"
for(i in seq_along(idx))
{
  id <- idx[i]
  onset <- if (is.na(locations[i])) NULL else locations[i]
  pp <- plot_dimension(RE = scores[id,], LSP = lsp[id,], 
                       onset = onset, increase = FALSE)
  nam <- paste0("johnson2005_", letters[i], "_id", id, ".pdf")
  # as tall as possible while all six plots still fitting on one page
  ggsave(filename = nam, plot = pp, device = "pdf", 
         path = dir_save, width = 0.49 * 7.5, height = 0.49 * 4.785)
}


## how many respondents have LSP values larger than that of if = 19660?
maxlsp <- sapply(seq_len(nrow(lsp)), function(i) max(lsp[i,]))
sum(maxlsp > max(lsp[19660,])) # 4
ecdf_LSP <- ecdf(x = maxlsp)
ecdf_LSP(max(lsp[19660,]))
# 0.9998095


### quantile of average RE of selected careless respondent throughout
# get mean (squared & scaled) reconstruction error for each individual
RE <- rowMeans(scores)

# index of one (or more) observations whose empirical quantile we are interested in
# here we choose an observation that was not flagged by us, but by Johnson (2005)
id <- idx[which(is.na(locations))]  # old selection: 8474

# value 
(val_id <- RE[id])
# 0.07938923

# order statistics
RE_ord <- order(RE, decreasing = TRUE)
which(RE_ord == id) 
# 7

# empirical distribution
ecdf. <- ecdf(x = RE)

# percentile of the observation's RE
(prctl <- ecdf.(val_id))
# 0.9997142

# how many observations (in %) have bigger RE?
1.0 - prctl
# 0.0002858096

library("ggplot2")
p <- ggplot(data.frame(RE = RE), aes(x = RE)) + 
  geom_histogram(aes(y = ..density..), colour = "black", fill = "gray", bins = 80)+
  geom_density(alpha = 0.2, fill = "orange") +
  #geom_vline(aes(xintercept = meanRE),
  #           color = "red", linetype = "dashed", size = 0.8) +
  geom_vline(aes(xintercept = val_id),
             color = "blue", linetype = "dotted", size = 0.8) +
  #geom_text(aes(x = meanRE, label= paste0("Mean = ", round(meanRE, 4)), y = 75), colour = "red", angle = 0, hjust = -0.1) +
  geom_text(aes(x = val_id, label= paste0(round(val_id, 3)), y = 52.5), colour = "blue", angle = 0, hjust = -0.15) +
  xlab("RE averaged across items") +
  ylab("Frequency") +
  theme_bw()

ggsave(filename = paste0("RE_ECDF_id", id, ".png"), plot = p, device = "png", path = "application/plots", dpi = 500, width = 7.5, height = 3.33)
