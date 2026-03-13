rm(list = ls()) ; cat("\014")

NA_type <- "all"

# load results
load(paste0("application/results/", "johnson2005",
            NA_type, "_results.Rdata"))

N <- 22448 # n_baseline

### table 1 ----
alphas <- rev(c("0.01", "0.005", "0.001"))
tab1 <- matrix(NA_real_, 3, 5)
colnames(tab1) <- c("abs", "rel", "min", "mean", "max")
rownames(tab1) <- alphas

## overview of number of flagged respondents
for(a in seq_along(alphas))
{
  alpha <- alphas[a]
  onset <- cp_LS_RE[[alpha]]
  flagged <- which(!is.na(onset))
  locations <- onset[flagged]
  num_flagged <- length(flagged)
  tab1[a,] <- c(num_flagged, num_flagged / N, min(locations), mean(locations), max(locations))
}

tab1
# > tab1
#       abs         rel min     mean max
# 0.001 184 0.008196721  19 143.6902 292
# 0.005 571 0.025436565  19 145.8074 294
# 0.01  856 0.038132573  15 147.8703 294