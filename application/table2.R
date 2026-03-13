rm(list = ls()) ; cat("\014")
NA_type <- "all"

# load onset results
load(paste0("application/results/", "johnson2005",
            NA_type, "_results.Rdata"))

# n_baseline
N <- 22448 

## load self-computed antonyms + reliability coefficients
load("application/results/johnson2005reliability_antonyms_flagged.Rdata")

## update naming 
idx_flagged_antonym <- idx_flagged_antonym_new
idx_flagged_reliability <- idx_flagged_reliability_new
idx_flagged_reliability_antonym <- intersect(idx_flagged_antonym, idx_flagged_reliability) 

### table 2 ----
alphas <- rev(c("0.01", "0.005", "0.001"))
tab2 <- matrix(NA_real_, 3, 5)
colnames(tab2) <- c("num", "rel", alphas)
rownames(tab2) <- c("reliability", "antonym", "both")

tab2[,"num"] <- c(length(idx_flagged_reliability),
                  length(idx_flagged_antonym),
                  length(idx_flagged_reliability_antonym))
tab2[,"rel"] <- round(100 * tab2[,"num"] / N, 1)


for(a in seq_along(alphas))
{
  alpha <- alphas[a]
  onset <- cp_LS_RE[[alpha]]
  flagged <- which(!is.na(onset))
  
  intersect_reliability <- length(intersect(idx_flagged_reliability, flagged))
  intersect_antonym     <- length(intersect(idx_flagged_antonym, flagged))
  intersect_both        <- length(intersect(idx_flagged_reliability_antonym, flagged))

  tab2["reliability", alpha] <- intersect_reliability
  tab2["antonym", alpha]     <- intersect_antonym
  tab2["both", alpha]        <- intersect_both
}

tab2
# > tab2
#             num rel 0.001 0.005 0.01
# reliability 140 0.6     2     8   10
# antonym     284 1.3     9    18   24
# both         11 0.0     1     1    1
