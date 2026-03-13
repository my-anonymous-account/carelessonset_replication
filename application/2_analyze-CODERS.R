rm(list = ls()) ; cat("\014")

NA_type <- "all"

# load visualizing functions
source("application/plotfuns.R")

# load antonyms + reliability coefficients as provided by Johnson
load("application/data/johnson2005data_preprocessed.Rdata")


# load results
load(paste0("application/results/johnson2005", 
            NA_type, "_results.Rdata"))

# load LSP and RE scores
load(paste0("application/results/johnson2005", 
            NA_type, "_results_LSP.Rdata"))
load(paste0("application/results/johnson2005", 
            NA_type, "_results_scores.Rdata"))


# how many samples were flagged?
methods <- c("cp_LS_RE", "cp_LS", "cp_RE")
flagged <- lapply(seq_along(methods), function(k){
  x <- get(methods[k])
  num_flagged <- lapply(seq_along(x), function(i){
    crlss <- !is.na(x[[i]])
    crlss_idx <- which(crlss)
    sm <- sum(crlss)
    oset <- x[[i]][crlss] 
    scores_att <- mean(scores[!crlss,])
    SC_crlss_mean_pre <- SC_crlss_mean_post <- SC_crlss_sd_pre <- SC_crlss_sd_post <- 
      LS_crlss_mean_pre <- LS_crlss_mean_post <- LS_crlss_sd_pre <- LS_crlss_sd_post <- 
      rep(NA_real_, sm)
    for(j in seq_along(crlss_idx))
    {
      scores_pre <- scores[crlss_idx[j], seq_len(oset[j]-1L)]
      scores_post <- scores[crlss_idx[j], oset[j]:num_items]
      lsp_pre <- lsp[crlss_idx[j], seq_len(oset[j]-1L)]
      lsp_post <- lsp[crlss_idx[j], oset[j]:num_items]
      SC_crlss_mean_pre[j] <- mean(scores_pre)
      SC_crlss_mean_post[j] <- mean(scores_post)
      SC_crlss_sd_pre[j] <- sqrt(var(scores_pre))
      SC_crlss_sd_post[j] <- sqrt(var(scores_post))
      LS_crlss_mean_pre[j] <- mean(lsp_pre)
      LS_crlss_mean_post[j] <- mean(lsp_post)
      LS_crlss_sd_pre[j] <- sqrt(var(lsp_pre))
      LS_crlss_sd_post[j] <- sqrt(var(lsp_post))
    }
    
    SC_att_mean <- apply(scores[!crlss,], 1, mean)
    SC_att_sd <- apply(scores[!crlss,], 1, function(u) sqrt(var(u)))
    LS_att_mean <- apply(lsp[!crlss,], 1, mean)
    LS_att_sd <- apply(lsp[!crlss,], 1, function(u) sqrt(var(u)))
    
    summary_stats <- 
      c(num = sm, prop = sm / n, onset = mean(oset), 
        SC_att_mean =  mean(SC_att_mean), 
        SC_att_sd = mean(SC_att_sd),
        SC_crlss_mean_pre = mean(SC_crlss_mean_pre),
        SC_crlss_sd_pre = mean(SC_crlss_sd_pre),
        SC_crlss_mean_post = mean(SC_crlss_mean_post),
        SC_crlss_sd_post = mean(SC_crlss_sd_post),
        LS_att_mean =  mean(LS_att_mean), 
        LS_att_sd = mean(LS_att_sd),
        LS_crlss_mean_pre = mean(LS_crlss_mean_pre),
        LS_crlss_sd_pre = mean(LS_crlss_sd_pre),
        LS_crlss_mean_post = mean(LS_crlss_mean_post),
        LS_crlss_sd_post = mean(LS_crlss_sd_post))
    
    individual_SC <- list(
      att = cbind(mean = SC_att_mean, sd = SC_att_sd),
      crlss_pre = cbind(mean = SC_crlss_mean_pre, sd = SC_crlss_sd_pre),
      crlss_post = cbind(mean = SC_crlss_mean_post, sd = SC_crlss_sd_post))
    
    individual_LS <- list(
      att = cbind(mean = LS_att_mean, sd = LS_att_sd),
      crlss_pre = cbind(mean = LS_crlss_mean_pre, sd = LS_crlss_sd_pre),
      crlss_post = cbind(mean = LS_crlss_mean_post, sd = LS_crlss_sd_post))
    
    individual <- list(idx_careless = crlss_idx, idx_att = which(!crlss),
                       RE = individual_SC,
                       LSP = individual_LS)
    list(summary = summary_stats, individual = individual)
    
  })
  names(num_flagged) <- names(x)
  num_flagged
})
names(flagged) <- methods



#### benchmark methods ----
# NOTE: the code below requires output from the script `application/0_antonym+reliability.R`, so run this first!
load("application/results/johnson2005reliability_antonyms_flagged.Rdata")


rel_stats <- c(num = length(idx_flagged_reliability_new),
               prop = length(idx_flagged_reliability_new) / n )
flagged$reliability <- list(summary = rel_stats,
                            idx_careless = idx_flagged_reliability_new)
ant_stats <- c(num = length(idx_flagged_antonym_new),
               prop = length(idx_flagged_antonym_new) / n )
flagged$antonym <- list(summary = ant_stats,
                        idx_careless = idx_flagged_antonym_new)

results <- flagged
save(results, file = "application/results.Rdata")



#### correlations
## NOTE: carelessness is indicated by low values of antonym and reliability, but large values of reconstruction errors. So multiply RMSRE by -1
cor(reliability_new, antonym_new)
rmsre <- -rowMeans(data)
cor(rmsre, reliability_new) # weak correlation (-0.1190495)
cor(rmsre, antonym_new) # -0.08675789
