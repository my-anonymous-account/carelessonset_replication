rm(list = ls()) ; cat("\014")
NA_type <- "all"


# load antonym+reliability indices
load("application/results/johnson2005reliability_antonyms_flagged.Rdata")

# load onset results
load(paste0("application/results/", "johnson2005",
            NA_type, "_results.Rdata"))

alphas <- rev(c("0.01", "0.005", "0.001"))

flagged_either <- lapply(seq_along(alphas), function(a){
  alpha <- alphas[a]
  onset <- cp_LS_RE[[alpha]]
  flagged <- which(!is.na(onset))
  return(c(union(flagged, union(idx_flagged_reliability_new, idx_flagged_antonym_new))))
})
names(flagged_either) <- alphas


# be very liberal in kicking out observations (so eps=0.01); round to 2nd decimal
# the raw data have already been recoded so that all items are positively worded
flagged <- flagged_either$`0.01`
cormat_clean_posword <- round(cor(responses_raw[-flagged , ], use = "pairwise.complete.obs"), 2)

# use the responses as given: that is, some are negatively worded
cormat_clean_negword <- round(cor(data[-flagged , ], use = "pairwise.complete.obs"), 2) 

matrixcalc::is.positive.definite(cormat_clean_posword) # TRUE
matrixcalc::is.positive.definite(cormat_clean_negword) # TRUE


cormat_clean_posword <- cormat_clean_posword[key_order, key_order]
cormat_clean_negword <- cormat_clean_negword[key_order, key_order]

p <- ncol(cormat_clean_posword)
df_posword <- matrix(NA_real_, p * p, 3L)
colnames(df_posword) <- c("x", "y", "Correlation")
df_negword <- df_posword

ct <- 1L
for(i in 1:p)
{
  for(j in 1:p)
  {
    df_posword[ct,] <- c(i,j,cormat_clean_posword[i,j])
    df_negword[ct,] <- c(i,j,cormat_clean_negword[i,j])
    ct <- ct + 1L
  }
}

df_negword <- data.frame(df_negword)
df_posword <- data.frame(df_posword)


library("ggplot2")
p_posword <- ggplot(df_posword, aes(x = x, y = y, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(high = "#075AFF",
                       mid = "white",
                       low = "#FF0000") +
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "right")

p_negword <- ggplot(df_negword, aes(x = x, y = y, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(high = "#075AFF",
                       mid = "white",
                       low = "#FF0000") +
  coord_fixed() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "right")

# ggsave("johnson2005_cleancormat_posword.pdf", plot = p_posword, device = "pdf", path = "application/plots", width = 0.75 * 7.5, height = 0.75 * 6)
ggsave("johnson2005_cleancormat_negword.pdf", plot = p_negword, device = "pdf", path = "application/plots", width = 0.75 * 7.5, height = 0.75 * 6)


## keep the clean data (as given, so no recoding of neg worded items)
data_given_clean <- data[-flagged,]

## same for recoded ones 
data_posword_clean <- responses_raw[-flagged, key_order]

## marginal response distributions of recoded data (also align order with order in cormat)
p <- ncol(data_posword_clean)
n <- nrow(data_posword_clean)

marginals <- t(sapply(seq_len(p), function(j){
  
  ## the positively worded responses, including NAs
  x_withNA <- data_posword_clean[,j]
  
  ## kick out NA responses
  x <- x_withNA[!is.na(x_withNA)]
  n_red <- length(x)
  
  ## calculate response probabilities
  dist <- as.numeric(table(factor(x, levels = 1:5)) / n_red)
  names(dist) <- 1:5
  return(dist)
}))
rownames(marginals) <- colnames(data_posword_clean)

# make plot of marginals
df_marginals <- reshape2::melt(marginals, id.vars = 1:5)
colnames(df_marginals) <- c("Item", "Response option", "Probability")
df_marginals$Item <- factor(df_marginals$Item, levels = rownames(marginals))
df_marginals$Item <- as.integer(df_marginals$Item)
df_marginals$`Response option` <- factor(df_marginals$`Response option`, levels = 1:5)


p_marginals <- ggplot(df_marginals, aes(x = Item, y = `Response option`, 
                                        fill = Probability)) +
  geom_tile() +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "white",
                       high = "black", 
                       limits = c(0,1), 
                       breaks = c(0.2, 0.4, 0.6, 0.8)) +
  theme(legend.position = "right") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0))

ggsave("johnson2005_marginals.pdf", plot = p_marginals, device = "pdf", path = "application/plots", width = 7.5, height = 2.67)


## save marginals data
save(marginals, data_given_clean, cormat_clean_posword, cormat_clean_negword, key_order, key_membership, facets_membership, key_sign, 
     file = "application/results/johnson2005_cleandata.Rdata")


## save covmat data
save(cormat_clean_posword, cormat_clean_negword, key_order, key_membership, facets_membership, key_sign, file = "application/results/johnson2005_cleancormat.Rdata")
