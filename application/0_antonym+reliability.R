rm(list = ls()) ; cat("\014")

# load antonyms + reliability coefficients
load("application/data/johnson2005data_preprocessed.Rdata")
load("application/data/johnson2005data_raw_preprocessed.Rdata")


#### antonym analysis (no reversely coded items here!)
antonym_new <- careless::psychant(responses_given)
idx_flagged_antonym_new <- which(antonym_new > -0.03)  # cutoff suggested by Johnson (p. 118); see also guidelines by Goldammer et al (2020; Tab 1)
length(idx_flagged_antonym_new) # 284


#### personal reliability analysis
## get facet membership of each item
# Responses to negative items have already been recoded to be positively keyed in the raw data, so we need to decode them again to obtain the given responses
keys_data <- readxl::read_xlsx("application/data/IPIP-NEO-300 scoring tool_2.xlsx",
                               sheet = "Input")
key_membership <- keys_data[,"Key"][[1]]
key_sign <- substr(keys_data[,"Sign"][[1]], 1, 1)
facets_membership <- keys_data[,"Facet"][[1]] # facet membership of each item
factor_membership <- substr(keys_data[,"Key"][[1]], 1, 1)

# 10 items per facet; 30 items in total
table(key_membership)
factors <- rep(10, 30)
key_order <- order(key_membership)


# calculate personal reliability
#  Package warning: "Computation of even-odd has changed for consistency of interpretation with other indices. This change occurred in version 1.2.0. A higher score now indicates a greater likelihood of careless responding"
reliability_new <- careless::evenodd(x = responses_raw[,key_order], factors = factors) 
idx_flagged_reliability_new <- which(reliability_new > 0.3) # cutoff suggested by Johnson;  see also guidelines by Goldammer et al (2020; Tab 1)
length(idx_flagged_reliability_new) # 140 


save(reliability_new, antonym_new, idx_flagged_antonym_new, idx_flagged_reliability_new, key_order, responses_raw, responses_given, key_membership, facets_membership, key_sign, file = paste0("application/results/johnson2005reliability_antonyms_flagged.Rdata"))
