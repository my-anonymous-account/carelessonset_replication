# read DAT20993.doc for information on the dataset

## import raw dataset
# besides responses to 300 IPIP-NEO-300 items, the data contain info on SEX and AGE, as well as 30 facet scores plus 5 domain scores (computed from the responses) of the IPIP Big 5. Furthermore, the variables RNDSPLIT (Jackson split-half intra-individual reliability coefficient) and RNDCORR (Goldberg Psychometric Antonyms) measure carelessness
# there are 20,993 respondents in the dataset. The original sample size reported in the paper is 23,994, but the authors eliminated participants with more than 10 missing responses (p. 116), leaving 20,993 respondents.
data_raw <- read.table("application/data/ipip20993.dat", header = TRUE, sep = "\t")
data_raw_nam <- colnames(data_raw)

## preprocess the responses
# 0 means NA
# Note that there is a single cell that is NA already (row 25, col 3), which is NA also in the SPSS versions of the dataset. So this is likely no error in the .dat version
responses_raw <- data_raw[,data_raw_nam %in% paste0("I", 1:300)]
responses_raw[responses_raw == 0] <- NA

table(unlist(responses_raw), useNA = "always")

## information on scoring
# Responses to negative items have already been recoded to be positively keyed in the raw data, so we need to decode them again to obtain the given responses
keys_data <- readxl::read_xlsx("application/data/IPIP-NEO-300 scoring tool_2.xlsx",
                               sheet = "Input")
key <- keys_data[,"Sign"][[1]]
key_positive <- startsWith(key, "+") # is the item positively worded?

## recoding
# the responses observed in this dataset have been recoded as follows:
# observed_response = num_likert + 1 - given_response
# rearrange this to recover the originally given responses
n <- nrow(responses_raw)
p <- ncol(responses_raw)
num_likert <- 5L
responses_given <- matrix(NA, nrow = n, ncol = p)
colnames(responses_given) <- colnames(responses_raw)
for(i in seq_len(p))
{
  # raw responses to i-th item
  raw_i <- as.integer(responses_raw[,i])
  
  if(key_positive[i])
  {
    # originally positively worded item, so no recoding necessary
    responses_given[,i] <- raw_i
  } else{
    # originally negatively worded item, so recoding necessary
    responses_given[,i] <- num_likert + 1L - raw_i
  } # IF
} # FOR


## carelessness scores
reliability <- data_raw$RNDSPLIT
antonym <- data_raw$RNDCORR

### playground: how were NAs handled?
# I suspect they were just dropped -> indeed
keys_facet <- sub(".", "", key)
facets_membership <- keys_data[,"Facet"][[1]] # facet membership of each item
facets <- cbind(key = keys_facet, facet = facets_membership)
facets <- facets[!duplicated(facets),] 
facets <- facets[order(facets[,"key"]),] # info on the 30 facets

# consistency in facet naming
scores_idx <- 303:332 # indices in data_raw_nam that denote facet scores
facet_nam_raw <- data_raw_nam[scores_idx]
facet_nam_rec <- facet_nam_raw
facet_nam_rec[facet_nam_rec == "ANXIETY"] <- "Anxiety"
facet_nam_rec[facet_nam_rec == "FRIENDLI"] <- "Friendliness"
facet_nam_rec[facet_nam_rec == "IMAGINAT"] <- "Imagination"
facet_nam_rec[facet_nam_rec == "TRUST"] <- "Trust"
facet_nam_rec[facet_nam_rec == "SELFEFFI"] <- "Self-Efficacy"
facet_nam_rec[facet_nam_rec == "ANGER"] <- "Anger"
facet_nam_rec[facet_nam_rec == "GREGAR"] <- "Gregariousness"
facet_nam_rec[facet_nam_rec == "ARTISTIC"] <- "Artistic Interests"
facet_nam_rec[facet_nam_rec == "MORALITY"] <- "Morality"
facet_nam_rec[facet_nam_rec == "ORDER"] <- "Orderliness"
facet_nam_rec[facet_nam_rec == "DEPRESS"] <- "Depression"
facet_nam_rec[facet_nam_rec == "ASSERTIV"] <- "Assertiveness"
facet_nam_rec[facet_nam_rec == "EMOTION"] <- "Emotionality"
facet_nam_rec[facet_nam_rec == "ALTRUISM"] <- "Altruism"
facet_nam_rec[facet_nam_rec == "DUTIFUL"] <- "Dutifulness"
facet_nam_rec[facet_nam_rec == "SELFCONS"] <- "Self-Consciousness"
facet_nam_rec[facet_nam_rec == "ACTIVITY"] <- "Activity Level"
facet_nam_rec[facet_nam_rec == "ADVENTUR"] <- "Adventurousness"
facet_nam_rec[facet_nam_rec == "COOPERAT"] <- "Cooperation"
facet_nam_rec[facet_nam_rec == "ACHIEVE"] <- "Achievement-Striving"
facet_nam_rec[facet_nam_rec == "IMMODERA"] <- "Immoderation"
facet_nam_rec[facet_nam_rec == "EXCITE"] <- "Excitement-Seeking"
facet_nam_rec[facet_nam_rec == "INTELLEC"] <- "Intellect"
facet_nam_rec[facet_nam_rec == "MODESTY"] <- "Modesty"
facet_nam_rec[facet_nam_rec == "SELFDISC"] <- "Self-Discipline"
facet_nam_rec[facet_nam_rec == "VULNERAB"] <- "Vulnerability"
facet_nam_rec[facet_nam_rec == "CHEERFUL"] <- "Cheerfulness"
facet_nam_rec[facet_nam_rec == "LIBERAL"] <- "Liberalism"
facet_nam_rec[facet_nam_rec == "SYMPATHY"] <- "Sympathy"
facet_nam_rec[facet_nam_rec == "CAUTIOUS"] <- "Cautiousness"

## recover facet scores (use raw data, as they are positively coded)
replicate_scores <- 
  sapply(seq_along(facet_nam_raw), function(i) {
    facet_i <- facet_nam_rec[i]
    responses_facet <- responses_raw[, facets_membership == facet_i]
    
    scores_reproduced <- rowSums(responses_facet, na.rm = TRUE)
    scores_reported <- data_raw[, data_raw_nam == facet_nam_raw[i]]
    all(scores_reproduced == scores_reported)
  })

# successful replication! The only FALSE is due to the one missing response that was not recoded as 0 (respondent 25, item 3)
replicate_scores

## save data
save(responses_given, antonym, reliability, 
     file = "application/data/johnson2005data_preprocessed.Rdata")


## recoded responses
responses_raw <- as.matrix(responses_raw)
save(responses_raw, file = "application/data/johnson2005data_raw_preprocessed.Rdata")
