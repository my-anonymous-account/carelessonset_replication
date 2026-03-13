# load preprocessed data
load("application/data/johnson2005data_preprocessed.Rdata")
# dim(responses_given) # 20,993 x 300
num_items <- ncol(responses_given)


# specify layers for auteoncoders
num_scales <- 30 # 30 facets

# significance levels
alpha  <- c(0.01, 0.005, 0.001) 


# number of cores
mc_cores <- 6L

# how to handle missing values>
NA_handling <- c("all", "completecases")
NA_handling = "all"


for(NA_type in NA_handling)
{
  
  if(NA_type == "all")
  {
    # there are NAs in the data that we shouldn't drop for compatibility with the paper, so assign them value 0
    # the non-NA responses therefore start with 0
    #mean(is.na(responses_given)) # 0.005349561
    data <- responses_given
    data[is.na(data)] <- 0L
    num_likert <- 5L + 1L # the NA answer category increases number of answer categories
    
  } else {
    data <- na.omit(responses_given)
    num_likert <- 5L
  } # IF
  
  # sample size
  n <- nrow(data)
  

  
  ## run method on 2-dim series
  crssonset <- carelessonset::carelessonset(responses = data, 
                                            num_scales = num_scales,
                                            num_likert = num_likert,
                                            longstring = TRUE,
                                            alpha = alpha,
                                            mc_cores = mc_cores,
                                            epochs = 100, 
                                            maxlen = 5L, # for consistency in paper: don't consider NAs here
                                            batch_size = 10,
                                            seed_tf = 12345, 
                                            seed_R = 12345)
  cp_LS_RE <- crssonset$changepoints
  
  
  ## run univariate CP detection on the individual series
  cp_RE <- carelessonset::changepoints_univariate(data = crssonset$series$RE, 
                                                  alpha = alpha, 
                                                  theta = "mean", 
                                                  mc_cores = mc_cores, 
                                                  matrix = FALSE,
                                                  teststat = TRUE)$changepoints
  
  cp_LS <- carelessonset::changepoints_univariate(data = crssonset$series$LSP_noise, 
                                                  alpha = alpha, 
                                                  theta = "mean", 
                                                  mc_cores = mc_cores, 
                                                  matrix = FALSE,
                                                  teststat = TRUE)$changepoints
  

  # save results
  save(cp_LS_RE, cp_LS, cp_RE, n, num_items, data,
       file = paste0("application/results/johnson2005", 
                     NA_type, "_results.Rdata"))
  
  # save scores and lsp separately for size reasons
  scores = crssonset$series$RE
  save(scores,
       file = paste0("application/results/johnson2005", 
                     NA_type, "_results_scores.Rdata"))
                     
  lsp = crssonset$series$LSP
  save(lsp,
       file = paste0("application/results/johnson2005", 
                     NA_type, "_results_LSP.Rdata"))
  
} # FOR 
