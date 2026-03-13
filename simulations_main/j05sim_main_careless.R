source("R/contamination_middling.R")
source("R/dgp_johnson2005_cleandata.R")
source("R/benchmarks.R")


# initialize
R              <- 100L
n              <- 500L ## can be changed to 248 and 1000
scale_size     <- 10L
num_scales     <- 30L
num_likert     <- 5L
batch_size     <- 10L
scale_sizes    <- rep(10L, 30L) # for the option 'factors' in reliability benchmark method
alpha          <- c(0.01, 0.005, 0.001) 
mc_cores       <- 40L 



# directory to save results in
dir_save <- "simulations_main/results"
suffix_type <- "j05sim-main"
suffix_location <- "late" ## change to `early` or `mixed`


# specify how many samples become careless
prop_crls_arr <- c(0.02, 0.04, seq(from = 0.1, to = 0.6, by = 0.1), 0.8, 1.0)


# number of items
num_items <- scale_size * num_scales

# number alphas
num_alphas <- length(alpha)

## randomly sample the carelessness cutoffs
# we deliberately sample too many (n instead of num_careless), this makes the later code easier
set.seed(3124851)

## randomly sample indices of careful respondents (stays fixed across repetitions)
careless_indices_ls <- sample_careless_idx_types(n = n, prop_crls = prop_crls_arr, exclude_type = "middling")
careless_start <- carelessness_start_fn(n = n, type = suffix_location)


for(a in seq_along(prop_crls_arr))
{
  
  prop_crls <- prop_crls_arr[a]
  careless_indices_ls_a <- careless_indices_ls[[a]]
  idx_pattern <- careless_indices_ls_a$pattern
  idx_random <- careless_indices_ls_a$random
  idx_straight <- careless_indices_ls_a$straight
  idx_extreme <- careless_indices_ls_a$extreme
  idx_careless <- careless_indices_ls_a$careless
  idx_careful <- careless_indices_ls_a$careful
  group_size_crls <- careless_indices_ls_a$num_careless_per_group
  
  # initialize output list for CODERS
  out_CP_LS_SC <- out_CP_LS <- out_CP_SC <- rep(list(NA), R)

  # initialize output for benchmarks
  longstring_mat <- antonym_mat <- reliability_mat <- matrix(NA_real_, n, R)
  
  tstart <- Sys.time()
  
  for(r in seq_len(R)){
    
    # sample uncontaminated data (ordered as given by johnson2005)
    data_ls <- generate_data_johnson05marginals(n)
    data <- data_ls$data
  
    # contaminate data
    # NOTE: we make the implicit assumption here that careless respondents will
    # always choose an answer category (so no response isn't an option). In simulated data, there are no NAs either
    data_random_tmp <- data_contamination(data_uncontam  = data[idx_random,,drop = FALSE],
                                          num_likert     = num_likert, 
                                          careless_start = careless_start[idx_random], 
                                          careless_end   = rep(num_items, group_size_crls), 
                                          type           = "random")
    data_straight_tmp <- data_contamination(data_uncontam  = data[idx_straight,,drop = FALSE],
                                            num_likert     = num_likert, 
                                            careless_start = careless_start[idx_straight], 
                                            careless_end   = rep(num_items, group_size_crls), 
                                            type           = "straightlining")
    data_extreme_tmp <- data_contamination(data_uncontam  = data[idx_extreme,,drop = FALSE],
                                           num_likert     = num_likert, 
                                           careless_start = careless_start[idx_extreme], 
                                           careless_end   = rep(num_items, group_size_crls), 
                                           type           = "extreme")
    data_pattern_tmp <- data_contamination(data_uncontam  = data[idx_pattern,,drop = FALSE],
                                           num_likert     = num_likert, 
                                           careless_start = careless_start[idx_pattern], 
                                           careless_end   = rep(num_items, group_size_crls), 
                                           type           = "pattern")
    
    
    # merge with uncontaminated data
    data_contam                <- data
    data_contam[idx_pattern,]  <- data_pattern_tmp
    data_contam[idx_random,]   <- data_random_tmp
    data_contam[idx_straight,] <- data_straight_tmp
    data_contam[idx_extreme,]  <- data_extreme_tmp
    
    ## run method on 2-dim series
    crssonset <- carelessonset::carelessonset(responses = data_contam, 
                                              num_scales = num_scales,
                                              num_likert = num_likert,
                                              longstring = TRUE,
                                              alpha = alpha,
                                              mc_cores = mc_cores,
                                              epochs = 100, 
                                              batch_size = batch_size,
                                              seed_tf = 12345, 
                                              seed_R = NULL)
    
    
    ## run univariate CP detection on the individual series
    cp_RE <- carelessonset::changepoints_univariate(
      data = crssonset$series$RE, 
      alpha = alpha, 
      theta = "mean", 
      mc_cores = mc_cores, 
      matrix = FALSE,
      teststat = TRUE, 
      CP_at_segment_end = FALSE)$changepoints
    
    cp_LS <- carelessonset::changepoints_univariate(
      data = crssonset$series$LSP_noise, 
      alpha = alpha, 
      theta = "mean", 
      mc_cores = mc_cores, 
      matrix = FALSE,
      teststat = TRUE, 
      CP_at_segment_end = FALSE)$changepoints


    ## organize results
    out_CP_LS_SC[[r]]    <- crssonset$changepoints
    out_CP_LS[[r]]       <- cp_LS
    out_CP_SC[[r]]       <- cp_RE
    
    
    ### run benchmark methods
    ## longstring
    longstring_mat[,r] <- longstring(data_contam)
    
    ## reliability requires all items to be positively worded...
    data_posword <- recode(data_contam, keys = data_ls$keys, xmax = 5L, xmin = 1L)
    
    ## ... and to be ordered by facets
    data_posword <- data_posword[,data_ls$idx_ordered]
    
    ## now run reliability method
    reliability_mat[,r] <- reliability(x = data_posword, scale_sizes = scale_sizes)
    
    
    ## psychometric antonym: doesn't require recoding or reordering of items
    colnames(data_contam) <- paste0("v", seq_len(num_items)) # to prevent bug in pkg 'careless'
    antonym_mat[,r] <- antonym(x = data_contam, min_corr = -0.6, min_pairs = 6L)
    
    print(paste0(r, " at ", Sys.time(), " for prop_crls=", prop_crls))
    
  } # FOR
  
  # save
  save(out_CP_LS_SC, out_CP_LS, out_CP_SC,
       careless_start, idx_careful, idx_careless, idx_pattern, idx_random, idx_straight, idx_extreme,
       prop_crls, alpha, n, num_items, mc_cores, reliability_mat, antonym_mat, longstring_mat,
       file = paste0(dir_save, "/", suffix_type, "_", suffix_location, "_prop-crls=", prop_crls, "_n=", n,  "_R=", R, ".Rdata"))
  print(paste0("Finished prop_crls=", prop_crls, " after ", Sys.time() - tstart))
  
} # FOR prop_crlss

