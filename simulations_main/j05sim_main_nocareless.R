source("R/contamination_middling.R")
source("R/dgp_johnson2005_cleandata.R")
source("R/benchmarks.R")


# initialize
R              <- 100L
n              <- 500L
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


# number of items
num_items <- scale_size * num_scales

# number alphas
num_alphas <- length(alpha)

## randomly sample the carelessness cutoffs
# we deliberately sample too many (n instead of num_careless), this makes the later code easier
set.seed(3124851)

## randomly sample indices of careful respondents (stays fixed across repetitions)
careless_indices_ls <- sample_careless_idx_types(n = n, 
                                                 prop_crls = c(0.02, 0.04, seq(from = 0.1, to = 0.6, by = 0.1), 0.8, 1.0),
                                                 exclude_type = "middling") # dead code here
careless_start <- carelessness_start_fn(n = n, type = "early") # dead code here


# initialize output list
out_CP_LS_SC <- out_CP_LS <- out_CP_SC <- rep(list(NA), R)

# initialize output for benchmarks
longstring_mat <- antonym_mat <- reliability_mat <- matrix(NA_real_, n, R)

tstart <- Sys.time()

for(r in seq_len(R)){
  
  # sample uncontaminated data (ordered as given by johnson2005)
  data_ls <- generate_data_johnson05marginals(n)
  data <- data_ls$data
  
  ## run method on 2-dim series
  crssonset <- carelessonset::carelessonset(responses = data, 
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
  longstring_mat[,r] <- longstring(data)
  
  ## reliability requires all items to be positively worded...
  data_posword <- recode(data, keys = data_ls$keys, xmax = 5L, xmin = 1L)
  
  ## ... and to be ordered by facets
  data_posword <- data_posword[,data_ls$idx_ordered]
  
  ## now run reliability method
  reliability_mat[,r] <- reliability(x = data_posword, scale_sizes = scale_sizes)
  
  
  ## psychometric antonym: doesn't require recoding or reordering of items
  colnames(data) <- paste0("v", seq_len(num_items)) # to prevent bug in pkg 'careless'
  antonym_mat[,r] <- antonym(x = data, min_corr = -0.6, min_pairs = 6L)
  
  print(paste0(r, " at ", Sys.time(), " for prop_crls=", 0.0))
  
} # FOR

# save
save(out_CP_LS_SC, out_CP_LS, out_CP_SC, alpha, mc_cores, reliability_mat, antonym_mat, longstring_mat,
     file = paste0(dir_save, "/", suffix_type, "_prop-crls=", 0.0, "_n=", n, "_R=", R, ".Rdata"))
print(paste0("Finished prop_crls=", 0.0, " after ", Sys.time() - tstart))

