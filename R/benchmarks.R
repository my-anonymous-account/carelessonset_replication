# caculate benchmark methods here. NB: NULL is returned in case of error, account for this possibility in the simulations!


## calculate psychometric antonyms (so that we allow for negatively correlated item pairs)
# 'x' is data matrix, 'min_corr' is minimum correlation for item pairs
# (the default is -0.6, which is also the default value in package 'careless'),
# 'min_pairs' is the minimum number of item pairs that should survive the
# correlation screening (the default value 6 is motivated by Goldammer et al.
# (2020; Table 9), who recommend who recommend using psychometric synonyms only with
# more than 5 item pairs of sufficient correlation)
antonym <- function(x, min_corr = -0.6, min_pairs = 6L) {
  
  nam <- colnames(x)
  
  if(!is.null(nam))
  {
    if(any(duplicated(nam)))
    {
      # prevent bug in careless package by assigning non-unique names 
      colnames(x) <- paste0("v", seq_len(ncol(x)))
    }
  }
  
  # in case of an error, return NAs
  tryCatch({
    # check whether we have enough item pairs with minimum correlation
    cor_mat <- cor(x)
    cor_ordered <- sort(cor_mat[lower.tri(cor_mat)], decreasing = FALSE)
    if (sum(cor_ordered < min_corr) < min_pairs) {
      # if we don't have enough item pairs,
      min_corr <- mean(cor_ordered[min_pairs + 0:1])
    }
    # compute psychometric antonyms
    careless::psychsyn(x = x, critval = min_corr, anto = TRUE, # anto = TRUE ensures use of antonyms
                       diag = FALSE, resample_na = TRUE)
  },
  error = function(e) rep(NA_real_, nrow(x)))
}


## calculate personal reliability (also called even-odd consistency)
# 'x' is a data matrix where the items are ordered according to the scales,
# 'scale_sizes' is a vector of integers specifying the length of each scale
# in the dataset
reliability <- function(x, scale_sizes) {
  # in case of an error, return NAs
  tryCatch({
    # omit uninformative warning reminding us that high values are indicative of carelessness
    suppressWarnings(
      careless::evenodd(x = x, factors = scale_sizes, diag = FALSE)
    )
  }, error = function(e) rep(NA_real_, nrow(x)))
}


## longstring index of johnson
# x is a data matrix
longstring <- function(x) {
  # in case of an error, return NAs
  tryCatch({
    careless::longstring(x = x, avg = FALSE)
  }, error = function(e) rep(NA_real_, nrow(x)))
}


## recode items to be positively worded
# x are response df
# keys contain keys for positive and negative wording of an item. Values in {-1,1}
# xmax is numeric value of highest response option, xmin that of the lowest; usually 5 and 1 in 5-point likert scales
recode <- function(x, keys, xmax = max(x), xmin = min(x))
{
  
  ## input checks
  stopifnot(all.equal(ncol(x), length(keys)))
  stopifnot(all(keys %in% c("1", "-1")))
  
  ## only the negatively worded items need to be recoded
  idx_recoded <- which(keys == "-1")
  x_out <- x
  
  for(i in idx_recoded)
  {
    ## extract responses to i-th item and recode them
    responses <- x[,i]
    responses <- xmax + xmin - responses
    x_out[,i] <- responses
  } # FOR
  
  return(x_out)
  
}
