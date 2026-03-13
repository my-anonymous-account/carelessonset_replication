
data_contamination <- function(data_uncontam, num_likert,
                               careless_start,
                               careless_end = rep(ncol(data_uncontam), nrow(data_uncontam)),
                               type = "random"){
  
  data_contam <- data_uncontam
  
  for(i in 1:nrow(data_contam)){
    
    interval <- careless_start[i]:careless_end[i]
    num_careless <- length(interval)
    
    if(type == "random")
    {
      tmp <- sample.int(num_likert, size = num_careless, replace = TRUE)
    } else if(type == "pattern"){
      # randomly select a combination of response categories of any subset size
      # (except single response categories to exclude staightliners)
      seq_likert <- seq_len(num_likert)
      combinations <- combn(seq_likert, 2, simplify = FALSE)
      keep <- sample.int(length(combinations), size = 1L)
      selected <- combinations[[keep]]
      # permute the selected combination (if more than one response category)
      if (length(selected) > 1L) selected <- sample(selected, replace = FALSE)
      # repeat the selected pattern
      tmp <- rep_len(selected, length.out = num_careless)
    } else stop("type needs to be either 'random' or 'pattern'")
    
    data_contam[i, interval] <- tmp
    
  } # FOR
  
  data_contam
  
} # FUN



#' Sample indices of careless respondents
#'
#' @param indices Vector of indices to sample careless respondents from
#' @param num_careless An integer vector denoting
#' the total numbers of careless respondents in the sample of size \code{length(indices)} that
#' we are interested in.
#'
#' @return A list of indices for each level of carelessness, as governed
#' by \code{num_careless}.
#'
#' @note The indices are nested, meaning that indices for a lower level of
#' carelessness will also be indices in all higher levels of carelessness.
sample_careless_idx <- function(indices, num_careless) {
  
  # randomly sample indices
  idx_careless <- sample(x = indices,
                         size = max(num_careless),
                         replace = FALSE)
  
  # construct list of sampled indices for different proportions of carelessness
  out <- lapply(num_careless, function(num) {
    if (num == 0L) integer()
    else idx_careless[seq_len(num)]
  })
  
  # add names and return indices
  names(out) <- names(num_careless)
  out
}


#' Sample indices of careless respondents by type (5 types in total)
#'
#' @param n sample size
#' @param prop_careless_arr A numeric vector with elements in [0,1] denoting
#' the proportions of careless respondents in the sample of size \code{n} that
#' we are interested in.
#' @param exclude_type carelessness type to be excluded
#'
#' @return A list of indices for each level of carelessness by type of carelessness, as governed
#' by \code{prop_careless}.
#'
#' @note The indices are nested, meaning that indices for a lower level of
#' carelessness will also be indices in all higher levels of carelessness.
sample_careless_idx_types <- function(n, prop_crls)
{
  idx_set <- seq_len(n)
  
  # specify group size 
  num_careless_per_group <- floor(prop_crls * n / 2)
  names(num_careless_per_group) <- prop_crls
  last_prop <- length(prop_crls)
  
  # pattern
  ls_pattern <- sample_careless_idx(indices = idx_set, num_careless = num_careless_per_group)
  idx_set     <- setdiff(idx_set, ls_pattern[[last_prop]])
  
  # random
  ls_random <- sample_careless_idx(indices = idx_set, num_careless = num_careless_per_group)
  idx_set    <- setdiff(idx_set, ls_random[[last_prop]])
  
  out <-
    lapply(seq_along(prop_crls), function(i){
      
      idx_random   <- ls_random[[i]]
      idx_pattern  <- ls_pattern[[i]]
      idx_careless <- c(idx_pattern, idx_random)
      idx_careful  <- setdiff(seq_len(n), idx_careless)
      
      return(list(random = idx_random, 
                  pattern = idx_pattern, 
                  careless = idx_careless,
                  careful = idx_careful,
                  num_careless_per_group = num_careless_per_group[i]))
    })
  
  names(out) <- prop_crls
  return(out)
}




## sample carelessness onset
carelessness_start_fn <- function(n, type, num_items = 300L)
{
  if(type == "early")
  {
    carelessness_start <- sample(seq(from = floor(0.1 * num_items), 
                                     to   = floor(0.5 * num_items)), size = n, replace = TRUE)
  } else if(type == "late")
  {
    carelessness_start <- sample(seq(from = floor(0.5 * num_items), 
                                     to   = floor(0.9 * num_items)), size = n, replace = TRUE)
  } else if(type == "mixed")
  {
    carelessness_start <- sample(seq(from = floor(0.1 * num_items), 
                                     to   = floor(0.9 * num_items)), size = n, replace = TRUE)
  } else stop("type not in c(early,late,mixed)")
  
  return(carelessness_start)
}


## the second changepoint; for the additioal simulations
carelessness_end_fun <- function(n)
{
  sample(seq(from = 151L, to = 270L), size = n, replace = TRUE)
} # FUN
