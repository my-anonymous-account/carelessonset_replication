## calculate Cronbach's alpha on population level

#' Calculate Cronbach's alpha for all scales
#' @param Rho correlation matrix
#' @param item_probabilities matrix of distributions (rows correspond to items,
#' columns to response categories)
#' @param num_scales number of scales
cronbachs_alpha <- function(Rho, item_probabilities, scale_size, num_scales) {
  # initializations
  seq_scale_size <- seq_len(scale_size)
  seq_num_scales <- seq_len(num_scales)
  # loop over scales and compute Cronbach's alpha
  sapply(seq_num_scales, function(i) {
    indices <- seq_scale_size + (i-1) * (scale_size)
    .cronbachs_alpha(Rho[indices, indices], item_probabilities[indices, ])
  })
}

#' Calculate Cronbach's alpha for a single scale
#' @param Rho correlation matrix
#' @param item_probabilities matrix of distributions (rows correspond to items,
#' columns to response categories)
.cronbachs_alpha <- function(Rho, item_probabilities) {
  # initializations
  num_items <- nrow(item_probabilities)
  num_likert <- ncol(item_probabilities)
  # compute variance of each item
  sigma2_items <- apply(item_probabilities, 1, function(p, x) {
    mu <- sum(p * x)
    sum(p * (x - mu)^2)
  }, x = seq_len(num_likert))
  # compute covariance matrix
  sigma_items <- sqrt(sigma2_items)
  Sigma <- tcrossprod(sigma_items) * Rho
  # compute variance of total score
  ones <- rep.int(1, num_items)
  sigma2_total <- drop(crossprod(ones, Sigma) %*% ones)
  # compute Cronbach's alpha
  num_items / (num_items-1) * (1 - sum(sigma2_items) / sigma2_total)
}
