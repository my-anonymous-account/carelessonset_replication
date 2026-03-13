# the R implementation below will be slow for many changepoints; consider rewriting in C++!

# length_series is scalar, changepoints a vector. There will be (num_changepoints + 1) segments
segmentation <- function(length_series = integer(), changepoints = integer())
{
    out <- integer()
    periods <- 1:length_series
    
    if(length(changepoints) == 0){
      
      # case 1: no changepoints are found. In this case, every period belongs to the same segment
      out <- rep(1L, length_series)
      
    } else{
      
      # initialize
      out <- rep(NA_integer_, length_series)
      
      # case 2: at least one changepoint is found
      num_cp <- length(changepoints)
      
      for(j in 1:num_cp)
      {
        if(j == 1){
          idx <- which(periods < changepoints[1])
          out[idx] <- 1L
          
        } else {
          idx <- which(changepoints[j-1] <= periods & periods < changepoints[j])
          out[idx] <- j
        } # IF
      } # FOR
      
      # assign last cluster
      idx <- which(periods >= changepoints[num_cp])
      out[idx] <- num_cp + 1L
      
    } # IF 
    
    # return
    return(out)
    
} # FUN

# get all n_ij 
get_n_ij <- function(length_series, changepoints_1, changepoints_2)
{
  membership <- cbind(method1 = segmentation(length_series, changepoints_1),
                      method2 = segmentation(length_series, changepoints_2))
  
  # number of segments identified by each method
  k1 <- length(changepoints_1) + 1L
  k2 <- length(changepoints_2) + 1L
  
  # initialize output matrix
  out <- matrix(NA_integer_, nrow = k1 * k2, ncol = 3)
  colnames(out) <- c("i_in_k1", "j_in_k2", "n_ij")
  ct <- 1L
  
  for(i in 1:k1){
    for(j in 1:k2){
      n_ij <- sum(membership[,1] == i & membership[,2] == j)
      out[ct,] <- c(i, j, n_ij)
      ct <- ct + 1L
    } # FOR
  } # FOR
  
  return(out)
  
} # FUN


#' Calculates adjusted Rand index (ARI) as in Morey and Agresti (1984), adapted for changepoint estimates
#' 
#' @param changepoints_1 Vector of first estimation of changepoints
#' @param changepoints_2 Vector of first estimation of changepoints
#' @param length_series length of the series whose changepoints we the two changepoint vectors estimate
#' 
#' @return the ARI (bounded by [0, 1]); higher ARI indicates more coherence between the two partitions induced by the two changepoint vectors
#' 
#' @export
ARI <- function(changepoints_1, changepoints_2, length_series)
{
  n_ij <- get_n_ij(length_series  = length_series, 
                   changepoints_1 = changepoints_1,
                   changepoints_2 = changepoints_2)
  
  n2 <- length_series^2
  
  # number of segments
  k1 <- length(changepoints_1) + 1L 
  k2 <- length(changepoints_2) + 1L 
  
  # get TSS
  sum_nij2 <- sum(n_ij[,"n_ij"]^2)
  
  # get the marginal sums
  n_i <- sapply(1:k1, function(i) sum(n_ij[n_ij[,"i_in_k1"] == i, "n_ij"]) )
  n_j <- sapply(1:k2, function(j) sum(n_ij[n_ij[,"j_in_k2"] == j, "n_ij"]) )
  sum_ni2 <- sum(n_i^2)
  sum_nj2 <- sum(n_j^2)
  
  # get omega tilde
  prod_scaled <- sum_ni2 * sum_nj2 / n2
  numer <- sum_nij2 - prod_scaled
  denom <- 0.5 * sum_ni2 + 0.5 * sum_nj2 - prod_scaled
  
  return(numer / denom)
  
} # FOR


d1 <- function(tau, tau_hat)
{
  m1 <- length(tau)
  m2 <- length(tau_hat)
  
  # in case tau_hat is empty, return min of tau
  stopifnot(m1 > 0)
  if(m2 == 0) return(min(abs(tau)))
  
  # initialize
  distances_tauhat <- rep(NA_real_, m2)
  
  for(i in 1:m2)
  {
    distances_tau <- rep(NA_real_, m1)
    
    for(j in 1:m1)
    {
      distances_tau[j] <- abs(tau_hat[i] - tau[j])
    } # FOR
    distances_tauhat[i] <- min(distances_tau)
  } # FOR
  
  return(max(distances_tauhat))

} # FUN


d2 <- function(tau, tau_hat)
{
  m1 <- length(tau)
  m2 <- length(tau_hat)
  
  # in case tau_hat is empty, return max of tau
  stopifnot(m1 > 0)
  if(m2 == 0) return(max(abs(tau)))
  
  # initialize
  distances_tau <- rep(NA_real_, m1)
  
  for(i in 1:m1)
  {
    distances_tauhat <- rep(NA_real_, m2)
    
    for(j in 1:m2)
    {
      distances_tauhat[j] <- abs(tau[i] - tau_hat[j])
    } # FOR
    distances_tau[i] <- min(distances_tauhat)
  } # FOR
  
  return(max(distances_tau))
  
} # FUN

# length_series is the length of the original series whose changepoints we try to estimate
# function uses relative changepoints, scaling controlled via length_series

#' Calculate Hausdorff distance of two changepoint estimates
#' 
#' @param tau Integer vector of 'true' changepoints
#' @param tau_hat Integer vector of estimated changepoints
#' @param scaling A scaling factor that the vectors will be scaled with. To obtain relative distances, set equal to the length of the series whose changepoints we are interested with.
#' 
#' @return A vector of d1 (over-segmentation error), d2 (under-segmentation error), and dH (Hausdorff distance).
#' 
#' @export
hausdorff <- function(tau, tau_hat, scaling = 1.0)
{
  # normalize 
  tau     <- tau / scaling
  tau_hat <- tau_hat / scaling
  
  # calculate distances
  distance1 <- d1(tau = tau, tau_hat = tau_hat)
  distance2 <- d2(tau = tau, tau_hat = tau_hat)
  
  return(c(d1 = distance1, 
           d2 = distance2,
           dH = max(distance1, distance2)))
  
} # FOR


# ARI with adjustment of Hubert and Arabie (1985; corrected-for-chance version of the Rand index)
ARI_hubert <- function(changepoints_true, changepoints_estimate, length_series) {
  
  true <- segmentation(length_series = length_series, changepoints = changepoints_true)
  estimate <- segmentation(length_series = length_series, changepoints = changepoints_estimate)
  
  # calculation below does not work if each item is its own cluster
  seq_p <- seq_along(true)
  if (identical(unname(true), seq_p) && identical(unname(true), seq_p)) {
    return(1.0)
  }
  # compute contingency table
  n_ij <- table(true, estimate)
  a_i <- rowSums(n_ij)
  b_j <- colSums(n_ij)
  n <- sum(n_ij)
  # compute sums of binomial coefficients
  sum_n_ij_2 <- sum(choose(n_ij, 2))
  sum_a_i_2 <- sum(choose(a_i, 2))
  sum_b_j_2 <- sum(choose(b_j, 2))
  n_2 <- choose(n, 2)
  # compute adjusted Rand index
  tmp <- sum_a_i_2 * sum_b_j_2 / n_2
  (sum_n_ij_2 - tmp) / (0.5 * (sum_a_i_2 + sum_b_j_2) - tmp)
}



cp_heatmap <- function(cp_mat, 
                         respondent_order = 1:nrow(cp_mat), # how to order respondents on y-axis?
                         true_start = list(Respondent = 1:nrow(cp_mat),
                                           Item = rep(1, nrow(cp_mat))),
                         true_stop = list(Respondent = 1:nrow(cp_mat),
                                          Item = rep(ncol(cp_mat), nrow(cp_mat))),
                         legend = TRUE,
                         frame_size = 0.3){
  
  library(ggplot2)
  df    <- expand.grid(Respondent = respondent_order, Item = 1:ncol(cp_mat))
  df$cp <- as.numeric(cp_mat)
  
  plt <- ggplot(df, aes(x = Respondent, y = Item, fill = cp)) + 
    geom_tile() +
    scale_y_continuous(expand = c(0,0))+
    scale_x_continuous(expand = c(0,0)) +
    scale_fill_gradient(low = "azure", high = "red")
  
  if(!is.null(true_start)){
    frames_start <- data.frame(Respondent = true_start$Respondent,
                               Item       = true_start$Item)
    plt <- plt +
      geom_rect(data = frames_start, size = frame_size,
                fill = NA, color = "blue",
                aes(xmin = Respondent - 0.5,
                    xmax = Respondent + 0.5,
                    ymin = Item - 0.5, 
                    ymax = Item + 0.5))
    
  }
  if(!is.null(true_stop)){
    frames_stop <- data.frame(Respondent = true_stop$Respondent,
                              Item       = true_stop$Item)
    plt <- plt +
      geom_rect(data = frames_stop, size = frame_size,
                fill = NA, color = "blue",
                aes(xmin = Respondent - 0.5,
                    xmax = Respondent + 0.5,
                    ymin = Item - 0.5, 
                    ymax = Item + 0.5))
  } # IF
  
  if(!legend) plt <- plt + theme(legend.position = "none") 
  
  # flip coordinates
  plt <- plt + coord_flip()
  
  return(plt)
} # FUN


# the rows of cp_mat only hold temporarily careless observations
careless_evaluate <- function(cp_mat, start_careless, stop_careless)
{
  # initialize
  n                           <- nrow(cp_mat)
  num_items                   <- ncol(cp_mat)
  perf_arr_careless           <- matrix(0.0, n, 11L)
  colnames(perf_arr_careless) <- c(paste0("d=", seq(-3,3,1)),
                                   "d1", "d2", "dH", "ARI")
  
  # input check
  # stopifnot(n == length(start_careless) & length(start_careless) == length(stop_careless))
  
  for(i in 1:n){
    
    tau_hat <- which(cp_mat[i,] == 1L)
    m_hat   <- length(tau_hat)
    
    # difference between number of changepoints
    tau     <- c(start_careless[i], stop_careless[i])
    m_true  <- length(tau)
    diff    <- m_hat - m_true
    
    # truncate if necessary
    if(diff > 3)  diff <- 3L
    if(diff < -3) diff <- -3L
    
    # add counter on difference
    perf_arr_careless[i, diff + 4L] <- 1L
    
    # calculate and append the performance measures
    d <- hausdorff(tau = tau, tau_hat = tau_hat, scaling = num_items)
    ari <- ARI(changepoints_1 = tau, 
               changepoints_2 = tau_hat,
               length_series  = num_items)
    perf_arr_careless[i, 8:10] <- d
    perf_arr_careless[i, 11]   <- ari
    
  } # FOR i in 1:n
  
  return(perf_arr_careless)
  
} # FUN


# counts number of flagged periods under the assumption that there are no true changepoints
careful_evaluate <- function(cp_mat)
{
  # initialize
  n                           <- nrow(cp_mat)
  num_items                   <- ncol(cp_mat)
  perf_arr_careful            <- matrix(0.0, n, 3L)
  colnames(perf_arr_careful)  <- paste0("mhat=", 0:2)
  
  for(i in 1:n){
    tau_hat <- which(cp_mat[i,] == 1L)
    m_hat   <- length(tau_hat)
    
    # truncate if necessary
    if(m_hat > 2L) m_hat <- 2L  
    
    # append counter
    perf_arr_careful[i, m_hat + 1L] <- 1L 
  } # FOR
  
  return(perf_arr_careful)
} # FUN


plot_scores <- function(scores, cp_true, cp_hat, time = FALSE)
{
  # all inputs are vectors
  stopifnot(is.vector(scores))
  if(time){
    plot(scores, type = "l", ylab = "response time", xlab = "item")
  } else{
    plot(scores, type = "l", ylab = "carelessness score", xlab = "item")
  }
  
  
  # add vertical lines for true changepoints
  for(j in cp_true){
    abline(v = j, col = "red")
  } # FOR
  
  # add vertical lines for estimates changepoints
  for(j in cp_hat){
    abline(v = j, lty = "longdash", col = "blue")
  } # FOR
} # FUN


# if careful = TRUE, then all elements in start_careless will be considered empty (required for ARI in case of 'true' segmentation being only one segment). WARNING: In this case, the ARI will always be zero if the estimated segmentation does not consist of only one unique value (that is, ARI is zero if the two segmentations aren't exactly equal to each other). Thus, ARI isn't practical for evaluation of careful respondents
evaluate_endsurvey <- function(cp_mat, start_careless, careful = FALSE)
{
  # initialize
  n         <- nrow(cp_mat)
  num_items <- ncol(cp_mat)
  perf      <- matrix(0.0, n, 3L)
  colnames(perf) <- c(paste0("m=", 0:1), "ARI")
  
  for(i in 1:n){
    
    tau_hat <- which(cp_mat[i,] == 1)
    m_hat   <- length(tau_hat)
    perf[i, m_hat + 1] <- 1.0
    
    if(careful){
      tau <- integer()
    } else{
      tau <- start_careless[i] 
    } # IF
    
    # calculate ARI
    true <- segmentation(num_items, tau)
    est  <- segmentation(num_items, tau_hat)
    perf[i, 3L] <- mclust::adjustedRandIndex(true, est)
  } # FOR
  return(perf)
} # FUN


# if careful = TRUE, then all elements in start_careless will be considered empty (required for ARI in case of 'true' segmentation being only one segment). WARNING: In this case, the ARI will always be zero if the estimated segmentation does not consist of only one unique value (that is, ARI is zero if the two segmentations aren't exactly equal to each other). Thus, ARI isn't practical for evaluation of careful respondents
evaluate_endsurvey_vec <- function(flagged, start_careless, num_items, careful = FALSE)
{
  # init
  ari <- diver <- rep(NA, length(flagged))

  # loop over respondents
  for(i in seq_along(flagged)){
    
    tau_hat <- flagged[i]
    if(is.na(tau_hat))
    {
      tau_hat <- integer()
    }
    
    if(careful){
      tau <- integer()
    } else{
      tau <- start_careless[i] 
    } # IF
    
    # I'd like to count the divergence between the flagged and true changepoint. 
    # However, it's not clear how to handle the case when there is a CP, but
    # nothing was flagged. In this case, we'll assign the value "Inf"
    # Conversely, if there is no CP, but a CP was flagged, we assign "-Inf"
    # if there is no CP and no CP is flagged, we assign value zero
    isempty_tau    <- length(tau) < 1
    isempty_tauhat <- length(tau_hat) < 1
    
    if(!isempty_tau & !isempty_tauhat)
    {
      d <- tau - tau_hat
    } else if(isempty_tau & !isempty_tauhat)
    {
      d <- -Inf
    } else if(!isempty_tau & isempty_tauhat)
    {
      d <- Inf
    } else{
      d <- 0L
    }
    
    # assign divergence
    diver[i] <- d 
    
    ## calculate ARI
    # If tau_hat = NA, but tau is nonempty (false negative), then ARI = 0
    true   <- segmentation(num_items, tau)
    est    <- segmentation(num_items, tau_hat)
    ari[i] <- mclust::adjustedRandIndex(true, est)
  } # FOR
  
  ## proportions of divergence
  d0 <- diver
  
  # censor at abs(5)
  d0[!is.infinite(d0) & d0 <= -5L] <- -5L
  d0[!is.infinite(d0) & d0 >=  5L] <- 5L
  
  # get frequencies
  tbl <- table(d0)
  nam_out <- c(-Inf, seq(from = -5, to = 5, by = 1), Inf)
  tab_out <- structure(rep(0L, length(nam_out)), names = nam_out)
  tab_out[names(tbl)] <- as.integer(tbl)
  
  # rename to account for censoring
  nam_out[nam_out %in% c(-5, 5)] <- c("leq-5", "geq5")
  names(tab_out) <- nam_out
  
  ## confusion
  cnf <- confusion(contam  = which(!is.na(start_careless)), 
                   flagged = which(!is.na(flagged)), 
                   sample_size = length(flagged))

  
  return(list(ARI = ari, 
              divergence = diver, 
              prop_divergence = tab_out / length(d0),
              confusion = cnf))
} # FUN


# contam and flagged are vector of indices, sample_size the sample size (=maximum index)
confusion <- function(contam, flagged, sample_size)
{
  # suppose that among the 21 flagged, 18 are correct and 3 are wrong -> FP = 3 / 21 = FDR and power = 18 / 21 = TPR
  # in case there is no contamination, then FPR = num_flagged / sample_size
  P  <- length(contam)
  tp <- sum(contam %in% flagged)
  fn <- P - tp
  
  # of the negatives in data, how many have been flagged and
  # how many have not been flagged? (N = FP + TN)
  negatives <- setdiff(seq_len(sample_size), contam)
  N  <- length(negatives)
  fp <- sum(negatives %in% flagged)
  tn <- N - fp
  
  # FDR and FOR
  FDR <- fp / (fp + tp)
  FOR <- fn / (fn + tn)
  
  return(c(
    TPR = tp / P, TNR = tn / N, FPR = fp / N, FNR = fn / P, FDR = FDR, FOR = FOR)
  )
}


# 'flagged' is a vector of length n that holds an integer of carelessness onset if identified and NA otherwise
# returns binary matrix intended as building block for the heatmap
get_cp_mat <- function(flagged, num_items)
{
  cp_mat <- matrix(0L, length(flagged), num_items)
  onset_subjects <- which(!is.na(flagged))
  onset_location <- flagged[onset_subjects]
  
  for(s in seq_along(onset_subjects))
  {
    cp_mat[onset_subjects[s], onset_location[s]] <- 1L
  }
  return(cp_mat)
}


## classification performance
classification <- function(idx_flagged, idx_careless, n)
{
  idx_unflagged <- setdiff(seq_len(n), idx_flagged)
  positives <- idx_careless
  negatives <- setdiff(seq_len(n), positives)
  num_pos <- length(positives)
  num_neg <- length(negatives)
  
  tp <- length(intersect(idx_flagged, positives))
  fp <- length(intersect(idx_flagged, negatives))
  tn <- length(intersect(idx_unflagged, negatives))
  fn <- length(intersect(idx_unflagged, positives))
  tpr <- tp / num_pos
  fnr <- fn / num_pos
  tnr <- tn / num_neg
  fpr <- fp / num_neg
  
  return(c(TPR = tpr, FNR = fnr, TNR = tnr, FPR = fpr))
}


#' calculate average difference between true and estimated onset location, as well as FNR
#' 
#' @param onset_estimated vector of size n that holds indices of estimated carelessness onset locations and NAs at subjects without an estimated onset
#' @param onset_true vector of size n that holds indices of true carelessness onset locations and NAs at subjects that never become careless
#' @param num_items number of items
#' @param onset2_true optional vector of 2nd onsets, same structure as \code{onset_estimated}. If supplied, report the difference to the closest true onset. Note: if there is an NA in a given position in \code{onset_true}, it must also be an NA in \code{onset2_true}.
#' @note don't forget that in our code, careless_start never has NAs, so to make it fit here we need to set careless_start[idx_careful] = NA. Also, don't run this code on truly attentive respondents only
onset_divergence <- function(onset_estimated, onset_true, num_items, onset2_true = NULL)
{
  stopifnot(length(onset_estimated) == length(onset_true))
  
  ### calculate difference between estimated and true onset
  # subjects that are careless and in which a changepoint was flagged
  isflagged     <- !is.na(onset_estimated)
  idx_flagged   <- which(isflagged)
  idx_positives <- which(!is.na(onset_true))
  idx_correctlyflagged <- intersect(idx_flagged, idx_positives)
  
  
  ## if there is a second true onset but only one changepoint can be flagged (like here),
  # then we calculate the difference to the nearest true onset point
  if(!is.null(onset2_true))
  {
    stopifnot(all(which(is.na(onset_true)) == which(is.na(onset2_true))))
    stopifnot(length(onset_true) == length(onset2_true))
    differences1 <- abs(onset_estimated[idx_correctlyflagged] - onset_true[idx_correctlyflagged])
    differences2 <- abs(onset_estimated[idx_correctlyflagged] - onset2_true[idx_correctlyflagged])
    
    # difference between estimated and true onset in `idx_correctlyflagged`
    differences  <- pmin(differences1, differences2)
    
  } else
  {
    ## no 2nd true onset
    # difference between estimated and true onset in `idx_correctlyflagged`
    differences <- abs(onset_estimated[idx_correctlyflagged] - onset_true[idx_correctlyflagged])
  } # IF
  
  
  ## average and median distance
  # NOTE: this is NaN if no changepoint is flagged
  differences_mean   <- mean(differences)
  differences_median <- median(differences)
  
  
  ### false negative rate
  fn <- intersect(idx_positives, which(!isflagged)) # set of false negatives
  fnr <- length(fn) / length(idx_positives)
    
  ## organize output
  out <- c(absdiff_mean   = differences_mean,
           absdiff_median = differences_median,
           reldiff_mean   = differences_mean / num_items,
           reldiff_median = differences_median / num_items,
           fnr            = fnr)
  return(out)
} # FUN


#' calculate FNR for the benchmark methods
#' 
#' @param idx_flagged indices of flagged respondents
#' @param idx_positives indices of true careless respondents
#' @param n sample size of ALL respondents
benchmarks_fnr <- function(idx_flagged, idx_positives, n, idx_subgroup = NULL)
{
  idx_notflagged <- setdiff(seq_len(n), idx_flagged) # set of unflagged respondents
  fn <- intersect(idx_positives, idx_notflagged) # set of false negatives
  return(length(fn) / length(idx_positives))
} # FUN
