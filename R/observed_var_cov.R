#' Compute observed covariances
#'
#' @param data Matrix. An n (time points) by d (variables) multivariate time series matrix
#' @param d Numeric. The number of variables
#' @param p Numeric. The VAR order. Default is 1
#' @param n Numeric. Time series length
#' @param corr Logical. Correlations or covariances
#' @return A d x d x (2p+1) array of observed covariance (or correlation) matrices
#'   at lags -p, ..., 0, ..., p.
#' @keywords internal

observed_var_cov <- function(data, d, p, n, corr) {

  cov_x_hat <- array(NA, dim = c(d, d, 2*p + 1))
  for (i in seq_len(d)) {
    for (j in seq_len(d)) {
      for (h in seq_len(2*p + 1)) {

        # set lag
        lag <- h - (p + 1)
        abs_lag <- abs(lag)

        # negative lags
        if (lag < 0) {
          idx1 <- (abs_lag + 1):n
          idx2 <- 1:(n - abs_lag)

        # positive lags
        } else if (lag > 0) {
          idx1 <- 1:(n - abs_lag)
          idx2 <- (abs_lag + 1):n

        # lag 0
        } else {
          idx1 <- idx2 <- 1:n
        }

        if (corr == FALSE) {
          cov_x_hat[i, j, h] <- cov(data[i, idx1], data[j, idx2]) * (n - abs_lag - 1) / (n - abs_lag)
        } else {
          cov_x_hat[i, j, h] <- cor(data[i, idx1], data[j, idx2]) * (n - abs_lag - 1) / (n - abs_lag) # check correlations
        }

      }
    }
  }
  return(cov_x_hat)
}



