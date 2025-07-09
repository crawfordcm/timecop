#' Get Yule-Walker estimates
#'
#' @param gamma_hat Matrix. Stacked covariances matrix
#' @param Gamma_hat Matrix. Toeplitz covariance matrix
#' @param d Numeric. The number of variables

fit_var <- function(gamma_hat, Gamma_hat, d) {

  # get Yule-Walker estimates
  gamma_hat_Z <- as.vector(gamma_hat)

  phi_hat <- matrix(t(as.matrix(bdiag(replicate(d, solve(Gamma_hat), simplify = FALSE)) %*% matrix(gamma_hat_Z, d^2, 1))), nrow = d)

  return(phi_hat)
}
