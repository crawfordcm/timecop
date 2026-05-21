#' Get Yule-Walker estimates
#'
#' @param gamma_hat Matrix. Stacked covariances matrix
#' @param Gamma_hat Matrix. Toeplitz covariance matrix
#' @param d Numeric. The number of variables
#' @return A d x d matrix of estimated VAR(1) coefficients.
#' @keywords internal

fit_var <- function(gamma_hat, Gamma_hat, d) {

  if (!is.matrix(gamma_hat)) stop("'gamma_hat' must be a matrix", call. = FALSE)
  if (!is.matrix(Gamma_hat)) stop("'Gamma_hat' must be a matrix", call. = FALSE)

  kappa_val <- kappa(Gamma_hat, exact = FALSE)
  if (kappa_val > 1e10) {
    warning(sprintf(
      "[timecop] Gamma_hat is ill-conditioned (kappa = %.2e). Yule-Walker estimates may be unreliable.",
      kappa_val
    ), call. = FALSE)
  }

  # get Yule-Walker estimates
  phi_hat <- t(solve(Gamma_hat) %*% t(gamma_hat))

  return(phi_hat)
}
