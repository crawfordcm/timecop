#' Check if covariance matrices are PD
#'
#' @param gamma_hat Matrix. Stacked covariance matrix
#' @param Gamma_hat Matrix. Toeplitz covariance matrix
#' @param d Numeric. The number of variables
#' @param eig.tol Numeric. Eigen tolerance
#' @return A list of two matrices: the (possibly adjusted) `gamma_hat` and
#'   `Gamma_hat`.
#' @importFrom matrixcalc is.positive.definite
#' @importFrom Matrix nearPD
#' @keywords internal

check_pd <- function(gamma_hat, Gamma_hat, d, eig.tol) {

  # build Toeplitz matrix
  cov_toep <- rbind(
    cbind(Gamma_hat, t(gamma_hat)),
    cbind(gamma_hat, Gamma_hat)
  )

  not_pd <- !is.positive.definite(cov_toep)

  if (not_pd) {
    posd.tol <- 1e-7
    pd_mat   <- nearPD(cov_toep, corr = TRUE, keepDiag = TRUE, eig.tol = eig.tol, posd.tol = posd.tol)
    cov_toep <- as.matrix(pd_mat$mat)

    gamma_hat <- cov_toep[(d+1):(2*d), 1:d]
    Gamma_hat <- cov_toep[1:d, 1:d]
  }

  return(list(gamma_hat, Gamma_hat))

}
