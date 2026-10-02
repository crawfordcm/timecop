#' Residual covariance for graphical VAR
#'
#' Computes S_eps(A) = S0_plus - A*S01 - S10*t(A) + A*S0_minus*t(A).
#'
#' @param A Matrix. d x d temporal coefficient matrix.
#' @param S0_plus Matrix. Lag-0 covariance of current-side variables.
#' @param S0_minus Matrix. Lag-0 covariance of past-side variables.
#' @param S01 Matrix. Lag-1 past-current cross-covariance.
#' @param S10 Matrix. Lag-1 current-past cross-covariance.
#' @return A d x d symmetric residual covariance matrix.
#' @keywords internal

gvar_resid_cov <- function(A, S0_plus, S0_minus, S01, S10) {
  S_eps <- S0_plus - A %*% S01 - S10 %*% t(A) + A %*% S0_minus %*% t(A)
  (S_eps + t(S_eps)) / 2   # symmetrize against numerical drift
}
