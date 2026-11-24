#' FISTA update for the temporal coefficient matrix
#'
#' For fixed Omega, minimizes the temporal subproblem
#'   0.5 * tr(S_eps(A) * Omega) + lambda_A * sum over i,j of W_A(i,j) * |a_ij|
#' using FISTA with monotone restart and a Lipschitz step size
#'   L = lambda_max(Omega) * lambda_max(S0_minus).
#'
#' When \code{W_A = NULL} (default) all weights are 1 and the penalty reduces
#' to the standard lasso.  For adaptive lasso, pass a d x d matrix of
#' entry-specific weights.
#'
#' @param A_init Matrix. d x d warm-start for A.
#' @param Omega Matrix. d x d current innovation precision matrix.
#' @param S0_plus Matrix. Lag-0 current-side covariance.
#' @param S0_minus Matrix. Lag-0 past-side covariance.
#' @param S01 Matrix. Lag-1 past-current covariance.
#' @param S10 Matrix. Lag-1 current-past covariance.
#' @param lambda_A Numeric. Temporal penalty parameter.
#' @param W_A Matrix or NULL. d x d matrix of entry-specific penalty weights
#'   for adaptive lasso. Default \code{NULL} (uniform weights).
#' @param max_iter Integer. Maximum FISTA iterations. Default 500.
#' @param tol Numeric. Relative objective convergence tolerance. Default 1e-7.
#' @return A list with elements \code{A} (updated matrix), \code{iter}
#'   (iterations used), and \code{obj} (final subproblem objective value,
#'   smooth part plus weighted L1 penalty).
#' @keywords internal

gvar_temporal_fista <- function(A_init, Omega, S0_plus, S0_minus, S01, S10,
                                lambda_A, W_A = NULL,
                                max_iter = 500L, tol = 1e-7) {

  # soft_thresh works element-wise; tau may be a scalar or a conformable matrix
  soft_thresh <- function(x, tau) sign(x) * pmax(abs(x) - tau, 0)

  d   <- nrow(A_init)
  tau <- if (is.null(W_A)) lambda_A else lambda_A * W_A   # scalar or matrix

  # full subproblem objective: smooth part + weighted L1 penalty.
  # sum(tau * abs(Amat)) handles scalar or matrix tau identically.
  obj_full <- function(Amat) {
    gvar_objective_smooth(Amat, Omega, S0_plus, S0_minus, S01, S10) +
      sum(tau * abs(Amat))
  }

  A <- A_init
  Y <- A_init
  q <- 1

  eig_Omega  <- max(eigen(Omega, symmetric = TRUE, only.values = TRUE)$values)
  eig_S0     <- max(eigen(S0_minus, symmetric = TRUE, only.values = TRUE)$values)
  L <- eig_Omega * eig_S0
  if (L < .Machine$double.eps) L <- 1e-8

  obj_old <- obj_full(A)

  for (iter in seq_len(max_iter)) {

    # gradient of smooth part at momentum point Y: Omega(Y*S0_minus - S10)
    grad_Y    <- Omega %*% (Y %*% S0_minus - S10)
    A_new     <- soft_thresh(Y - grad_Y / L, tau / L)
    obj_new   <- obj_full(A_new)

    # monotone restart: if the full objective increased, fall back to an ISTA
    # step from A (guarantees the full objective is non-increasing)
    if (obj_new > obj_old) {
      grad_A  <- Omega %*% (A %*% S0_minus - S10)
      A_new   <- soft_thresh(A - grad_A / L, tau / L)
      obj_new <- obj_full(A_new)
      q       <- 1
      Y       <- A_new
    } else {
      q_new <- (1 + sqrt(1 + 4 * q^2)) / 2
      Y     <- A_new + ((q - 1) / q_new) * (A_new - A)
      q     <- q_new
    }

    rel_change <- abs(obj_old - obj_new) / (abs(obj_old) + 1e-12)
    A       <- A_new
    obj_old <- obj_new

    if (rel_change < tol) break
  }

  list(A = A, iter = iter, obj = obj_old)
}
