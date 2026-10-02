#' Graphical lasso precision update for graphical VAR
#'
#' For fixed A, updates Omega by running graphical lasso on the residual
#' covariance S_eps(A).
#'
#' When \code{W_Omega = NULL} (default) a uniform penalty is applied to all
#' off-diagonal entries.  For adaptive glasso, pass a d x d symmetric matrix
#' of entry-specific weights; the diagonal is ignored regardless.
#'
#' @param A Matrix. d x d current temporal coefficient matrix.
#' @param S0_plus Matrix. Lag-0 current-side covariance.
#' @param S0_minus Matrix. Lag-0 past-side covariance.
#' @param S01 Matrix. Lag-1 past-current covariance.
#' @param S10 Matrix. Lag-1 current-past covariance.
#' @param lambda_Omega Numeric. Off-diagonal precision penalty parameter.
#' @param W_Omega Matrix or NULL. d x d symmetric matrix of entry-specific
#'   penalty weights for adaptive glasso. Default \code{NULL} (uniform weights).
#' @return A list with elements \code{Omega} (d x d sparse precision matrix)
#'   and \code{Sigma} (d x d corresponding covariance matrix).
#' @keywords internal

gvar_precision_update <- function(A, S0_plus, S0_minus, S01, S10,
                                  lambda_Omega, W_Omega = NULL) {

  S_eps <- gvar_resid_cov(A, S0_plus, S0_minus, S01, S10)

  # Ensure S_eps is positive definite before passing to glasso; the alternating
  # iterations can push A to values that make S_eps indefinite
  min_eig <- min(eigen(S_eps, symmetric = TRUE, only.values = TRUE)$values)
  if (min_eig < 1e-6) {
    S_eps <- S_eps + (1e-6 - min_eig) * diag(nrow(S_eps))
  }

  # rho may be a scalar (standard glasso) or a weight matrix (adaptive glasso);
  # glasso ignores the diagonal when penalize.diagonal = FALSE
  rho <- if (is.null(W_Omega)) lambda_Omega else lambda_Omega * W_Omega

  fit <- glasso::glasso(s = S_eps, rho = rho,
                        penalize.diagonal = FALSE, thr = 1e-8)

  list(Omega = fit$wi, Sigma = fit$w)
}
