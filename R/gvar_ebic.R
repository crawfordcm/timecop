#' EBIC criterion for graphical VAR model selection
#'
#' Computes the Extended Bayesian Information Criterion (EBIC) for a fitted
#' graphical VAR model, following the formulation of Foygel & Drton (2010).
#'
#' The EBIC is:
#'   EBIC = -2 * l + k * log(n) + 2 * gamma * (k_A * log(d^2) + k_Omega * log(d*(d-1)/2))
#'
#' where l is the conditional log-likelihood, k = k_A + k_Omega is the total
#' number of selected edges, k_A is the number of nonzero entries in A,
#' k_Omega is the number of nonzero off-diagonal entries in Omega (upper
#' triangle), d is the number of variables, and n is the effective sample size.
#' The coefficient is 2 (not 4) because the log terms are taken over the number
#' of candidate parameters (possible edges) rather than the number of nodes;
#' since log(possible edges) is approximately 2 * log(nodes), the two forms are
#' equivalent.
#'
#' @param A Matrix. d x d estimated temporal coefficient matrix.
#' @param Omega Matrix. d x d estimated precision matrix.
#' @param S0_plus Matrix. Lag-0 current-side covariance.
#' @param S0_minus Matrix. Lag-0 past-side covariance.
#' @param S01 Matrix. Lag-1 past-current covariance.
#' @param S10 Matrix. Lag-1 current-past covariance.
#' @param n_eff Integer. Effective sample size (time series length minus 1).
#' @param gamma_ebic Numeric. EBIC hyperparameter between 0 and 1. Default 0.5.
#' @return A list with elements \code{ebic}, \code{bic}, \code{log_lik},
#'   \code{df_A}, and \code{df_Omega}.
#' @keywords internal

gvar_ebic <- function(A, Omega, S0_plus, S0_minus, S01, S10,
                      n_eff, gamma_ebic = 0.5) {

  d     <- nrow(A)
  S_eps <- gvar_resid_cov(A, S0_plus, S0_minus, S01, S10)

  log_det <- as.numeric(determinant(Omega, logarithm = TRUE)$modulus)
  log_lik <- n_eff / 2 * (log_det - sum(S_eps * Omega))

  # degrees of freedom
  df_A     <- sum(abs(A) > 1e-8)
  df_Omega <- sum(abs(Omega[upper.tri(Omega)]) > 1e-8)
  df       <- df_A + df_Omega

  # possible edges
  p_A     <- d^2
  p_Omega <- d * (d - 1L) / 2

  bic  <- -2 * log_lik + df * log(n_eff)
  ebic <- bic + 2 * gamma_ebic * (
    df_A     * log(max(1, p_A)) +
    df_Omega * log(max(1, p_Omega))
  )

  list(ebic = ebic, bic = bic, log_lik = log_lik,
       df_A = df_A, df_Omega = df_Omega)
}
