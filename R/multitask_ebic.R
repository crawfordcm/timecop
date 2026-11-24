#' EBIC criterion for multitask (common-plus-individual) VAR model selection
#'
#' Computes the Extended Bayesian Information Criterion (EBIC, Foygel & Drton
#' 2010) for a fitted multitask VAR, used to select the penalties
#' (lambda_mu, lambda_delta) over a grid.
#'
#' The multitask objective is plain least squares (no innovation precision
#' matrix), so the innovation covariance enters only here, through the Gaussian
#' log-likelihood. A single pooled covariance is used:
#'   Sigma = (1/N) sum_k N_k * S_eps_k,    N = sum_k N_k,
#' where S_eps_k is subject k's model-implied residual covariance at B_k. The
#' pooled-covariance conditional log-likelihood then reduces to
#'   l = -N/2 * (log det Sigma + d).
#'
#' The EBIC is
#'   EBIC = -2*l + df*log(N)
#'          + 2*gamma * (df_mu*log(d^2) + df_delta*log(K*d^2)),
#' where df_mu is the number of nonzero entries in mu, df_delta is the total
#' number of nonzero entries across all Delta_k, and df = df_mu + df_delta. The
#' log terms are over the number of candidate parameters: d^2 for mu and K*d^2
#' for the stacked deviations.
#'
#' @param mu Matrix. d x d estimated common transition matrix.
#' @param delta List. Length-K list of d x d estimated deviations Delta_k.
#' @param B List. Length-K list of person-specific transition matrices
#'   B_k = mu + Delta_k.
#' @param S0 List. Length-K list of d x d lag-0 latent covariances (one per
#'   subject); used for both the past and current side under stationarity.
#' @param S10 List. Length-K list of d x d lag-1 current-past latent
#'   cross-covariances (one per subject).
#' @param N Numeric. Length-K vector of usable observation counts (n_k - p).
#' @param gamma_ebic Numeric. EBIC hyperparameter between 0 and 1. Default 0.5.
#' @param tol_zero Numeric. Magnitude below which a coefficient is treated as
#'   zero when counting degrees of freedom. Default 1e-8.
#' @return A list with elements \code{ebic}, \code{bic}, \code{log_lik},
#'   \code{df_mu}, \code{df_delta}, and \code{Sigma} (the pooled covariance).
#' @keywords internal

multitask_ebic <- function(mu, delta, B, S0, S10, N,
                           gamma_ebic = 0.5, tol_zero = 1e-8) {

  K <- length(B)
  d <- nrow(mu)
  N_total <- sum(N)

  # pooled residual covariance: (1/N) sum_k N_k * S_eps_k
  Sigma <- matrix(0, d, d)
  for (k in seq_len(K)) {
    S0k  <- S0[[k]]
    S10k <- S10[[k]]
    S01k <- t(S10k)
    S_eps_k <- gvar_resid_cov(B[[k]], S0k, S0k, S01k, S10k)
    Sigma   <- Sigma + N[k] * S_eps_k
  }
  Sigma <- Sigma / N_total

  # degrees of freedom
  df_mu    <- sum(abs(mu) > tol_zero)
  df_delta <- sum(vapply(delta, function(D) sum(abs(D) > tol_zero), numeric(1)))
  df       <- df_mu + df_delta

  # validity guard: the Gaussian log-likelihood is meaningful only for a PD
  # pooled covariance (determinant()$modulus silently drops the sign for
  # indefinite matrices). An invalid model scores Inf so it cleanly loses
  # every grid comparison.
  if (min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    return(list(ebic = Inf, bic = Inf, log_lik = -Inf,
                df_mu = df_mu, df_delta = df_delta, Sigma = Sigma))
  }

  log_det <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
  log_lik <- -N_total / 2 * (log_det + d)

  # candidate parameters
  p_mu    <- d^2
  p_delta <- K * d^2

  bic  <- -2 * log_lik + df * log(N_total)
  ebic <- bic + 2 * gamma_ebic * (
    df_mu    * log(max(1, p_mu)) +
    df_delta * log(max(1, p_delta))
  )
  if (!is.finite(ebic)) ebic <- Inf

  list(ebic = ebic, bic = bic, log_lik = log_lik,
       df_mu = df_mu, df_delta = df_delta, Sigma = Sigma)
}
