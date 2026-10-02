#' Joint EBIC for the multitask graphical VAR
#'
#' Computes the Extended Bayesian Information Criterion (Foygel & Drton 2010)
#' for a fitted multitask graphical VAR, scoring the complete model: both the
#' temporal networks (mu, Delta_k) and the innovation precision networks
#' (M, E_k). The Gaussian conditional log-likelihood uses the estimated
#' precision matrices Omega_k and the model-implied innovation covariances
#' S_eps,k(B_k):
#'   l = sum_k (N_k / 2) ( log det Omega_k - tr(S_eps,k Omega_k) ).
#'
#' Degrees of freedom count nonzeros in all four networks: all entries of mu and
#' the Delta_k (temporal edges include the autoregressive diagonal) and the
#' off-diagonal entries of M and the E_k (contemporaneous network edges). The
#' EBIC is
#'   EBIC = -2 l + df log(N)
#'          + 2 gamma ( df_mu   log(d^2)      + df_delta log(K d^2)
#'                    + df_M    log(d(d-1)/2) + df_E     log(K d(d-1)/2) ),
#' with N = sum_k N_k.
#'
#' @param mu Matrix. d x d common transition matrix.
#' @param delta List. Length-K list of d x d temporal deviations.
#' @param B List. Length-K list of transition matrices mu + Delta_k.
#' @param Omega List. Length-K list of estimated precision matrices.
#' @param M Matrix. Shared off-diagonal precision network.
#' @param E List. Length-K list of unique precision components.
#' @param S0 List. Length-K normalized latent lag-0 covariance blocks.
#' @param S10 List. Length-K normalized latent lag-1 covariance blocks.
#' @param N Numeric. Length-K vector of usable observation counts.
#' @param gamma_ebic Numeric. EBIC hyperparameter between 0 and 1. Default 0.5.
#' @param tol_zero Numeric. Magnitude below which a coefficient counts as zero.
#'   Default 1e-8.
#' @return A list with \code{ebic}, \code{bic}, \code{log_lik}, and the four
#'   degrees of freedom \code{df_mu}, \code{df_delta}, \code{df_M}, \code{df_E}.
#' @keywords internal

multitask_gvar_ebic <- function(mu, delta, B, Omega, M, E, S0, S10, N,
                                gamma_ebic = 0.5, tol_zero = 1e-8) {

  K <- length(B)
  d <- nrow(mu)
  N_total <- sum(N)

  # degrees of freedom (temporal: all entries; precision: off-diagonal edges)
  df_mu    <- sum(abs(mu) > tol_zero)
  df_delta <- sum(vapply(delta, function(D) sum(abs(D) > tol_zero), numeric(1)))
  df_M     <- sum(abs(M[upper.tri(M)]) > tol_zero)
  df_E     <- sum(vapply(E, function(Ek) sum(abs(Ek[upper.tri(Ek)]) > tol_zero),
                         numeric(1)))
  df <- df_mu + df_delta + df_M + df_E

  # log-likelihood from estimated precisions and residual covariances, with a
  # validity guard: a non-PSD residual covariance makes the Gaussian likelihood
  # unbounded (finite here, but meaninglessly large), so such models score Inf
  # and cleanly lose every grid comparison.
  log_lik <- 0
  for (k in seq_len(K)) {
    S_eps <- gvar_resid_cov(B[[k]], S0[[k]], S0[[k]], t(S10[[k]]), S10[[k]])
    if (min(eigen(S_eps, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
      return(list(ebic = Inf, bic = Inf, log_lik = -Inf,
                  df_mu = df_mu, df_delta = df_delta, df_M = df_M, df_E = df_E))
    }
    ldet  <- as.numeric(determinant(Omega[[k]], logarithm = TRUE)$modulus)
    log_lik <- log_lik + (N[k] / 2) * (ldet - sum(S_eps * Omega[[k]]))
  }

  # candidate parameter counts
  p_mu    <- d^2
  p_delta <- K * d^2
  p_M     <- d * (d - 1) / 2
  p_E     <- K * d * (d - 1) / 2

  bic  <- -2 * log_lik + df * log(N_total)
  ebic <- bic + 2 * gamma_ebic * (
    df_mu    * log(max(1, p_mu)) +
    df_delta * log(max(1, p_delta)) +
    df_M     * log(max(1, p_M)) +
    df_E     * log(max(1, p_E))
  )
  if (!is.finite(ebic)) ebic <- Inf

  list(ebic = ebic, bic = bic, log_lik = log_lik,
       df_mu = df_mu, df_delta = df_delta, df_M = df_M, df_E = df_E)
}
