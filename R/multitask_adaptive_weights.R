#' Adaptive penalty weights for the multitask (graphical) VAR
#'
#' Computes entry-specific adaptive-lasso weights for all four networks from a
#' per-subject pilot estimate (B_k, Omega_k). The pilot comes either from an
#' unpenalized fit (Yule-Walker transition + inverse residual covariance) or
#' from a penalized single-subject graphical VAR fit. Weights are built from the
#' common-plus-unique split of the pilots, mirroring the estimator:
#'   mu0 = mean_k B_k,                Delta0_k = B_k - mu0,
#'   M0  = mean_k offdiag(Omega_k),   E0_k     = Omega_k - M0,
#' with W = 1 / (|.| + eps)^gamma (via gvar_weights_A / gvar_weights_Omega).
#'
#' @param subjects List. The K timecop subject objects (for the penalized pilot).
#' @param S0,S10 Lists. Length-K normalized latent lag-0 and lag-1 covariances.
#' @param source Character. "unpenalized" (Yule-Walker + inverse residual) or
#'   "penalized" (per-subject fit_graphical_var). Default "unpenalized".
#' @param gamma Numeric. Adaptive-weight exponent. Default 1.
#' @param n_lambda_init Integer. Grid size for the per-subject penalized pilot
#'   fits. Only used when \code{source = "penalized"}. Default 10.
#' @return A list with \code{W_mu} (d x d), \code{W_delta} (length-K list),
#'   \code{W_M} (d x d), and \code{W_E} (length-K list).
#' @keywords internal

multitask_adaptive_weights <- function(subjects, S0, S10,
                                       source = "unpenalized", gamma = 1,
                                       n_lambda_init = 10L) {

  K <- length(S0)
  d <- nrow(S0[[1]])

  # ---- per-subject pilot estimates (B_k, Omega_k) ----------------------
  if (source == "unpenalized") {
    B <- lapply(seq_len(K), function(k) S10[[k]] %*% solve(S0[[k]]))
    Omega <- lapply(seq_len(K), function(k) {
      S_eps <- gvar_resid_cov(B[[k]], S0[[k]], S0[[k]], t(S10[[k]]), S10[[k]])
      solve(S_eps + diag(1e-3, d))
    })
  } else if (source == "penalized") {
    fits  <- lapply(seq_len(K), function(k)
      fit_graphical_var(subjects[[k]], n_lambda_A = n_lambda_init,
                        n_lambda_Omega = n_lambda_init))
    B     <- lapply(fits, function(f) f$A_hat)
    Omega <- lapply(fits, function(f) f$Omega_hat)
  } else {
    stop("'source' must be 'unpenalized' or 'penalized'", call. = FALSE)
  }

  # ---- common-plus-unique split ----------------------------------------
  mu0      <- Reduce(`+`, B) / K
  delta0   <- lapply(B, function(Bk) Bk - mu0)
  M0       <- Reduce(`+`, Omega) / K
  diag(M0) <- 0                                   # shared network is off-diagonal
  E0       <- lapply(Omega, function(Ok) Ok - M0)

  # ---- adaptive weights ------------------------------------------------
  W_mu    <- gvar_weights_A(mu0, penalty = "adaptive", gamma = gamma)
  W_delta <- lapply(delta0, function(Dk)
    gvar_weights_A(Dk, penalty = "adaptive", gamma = gamma))
  W_M     <- gvar_weights_Omega(M0, penalty = "adaptive", gamma = gamma)
  W_E     <- lapply(E0, function(Ek) {
    Ek_off <- Ek; diag(Ek_off) <- 0              # weight scale from off-diagonals
    gvar_weights_Omega(Ek_off, penalty = "adaptive", gamma = gamma)
  })

  list(W_mu = W_mu, W_delta = W_delta, W_M = W_M, W_E = W_E)
}
