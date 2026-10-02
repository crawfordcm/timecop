#' Alternating fit for the multitask graphical VAR
#'
#' Estimates the joint common-plus-unique temporal networks (B_k = mu + Delta_k)
#' and innovation precision networks (Omega_k = M + E_k) by alternating two
#' conditional blocks at fixed penalties:
#'   - Temporal block: precision-weighted proximal gradient ([multitask_pgd()])
#'     for mu and the Delta_k, holding the Omega_k fixed.
#'   - Precision block: multitask graphical lasso ADMM
#'     ([multitask_glasso_admm()]) for Omega_k, M, E_k, holding B_k fixed.
#' The two are coupled through the innovation covariance
#'   S_eps,k(B_k) = C_k - H_k B_k' - B_k H_k' + B_k G_k B_k'  ([gvar_resid_cov()]),
#' recomputed from the current B_k before each precision update.
#'
#' The blocks \code{G}, \code{H} are the N_k-weighted cross-products used by the
#' temporal solver; \code{S0}, \code{S10} are the normalized (covariance-scale)
#' latent lag-0 and lag-1 blocks used to form S_eps,k. Under stationarity the
#' lag-0 current and past covariances coincide, so \code{S0} serves both.
#'
#' @param G,H Lists. Length-K N_k-weighted cross-product blocks for the temporal
#'   solver.
#' @param S0,S10 Lists. Length-K normalized latent lag-0 and lag-1 covariance
#'   blocks.
#' @param N Numeric. Length-K vector of usable observation counts.
#' @param lambda_mu,lambda_delta Numeric. Temporal penalties.
#' @param lambda_M,lambda_E Numeric or NULL. Precision penalties. When
#'   \code{NULL}, defaulted to 0.1 times the KKT zeroing scales of the initial
#'   innovation covariances (see [multitask_prec_anchors()]).
#' @param W_mu,W_delta Matrix/list or NULL. Adaptive temporal weights passed to
#'   the temporal solver.
#' @param W_M,W_E Matrix/list or NULL. Adaptive precision weights passed to the
#'   graphical-lasso ADMM (shared network M and unique components E_k).
#' @param mu_init,delta_init Matrix/list or NULL. Warm starts for the temporal
#'   parameters. When both supplied, the initial unweighted temporal fit is
#'   skipped and these are used to seed the alternation (useful for grid
#'   warm-starting). Default \code{NULL}.
#' @param prec_init List or NULL. Warm start for the precision ADMM state (a
#'   list with \code{Omega}, \code{M}, \code{E}, \code{U}) from a nearby fit,
#'   e.g. the previous grid point. Default \code{NULL} (cold on the first outer
#'   iteration).
#' @param rho Numeric. ADMM penalty parameter. Default 1.
#' @param max_outer Integer. Maximum outer alternations. Default 50.
#' @param tol Numeric. Outer convergence tolerance on the max parameter change.
#'   Default 1e-5.
#' @param pgd_max,pgd_tol Temporal solver iteration cap and tolerance.
#' @param admm_max Integer. Precision ADMM iteration cap. Default 1000.
#' @param verbose Logical. Print outer-iteration progress. Default FALSE.
#' @return A list with \code{mu}, \code{delta}, \code{B}, \code{Omega},
#'   \code{M}, \code{E}, the four penalties used, \code{outer_iter}, and the
#'   final penalized \code{obj}.
#' @keywords internal

multitask_gvar_alternate <- function(G, H, S0, S10, N,
                                     lambda_mu, lambda_delta,
                                     lambda_M = NULL, lambda_E = NULL,
                                     W_mu = NULL, W_delta = NULL,
                                     W_M = NULL, W_E = NULL,
                                     mu_init = NULL, delta_init = NULL,
                                     prec_init = NULL,
                                     rho = 1, max_outer = 50L, tol = 1e-5,
                                     pgd_max = 1000L, pgd_tol = 1e-7,
                                     admm_max = 1000L, verbose = FALSE) {

  K <- length(G)
  weights    <- N / 2                                   # note w_k = N_k/2
  offdiag_l1 <- function(A) sum(abs(A[row(A) != col(A)]))
  resid_cov  <- function(Bk, k)
    gvar_resid_cov(Bk, S0[[k]], S0[[k]], t(S10[[k]]), S10[[k]])

  # ---- initial temporal parameters -------------------------------------
  # warm start if supplied; otherwise an unweighted (Omega = I) temporal fit
  if (!is.null(mu_init) && !is.null(delta_init)) {
    mu    <- mu_init
    delta <- delta_init
  } else {
    tfit  <- multitask_pgd(G, H, lambda_mu, lambda_delta,
                           W_mu = W_mu, W_delta = W_delta,
                           max_iter = pgd_max, tol = pgd_tol)
    mu    <- tfit$mu
    delta <- tfit$delta
  }

  # ---- data-driven default precision penalties -------------------------
  # 0.1 x the KKT zeroing scales (which carry the likelihood weights N_k/2);
  # see multitask_prec_anchors().
  if (is.null(lambda_M) || is.null(lambda_E)) {
    B0    <- lapply(seq_len(K), function(k) mu + delta[[k]])
    Seps0 <- lapply(seq_len(K), function(k) resid_cov(B0[[k]], k))
    anch  <- multitask_prec_anchors(Seps0, weights, W_M = W_M, W_E = W_E)
    if (is.null(lambda_M)) lambda_M <- 0.1 * anch$lam_M_max
    if (is.null(lambda_E)) lambda_E <- 0.1 * anch$lam_E_max
  }

  # ---- outer alternation ------------------------------------------------
  # prec carries the ADMM state; each solve warm-starts from the previous one
  # (prec = NULL on the first iteration, unless seeded via prec_init).
  prec <- prec_init
  for (outer in seq_len(max_outer)) {

    mu_prev    <- mu
    delta_prev <- delta

    # current transitions and innovation covariances
    B     <- lapply(seq_len(K), function(k) mu + delta[[k]])
    S_eps <- lapply(seq_len(K), function(k) resid_cov(B[[k]], k))

    # precision block (multitask graphical lasso ADMM), warm-started
    prec <- multitask_glasso_admm(S_eps, lambda_M, lambda_E,
                                  weights = weights, W_M = W_M, W_E = W_E,
                                  Omega_init = prec$Omega, M_init = prec$M,
                                  E_init = prec$E, U_init = prec$U,
                                  rho = rho, max_iter = admm_max)

    # temporal block (precision-weighted PGD, warm-started)
    tfit  <- multitask_pgd(G, H, lambda_mu, lambda_delta,
                           W_mu = W_mu, W_delta = W_delta,
                           mu_init = mu, delta_init = delta,
                           Omega = prec$Omega,
                           max_iter = pgd_max, tol = pgd_tol)
    mu    <- tfit$mu
    delta <- tfit$delta

    change <- max(abs(mu - mu_prev),
                  vapply(seq_len(K),
                         function(k) max(abs(delta[[k]] - delta_prev[[k]])),
                         numeric(1)))
    if (verbose) cat(sprintf("  outer %2d: max param change = %.3e\n", outer, change))
    if (change < tol) break
  }

  # ---- final objective (note section 5) --------------------------------
  B     <- lapply(seq_len(K), function(k) mu + delta[[k]])
  S_eps <- lapply(seq_len(K), function(k) resid_cov(B[[k]], k))
  obj <- 0
  for (k in seq_len(K)) {
    ldet <- as.numeric(determinant(prec$Omega[[k]], logarithm = TRUE)$modulus)
    obj  <- obj + weights[k] * (-ldet + sum(S_eps[[k]] * prec$Omega[[k]]))
  }
  obj <- obj +
    lambda_mu    * sum(abs(mu)) +
    lambda_delta * sum(vapply(delta, function(D) sum(abs(D)), numeric(1))) +
    lambda_M     * offdiag_l1(prec$M) +
    lambda_E     * sum(vapply(prec$E, offdiag_l1, numeric(1)))

  list(mu = mu, delta = delta, B = B,
       Omega = prec$Omega, M = prec$M, E = prec$E, U = prec$U,
       lambda_mu = lambda_mu, lambda_delta = lambda_delta,
       lambda_M = lambda_M, lambda_E = lambda_E,
       outer_iter = outer, obj = obj)
}
