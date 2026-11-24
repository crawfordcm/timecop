#' Fit a latent graphical VAR model
#'
#' Estimates sparse temporal (A) and contemporaneous precision (Omega) networks
#' from a [timecop-class] object using the MRCE (Multivariate Regression with
#' Covariance Estimation) objective in covariance-block form. The temporal
#' update uses FISTA (proximal gradient) and the precision update uses the
#' graphical lasso. Tuning parameters are selected over a
#' \code{n_lambda_A x n_lambda_Omega} grid using EBIC.
#'
#' The estimation operates entirely on the latent covariance matrices already
#' stored in the [timecop-class] object, so the marginal type (Bernoulli,
#' Poisson, Gaussian) is accounted for through the link function pipeline in
#' [timecop()].
#'
#' @param object A [timecop-class] object built using [timecop()].
#' @param n_lambda_A Integer. Number of temporal penalty values in the
#'   data-driven log-spaced grid. Ignored if \code{lambda_A} is supplied.
#'   Default 20.
#' @param n_lambda_Omega Integer. Number of precision penalty values in the
#'   data-driven log-spaced grid. Ignored if \code{lambda_Omega} is supplied.
#'   Default 20.
#' @param lambda_A Numeric or NULL. Custom temporal penalty values. A single
#'   value fits one model at that penalty; a vector uses those values as the
#'   grid. When \code{NULL} (default) a data-driven log-spaced grid of length
#'   \code{n_lambda_A} is used.
#' @param lambda_Omega Numeric or NULL. Custom precision penalty values. A
#'   single value fits one model at that penalty; a vector uses those values
#'   as the grid. When \code{NULL} (default) a data-driven log-spaced grid of
#'   length \code{n_lambda_Omega} is used.
#' @param gamma_ebic Numeric. EBIC hyperparameter between 0 and 1. Larger values
#'   favor sparser networks. Default 0.5.
#' @param max_iter Integer. Maximum outer alternating iterations at each grid
#'   point. Default 200.
#' @param tol Numeric. Relative objective convergence tolerance for the outer
#'   loop. Default 1e-6.
#' @param penalty Character. Penalty type for the temporal and contemporaneous
#'   networks. One of \code{"lasso"} (standard, uniform weights),
#'   \code{"adaptive"} (adaptive lasso / adaptive glasso — weights fixed from
#'   an initial Yule-Walker estimate), or \code{"scad"} (SCAD via local linear
#'   approximation — weights updated at each alternating iteration from the
#'   current estimate). Default \code{"lasso"}.
#' @param penalty_gamma Numeric. Exponent for adaptive lasso weights. Commonly
#'   1 or 2. Only used when \code{penalty = "adaptive"}. Default 1.
#' @param scad_a Numeric. Shape parameter for the SCAD penalty. Default 3.7
#'   (Fan & Li 2001). Only used when \code{penalty = "scad"}.
#' @param verbose Logical. Print progress through the lambda grid. Default
#'   \code{FALSE}.
#'
#' @return An object of class \code{"timecop_gvar"} (a named list) with
#'   elements:
#'   \item{A_hat}{d x d EBIC-selected temporal coefficient matrix.}
#'   \item{Omega_hat}{d x d EBIC-selected innovation precision matrix.}
#'   \item{lambda_A}{Selected temporal penalty.}
#'   \item{lambda_Omega}{Selected precision penalty.}
#'   \item{ebic}{EBIC value at the selected model.}
#'   \item{ebic_grid}{n_lambda_A x n_lambda_Omega matrix of EBIC values.}
#'   \item{lambda_A_seq}{Temporal penalty sequence used.}
#'   \item{lambda_Omega_seq}{Precision penalty sequence used.}
#'   \item{gamma_ebic}{EBIC hyperparameter used.}
#'   \item{penalty}{Penalty type used (\code{"lasso"}, \code{"adaptive"}, or
#'     \code{"scad"}).}
#'   \item{obj}{The original \code{timecop} object.}
#'
#' @examples
#' \donttest{
#' sim <- latent_var_sim(
#'   d = 3, n = 200, p = 1,
#'   param = list(0.5, 0.5, 0.5),
#'   phi_lv = matrix(c(0.4, 0.1, 0.0,
#'                     0.0, 0.3, 0.1,
#'                     0.1, 0.0, 0.4), 3, 3, byrow = TRUE),
#'   family = list("Bernoulli", "Bernoulli", "Bernoulli"),
#'   omega_lv = diag(3)
#' )
#' obj  <- timecop(data = t(sim$X_t), family = list("Bernoulli", "Bernoulli", "Bernoulli"))
#' gvar <- fit_graphical_var(obj, n_lambda_A = 10, n_lambda_Omega = 10)
#' gvar$A_hat
#' gvar$Omega_hat
#' }
#'
#' @usage \S4method{fit_graphical_var}{timecop}(object, n_lambda_A = 20,
#'   n_lambda_Omega = 20, lambda_A = NULL, lambda_Omega = NULL, gamma_ebic = 0.5,
#'   penalty = "lasso", penalty_gamma = 1, scad_a = 3.7,
#'   max_iter = 200, tol = 1e-06, verbose = FALSE)
#' @importFrom glasso glasso
#' @include timecopObjectClass.R
#' @export
#' @aliases fit_graphical_var,timecop-method

setGeneric(
  name = "fit_graphical_var",
  def  = function(object, ...) standardGeneric("fit_graphical_var")
)

setMethod(
  f = "fit_graphical_var",
  signature = "timecop",
  definition = function(object,
                        n_lambda_A     = 20L,
                        n_lambda_Omega = 20L,
                        lambda_A       = NULL,
                        lambda_Omega   = NULL,
                        gamma_ebic     = 0.5,
                        penalty        = "lasso",
                        penalty_gamma  = 1,
                        scad_a         = 3.7,
                        max_iter       = 200L,
                        tol            = 1e-6,
                        verbose        = FALSE) {

  if (!is.numeric(n_lambda_A) || n_lambda_A < 1) {
    stop("'n_lambda_A' must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(n_lambda_Omega) || n_lambda_Omega < 1) {
    stop("'n_lambda_Omega' must be a positive integer", call. = FALSE)
  }
  if (!is.null(lambda_A) && (!is.numeric(lambda_A) || any(lambda_A < 0))) {
    stop("'lambda_A' must be a non-negative numeric scalar or vector", call. = FALSE)
  }
  if (!is.null(lambda_Omega) && (!is.numeric(lambda_Omega) || any(lambda_Omega < 0))) {
    stop("'lambda_Omega' must be a non-negative numeric scalar or vector", call. = FALSE)
  }
  if (!is.numeric(gamma_ebic) || gamma_ebic < 0 || gamma_ebic > 1) {
    stop("'gamma_ebic' must be in [0, 1]", call. = FALSE)
  }
  if (!is.character(penalty) || length(penalty) != 1 ||
      !penalty %in% c("lasso", "adaptive", "scad")) {
    stop("'penalty' must be one of 'lasso', 'adaptive', or 'scad'", call. = FALSE)
  }
  if (!is.numeric(penalty_gamma) || penalty_gamma <= 0) {
    stop("'penalty_gamma' must be a positive number", call. = FALSE)
  }
  if (!is.numeric(scad_a) || scad_a <= 2) {
    stop("'scad_a' must be a number greater than 2", call. = FALSE)
  }

  d     <- object@d
  n_eff <- object@n - 1L   # T-1 usable observation pairs

  # ---- Extract latent covariance blocks --------------------------------
  # cov_z_hat[,,p]   = S10 = Cov(z_t, z_{t-1})   (lag -1, index p)
  # cov_z_hat[,,p+1] = S0  = Cov(z_t, z_t)        (lag 0,  index p+1)
  p       <- object@p
  S0      <- object@cov_z_hat[,, p + 1L]
  S10     <- object@cov_z_hat[,, p]
  S01     <- t(S10)
  S0_minus <- S0
  S0_plus  <- S0       # equal under stationarity

  # ---- Lambda grids ----------------------------------------------------
  A_yw       <- S10 %*% solve(S0_minus)
  S_eps_init <- gvar_resid_cov(A_yw, S0_plus, S0_minus, S01, S10)

  if (is.null(lambda_A)) {
    lambda_A_max <- max(abs(S10))
    lambda_A_seq <- exp(seq(log(lambda_A_max),
                            log(lambda_A_max * 1e-3),
                            length.out = n_lambda_A))
  } else {
    lambda_A_seq <- sort(lambda_A, decreasing = TRUE)
  }

  if (is.null(lambda_Omega)) {
    lambda_Omega_max <- max(abs(S_eps_init[upper.tri(S_eps_init)]))
    lambda_Omega_seq <- exp(seq(log(lambda_Omega_max),
                                log(lambda_Omega_max * 1e-3),
                                length.out = n_lambda_Omega))
  } else {
    lambda_Omega_seq <- sort(lambda_Omega, decreasing = TRUE)
  }

  # ---- Adaptive weights (computed once from initial estimates) ---------
  # For "lasso":    W_A = W_Omega = NULL throughout (no weighting)
  # For "adaptive": compute once here from Yule-Walker / unpenalized Omega
  # For "scad":     initialize to NULL (uniform); updated inside the loop
  W_A     <- NULL
  W_Omega <- NULL

  if (penalty == "adaptive") {
    S_eps_init <- gvar_resid_cov(A_yw, S0_plus, S0_minus, S01, S10)
    Omega_init <- tryCatch(solve(S_eps_init), error = function(e) diag(d))

    W_A     <- gvar_weights_A(A_yw, penalty = "adaptive", gamma = penalty_gamma)
    W_Omega <- gvar_weights_Omega(Omega_init, penalty = "adaptive", gamma = penalty_gamma)
  }

  # ---- Grid search with warm starts ------------------------------------
  nA <- length(lambda_A_seq)
  nO <- length(lambda_Omega_seq)
  ebic_grid <- matrix(NA_real_, nA, nO)

  best_ebic  <- Inf
  best_A     <- matrix(0, d, d)
  best_Omega <- diag(d)
  best_lA    <- NA_real_
  best_lO    <- NA_real_

  # warm-start seeds
  A_row_start     <- A_yw
  Omega_row_start <- diag(d)

  for (i in seq_len(nA)) {

    lA    <- lambda_A_seq[i]
    A_i   <- A_row_start
    Omega_i <- Omega_row_start

    for (j in seq_len(nO)) {

      lO <- lambda_Omega_seq[j]

      if (verbose) {
        cat(sprintf("  [%2d, %2d]  lambda_A = %.4f  lambda_Omega = %.4f\n",
                    i, j, lA, lO))
      }

      # Alternating MRCE
      obj_prev <- gvar_objective(A_i, Omega_i, S0_plus, S0_minus,
                                  S01, S10, lA, lO,
                                  penalty = penalty, W_A = W_A,
                                  W_Omega = W_Omega, scad_a = scad_a)
      for (it in seq_len(max_iter)) {

        # For SCAD: recompute weights from current estimates before each step
        if (penalty == "scad") {
          W_Omega <- gvar_weights_Omega(Omega_i, lambda = lO,
                                        penalty = "scad", scad_a = scad_a)
          W_A     <- gvar_weights_A(A_i, lambda = lA,
                                    penalty = "scad", scad_a = scad_a)
        }

        # Step 1: update Omega given A
        prec_res <- gvar_precision_update(A_i, S0_plus, S0_minus, S01, S10,
                                           lO, W_Omega = W_Omega)
        Omega_i  <- prec_res$Omega

        # Step 2: update A given Omega
        fista_res <- gvar_temporal_fista(A_i, Omega_i, S0_plus, S0_minus,
                                          S01, S10, lA, W_A = W_A,
                                          max_iter = 500L, tol = tol * 0.1)
        A_i <- fista_res$A

        obj_curr   <- gvar_objective(A_i, Omega_i, S0_plus, S0_minus,
                                     S01, S10, lA, lO,
                                     penalty = penalty, W_A = W_A,
                                     W_Omega = W_Omega, scad_a = scad_a)
        rel_change <- abs(obj_prev - obj_curr) / (abs(obj_prev) + 1e-12)
        obj_prev   <- obj_curr
        if (rel_change < tol) break
      }

      # EBIC at this grid point
      eb <- gvar_ebic(A_i, Omega_i, S0_plus, S0_minus, S01, S10,
                      n_eff, gamma_ebic)
      ebic_grid[i, j] <- eb$ebic

      if (eb$ebic < best_ebic) {
        best_ebic  <- eb$ebic
        best_A     <- A_i
        best_Omega <- Omega_i
        best_lA    <- lA
        best_lO    <- lO
      }

      # warm start along lambda_Omega path (Omega only; A stays)
      Omega_row_start <- Omega_i
    }

    # warm start A for next lambda_A row
    A_row_start     <- A_i
    Omega_row_start <- diag(d)   # reset Omega warm start for new lambda_A row
  }

  results <- list(
    A_hat            = best_A,
    Omega_hat        = best_Omega,
    lambda_A         = best_lA,
    lambda_Omega     = best_lO,
    ebic             = best_ebic,
    ebic_grid        = ebic_grid,
    lambda_A_seq     = lambda_A_seq,
    lambda_Omega_seq = lambda_Omega_seq,
    gamma_ebic       = gamma_ebic,
    penalty          = penalty,
    obj              = object
  )
  class(results) <- c("timecop_gvar", "list")

  results
})
