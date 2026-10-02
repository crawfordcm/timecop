#' Fit a multitask (common-plus-individual) latent VAR model
#'
#' Estimates a shared common transition matrix mu and sparse person-specific
#' deviations Delta_k from a [timecop_multitask-class] object, following the
#' Fisher-style common-plus-individual VAR. Each subject's transition matrix is
#' B_k = mu + Delta_k. Estimation uses the covariance-only proximal-gradient
#' solver ([multitask_pgd()]) on subject-specific (N_k-weighted) latent
#' cross-products, and the penalties (lambda_mu, lambda_delta) are selected over
#' a \code{n_lambda_mu x n_lambda_delta} grid using EBIC ([multitask_ebic()]).
#'
#' The estimation operates on the latent covariance matrices already stored in
#' each subject's [timecop-class] object, so the marginal type (Bernoulli,
#' Poisson, Gaussian) is accounted for through the link-function pipeline in
#' [timecop()]. Subjects are weighted by their usable observation count N_k, so
#' longer series contribute proportionally more (the common-plus-individual
#' analogue of pooled least squares).
#'
#' @param object A [timecop_multitask-class] object built using
#'   [timecop_multitask()].
#' @param n_lambda_mu Integer. Number of common-penalty values in the
#'   data-driven log-spaced grid. Ignored if \code{lambda_mu} is supplied.
#'   Default 20.
#' @param n_lambda_delta Integer. Number of deviation-penalty values in the
#'   data-driven log-spaced grid. Ignored if \code{lambda_delta} is supplied.
#'   Default 20.
#' @param lambda_mu Numeric or NULL. Custom common-penalty values. A single
#'   value fits one model at that penalty; a vector uses those values as the
#'   grid. When \code{NULL} (default) a data-driven log-spaced grid of length
#'   \code{n_lambda_mu} is used.
#' @param lambda_delta Numeric or NULL. Custom deviation-penalty values. As for
#'   \code{lambda_mu}. Default \code{NULL}.
#' @param gamma_ebic Numeric. EBIC hyperparameter between 0 and 1. Larger values
#'   favor sparser models. Default 0.5.
#' @param penalty Character. One of \code{"lasso"} (uniform weights),
#'   \code{"adaptive"} (adaptive lasso, weights fixed from an initial
#'   per-subject Yule-Walker estimate), or \code{"scad"} (SCAD via local linear
#'   approximation, weights updated across LLA iterations). Default
#'   \code{"lasso"}.
#' @param penalty_gamma Numeric. Exponent for adaptive lasso weights. Only used
#'   when \code{penalty = "adaptive"}. Default 1.
#' @param scad_a Numeric. SCAD shape parameter. Default 3.7 (Fan & Li 2001).
#'   Only used when \code{penalty = "scad"}.
#' @param adaptive_weights_from Character. Source of the pilot estimate used to
#'   build the adaptive-lasso weights (temporal and precision). Either
#'   \code{"unpenalized"} (per-subject Yule-Walker transition + inverse residual
#'   covariance) or \code{"penalized"} (per-subject EBIC-selected
#'   [fit_graphical_var()]). Only used when \code{penalty = "adaptive"}. Default
#'   \code{"unpenalized"}.
#' @param gvar Logical. If \code{TRUE}, fit the multitask graphical VAR: jointly
#'   estimate the temporal networks and common-plus-unique innovation precision
#'   networks (Omega_k = M + E_k) by alternating precision-weighted temporal PGD
#'   with a multitask graphical-lasso ADMM. The four penalties
#'   (lambda_mu, lambda_delta, lambda_M, lambda_E) are selected over the full
#'   4-D grid by joint EBIC; a supplied scalar penalty collapses its axis. Grid
#'   points are warm-started from the previous fit. Default \code{FALSE}
#'   (ordinary multitask VAR). Note: the full grid can be large, so use small
#'   \code{n_lambda_*} values.
#' @param lambda_M Numeric or NULL. Penalty on the shared off-diagonal precision
#'   network M. Only used when \code{gvar = TRUE}. When \code{NULL} (default) a
#'   log-spaced grid of length \code{n_lambda_prec} is used, anchored at the KKT
#'   zeroing scale (the smallest penalty giving a fully sparse M); when a
#'   scalar, held fixed.
#' @param lambda_E Numeric or NULL. Penalty on the task-specific off-diagonal
#'   precision deviations E_k. Only used when \code{gvar = TRUE}. When
#'   \code{NULL} (default) a log-spaced grid of length \code{n_lambda_prec} is
#'   used, anchored at its own KKT zeroing scale (M and E have different
#'   anchors); when a scalar, held fixed.
#' @param n_lambda_prec Integer. Number of penalty values per precision axis
#'   (lambda_M and lambda_E) in the 4-D EBIC grid. Only used when
#'   \code{gvar = TRUE}. Default 8.
#' @param search Character. How the 4-D penalty lattice is explored when
#'   \code{gvar = TRUE}. \code{"grid"} (default) evaluates every combination
#'   (global lattice optimum; cost is the product of the four grid sizes).
#'   \code{"coordinate"} alternates 2-D sweeps of the temporal and precision
#'   penalty pairs until the selection stops moving (a blockwise lattice
#'   optimum; cost is a few times the SUM of the two 2-D faces). Every
#'   evaluation fits the full joint model either way; only the walk through the
#'   lattice differs. With \code{"coordinate"}, unvisited entries of
#'   \code{ebic_grid} are \code{NA}.
#' @param max_sweeps Integer. Maximum number of coordinate sweeps (one sweep =
#'   a temporal face plus a precision face) before stopping. Only used when
#'   \code{gvar = TRUE} and \code{search = "coordinate"}; convergence is
#'   typically reached in 2-3 sweeps. Default 10.
#' @param rho Numeric. ADMM penalty parameter for the precision block. Only used
#'   when \code{gvar = TRUE}. Default 1.
#' @param max_outer Integer. Maximum outer alternations. Only used when
#'   \code{gvar = TRUE}. Default 50.
#' @param max_iter Integer. Maximum solver iterations per grid point. Default
#'   1000.
#' @param tol Numeric. Relative objective convergence tolerance. Default 1e-7.
#' @param verbose Logical. Show search progress: a progress bar for
#'   \code{search = "grid"}, or an in-place status line (sweep, face position,
#'   fits, running best EBIC) plus per-sweep summaries for
#'   \code{search = "coordinate"}. Default \code{FALSE}.
#'
#' @return When \code{gvar = FALSE}, an object of class
#'   \code{"timecop_multitask_fit"} (a named list) with elements:
#'   \item{mu_hat}{d x d EBIC-selected common transition matrix.}
#'   \item{delta_hat}{Length-K list of d x d selected deviation matrices.}
#'   \item{B_hat}{Length-K list of person-specific transition matrices
#'     mu + Delta_k.}
#'   \item{lambda_mu}{Selected common penalty.}
#'   \item{lambda_delta}{Selected deviation penalty.}
#'   \item{ebic}{EBIC value at the selected model.}
#'   \item{ebic_grid}{n_lambda_mu x n_lambda_delta matrix of EBIC values.}
#'   \item{lambda_mu_seq}{Common penalty sequence used.}
#'   \item{lambda_delta_seq}{Deviation penalty sequence used.}
#'   \item{sigma_hat}{Pooled innovation covariance at the selected model.}
#'   \item{gamma_ebic}{EBIC hyperparameter used.}
#'   \item{penalty}{Penalty type used.}
#'   \item{obj}{The original \code{timecop_multitask} object.}
#'
#'   When \code{gvar = TRUE}, an object of class
#'   \code{"timecop_multitask_gvar_fit"} with \code{mu_hat}, \code{delta_hat},
#'   \code{B_hat}, the precision networks \code{Omega_hat}, \code{M_hat},
#'   \code{E_hat}, the four EBIC-selected penalties, \code{ebic}, the 4-D
#'   \code{ebic_grid}, the four penalty sequences, the \code{df} of each network,
#'   \code{outer_iter}, \code{obj}, and \code{object}.
#'
#' @usage \S4method{fit_multitask}{timecop_multitask}(object, n_lambda_mu = 20,
#'   n_lambda_delta = 20, lambda_mu = NULL, lambda_delta = NULL, gamma_ebic = 0.5,
#'   penalty = "lasso", penalty_gamma = 1, scad_a = 3.7,
#'   adaptive_weights_from = "unpenalized",
#'   gvar = FALSE, lambda_M = NULL, lambda_E = NULL, n_lambda_prec = 8,
#'   search = "grid", max_sweeps = 10, rho = 1, max_outer = 50,
#'   max_iter = 1000, tol = 1e-07, verbose = FALSE)
#' @include timecopObjectClass.R
#' @export
#' @aliases fit_multitask,timecop_multitask-method

setGeneric(
  name = "fit_multitask",
  def  = function(object, ...) standardGeneric("fit_multitask")
)

setMethod(
  f = "fit_multitask",
  signature = "timecop_multitask",
  definition = function(object,
                        n_lambda_mu    = 20L,
                        n_lambda_delta = 20L,
                        lambda_mu      = NULL,
                        lambda_delta   = NULL,
                        gamma_ebic     = 0.5,
                        penalty        = "lasso",
                        penalty_gamma  = 1,
                        scad_a         = 3.7,
                        adaptive_weights_from = "unpenalized",
                        gvar           = FALSE,
                        lambda_M       = NULL,
                        lambda_E       = NULL,
                        n_lambda_prec  = 8L,
                        search         = "grid",
                        max_sweeps     = 10L,
                        rho            = 1,
                        max_outer      = 50L,
                        max_iter       = 1000L,
                        tol            = 1e-7,
                        verbose        = FALSE) {

  if (!is.numeric(n_lambda_mu) || n_lambda_mu < 1) {
    stop("'n_lambda_mu' must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(n_lambda_delta) || n_lambda_delta < 1) {
    stop("'n_lambda_delta' must be a positive integer", call. = FALSE)
  }
  if (!is.null(lambda_mu) && (!is.numeric(lambda_mu) || any(lambda_mu < 0))) {
    stop("'lambda_mu' must be a non-negative numeric scalar or vector", call. = FALSE)
  }
  if (!is.null(lambda_delta) && (!is.numeric(lambda_delta) || any(lambda_delta < 0))) {
    stop("'lambda_delta' must be a non-negative numeric scalar or vector", call. = FALSE)
  }
  if (!is.numeric(gamma_ebic) || gamma_ebic < 0 || gamma_ebic > 1) {
    stop("'gamma_ebic' must be in [0, 1]", call. = FALSE)
  }
  if (!is.character(penalty) || length(penalty) != 1 ||
      !penalty %in% c("lasso", "adaptive", "scad")) {
    stop("'penalty' must be one of 'lasso', 'adaptive', or 'scad'", call. = FALSE)
  }
  if (!is.character(adaptive_weights_from) || length(adaptive_weights_from) != 1 ||
      !adaptive_weights_from %in% c("unpenalized", "penalized")) {
    stop("'adaptive_weights_from' must be 'unpenalized' or 'penalized'", call. = FALSE)
  }
  if (!is.numeric(penalty_gamma) || penalty_gamma <= 0) {
    stop("'penalty_gamma' must be a positive number", call. = FALSE)
  }
  if (!is.numeric(scad_a) || scad_a <= 2) {
    stop("'scad_a' must be a number greater than 2", call. = FALSE)
  }
  if (!is.logical(gvar) || length(gvar) != 1) {
    stop("'gvar' must be a single logical value", call. = FALSE)
  }
  if (!is.null(lambda_M) && (!is.numeric(lambda_M) || length(lambda_M) != 1 || lambda_M < 0)) {
    stop("'lambda_M' must be NULL or a non-negative scalar", call. = FALSE)
  }
  if (!is.null(lambda_E) && (!is.numeric(lambda_E) || length(lambda_E) != 1 || lambda_E < 0)) {
    stop("'lambda_E' must be NULL or a non-negative scalar", call. = FALSE)
  }
  if (!is.numeric(rho) || length(rho) != 1 || rho <= 0) {
    stop("'rho' must be a positive scalar", call. = FALSE)
  }
  if (!is.numeric(max_outer) || max_outer < 1) {
    stop("'max_outer' must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(n_lambda_prec) || n_lambda_prec < 1) {
    stop("'n_lambda_prec' must be a positive integer", call. = FALSE)
  }
  if (!is.character(search) || length(search) != 1 ||
      !search %in% c("grid", "coordinate")) {
    stop("'search' must be 'grid' or 'coordinate'", call. = FALSE)
  }
  if (!is.numeric(max_sweeps) || max_sweeps < 1) {
    stop("'max_sweeps' must be a positive integer", call. = FALSE)
  }

  K <- object@K
  d <- object@d
  p <- object@p
  N <- object@N

  # ---- per-subject latent blocks ---------------------------------------
  # cov_z_hat[,,p+1] = lag-0 covariance (G_k role / Gamma^{--})
  # cov_z_hat[,,p]   = lag-1 Cov(z_t, z_{t-1}) (current-past); H_k = t(.)
  S0  <- lapply(object@subjects, function(s) s@cov_z_hat[,, p + 1L])
  S10 <- lapply(object@subjects, function(s) s@cov_z_hat[,, p])

  # N_k-weighted cross-products for the solver
  G <- lapply(seq_len(K), function(k) N[k] * S0[[k]])
  H <- lapply(seq_len(K), function(k) N[k] * t(S10[[k]]))

  # ---- adaptive weights (fixed once from a per-subject pilot estimate) --
  # All four weight sets (temporal W_mu/W_delta and precision W_M/W_E) come from
  # one pilot per subject, chosen by 'adaptive_weights_from'. Computed BEFORE
  # the lambda grids so the auto grids can be calibrated to the weighted
  # (adaptive) KKT conditions.
  W_mu <- NULL; W_delta <- NULL; W_M <- NULL; W_E <- NULL
  if (penalty == "adaptive") {
    aw <- multitask_adaptive_weights(object@subjects, S0, S10,
                                     source = adaptive_weights_from,
                                     gamma = penalty_gamma)
    W_mu    <- aw$W_mu
    W_delta <- aw$W_delta
    W_M     <- aw$W_M
    W_E     <- aw$W_E
  }

  # ---- lambda grids (KKT-anchored; see multitask_lambda_grids) ----------
  grids <- multitask_lambda_grids(H, S0, S10, N,
                                  n_lambda_mu = n_lambda_mu,
                                  n_lambda_delta = n_lambda_delta,
                                  n_lambda_prec = n_lambda_prec,
                                  lambda_mu = lambda_mu,
                                  lambda_delta = lambda_delta,
                                  lambda_M = lambda_M, lambda_E = lambda_E,
                                  W_mu = W_mu, W_delta = W_delta,
                                  W_M = W_M, W_E = W_E,
                                  gvar = gvar)
  lambda_mu_seq    <- grids$lambda_mu_seq
  lambda_delta_seq <- grids$lambda_delta_seq

  # ---- multitask graphical VAR: full 4-D EBIC grid ---------------------
  # Evaluate the joint EBIC over the complete (lambda_mu, lambda_delta,
  # lambda_M, lambda_E) grid. Every grid point runs the full alternating fit
  # (precision-weighted temporal PGD + precision ADMM) and is scored by the
  # joint EBIC. A supplied scalar penalty collapses that axis to one value.
  # Grid points are warm-started from the previous fit (continuation order).
  if (gvar) {

    # fail fast on non-PD joint latent covariances: if the 2d x 2d block matrix
    # [[S0, S10'], [S10, S0]] is not PSD, the innovation covariance S_eps(B)
    # can be indefinite, the Gaussian likelihood is unbounded, and both the
    # selection and the solvers degrade badly. The repair is pd_approx = TRUE
    # at construction.
    bad <- which(vapply(seq_len(K), function(k) {
      J <- rbind(cbind(S0[[k]], t(S10[[k]])), cbind(S10[[k]], S0[[k]]))
      min(eigen(J, symmetric = TRUE, only.values = TRUE)$values) < -1e-8
    }, logical(1)))
    if (length(bad)) {
      # immediate. = TRUE: by default R defers warning display until the
      # top-level call returns, which for a long grid search means the user
      # would only see this diagnosis AFTER the slow fit it warns about
      warning(sprintf(paste0(
        "[fit_multitask] the joint lag-0/lag-1 latent covariance is not positive ",
        "definite for subject(s) %s. The likelihood is unbounded at such points, ",
        "which degrades selection and slows the solvers; rebuild the object with ",
        "timecop_multitask(..., pd_approx = TRUE)."),
        paste(bad, collapse = ", ")), call. = FALSE, immediate. = TRUE)
    }

    lambda_M_seq <- grids$lambda_M_seq
    lambda_E_seq <- grids$lambda_E_seq

    nMu <- length(lambda_mu_seq); nD <- length(lambda_delta_seq)
    nMg <- length(lambda_M_seq);  nEg <- length(lambda_E_seq)
    total <- nMu * nD * nMg * nEg
    if (search == "grid" && total > 1000) {
      message(sprintf(paste0("[fit_multitask] gvar full grid has %d points ",
                             "(%d x %d x %d x %d); this may be slow. Reduce ",
                             "n_lambda_mu / n_lambda_delta / n_lambda_prec, ",
                             "or use search = \"coordinate\"."),
                      total, nMu, nD, nMg, nEg))
    }

    ebic_grid <- array(NA_real_, dim = c(nMu, nD, nMg, nEg))
    best_ebic <- Inf
    best      <- NULL
    best_idx  <- rep(NA_integer_, 4L)
    n_fit     <- 0L                            # actual model fits (cache misses)
    prev_mu <- NULL; prev_delta <- NULL       # temporal warm-start (continuation)
    prev_prec <- NULL                          # precision ADMM warm-start

    # One full-model evaluation at lattice point (i, j, a, b): complete
    # alternation (all four penalties active) scored by the joint EBIC.
    # Warm-started from the previous evaluation; skipped if already evaluated.
    eval_point <- function(i, j, a, b) {
      if (!is.na(ebic_grid[i, j, a, b])) return(invisible(NULL))
      gfit <- multitask_gvar_alternate(
        G, H, S0, S10, N,
        lambda_mu = lambda_mu_seq[i], lambda_delta = lambda_delta_seq[j],
        lambda_M = lambda_M_seq[a], lambda_E = lambda_E_seq[b],
        W_mu = W_mu, W_delta = W_delta, W_M = W_M, W_E = W_E,
        mu_init = prev_mu, delta_init = prev_delta, prec_init = prev_prec,
        rho = rho, max_outer = max_outer,
        pgd_max = max_iter, pgd_tol = tol
      )
      eb <- multitask_gvar_ebic(gfit$mu, gfit$delta, gfit$B,
                                gfit$Omega, gfit$M, gfit$E,
                                S0, S10, N, gamma_ebic = gamma_ebic)
      n_fit <<- n_fit + 1L
      ebic_grid[i, j, a, b] <<- eb$ebic
      # only warm-start from valid fits (an Inf-scoring degenerate fit would
      # contaminate the next grid point's starting state)
      if (is.finite(eb$ebic)) {
        prev_mu <<- gfit$mu; prev_delta <<- gfit$delta
        prev_prec <<- list(Omega = gfit$Omega, M = gfit$M, E = gfit$E, U = gfit$U)
      }
      if (eb$ebic < best_ebic) {
        best_ebic <<- eb$ebic
        best      <<- list(gfit = gfit, eb = eb)
        best_idx  <<- c(i, j, a, b)
      }
      invisible(NULL)
    }

    if (search == "grid") {
      # ---- exhaustive 4-D grid ------------------------------------------
      if (verbose) {
        pb   <- utils::txtProgressBar(min = 0, max = total, style = 3)
        done <- 0L
      }
      for (i in seq_len(nMu)) {
        for (j in seq_len(nD)) {
          for (a in seq_len(nMg)) {
            for (b in seq_len(nEg)) {
              eval_point(i, j, a, b)
              if (verbose) { done <- done + 1L; utils::setTxtProgressBar(pb, done) }
            }
          }
        }
      }
      if (verbose) close(pb)
    } else {
      # ---- coordinate search: alternating 2-D sweeps ---------------------
      # Every evaluation is the full joint model; only the walk through the
      # lattice is coordinate-wise. Sweep the temporal face holding the
      # precision penalties fixed, then the precision face holding the
      # temporal winners, until the selected quadruple stops moving
      # (a blockwise optimum on the lattice).
      n_sweeps <- as.integer(max_sweeps)

      # in-place status line: sweep, position within the current face, actual
      # model fits so far (cached points tick the position but not the fits),
      # and the running best EBIC
      status <- function(sweep, face, pos, tot) {
        cat(sprintf("\rsweep %d/%d | %-9s face %4d/%d | fits %4d | best EBIC %.2f   ",
                    sweep, n_sweeps, face, pos, tot, n_fit, best_ebic))
        utils::flush.console()
      }
      clear_status <- function() cat("\r", strrep(" ", 78), "\r", sep = "")

      idx <- c(1L, 1L, max(1L, ceiling(nMg / 2)), max(1L, ceiling(nEg / 2)))
      for (sweep in seq_len(n_sweeps)) {
        idx_old <- idx

        # temporal sweep at fixed (lambda_M, lambda_E)
        pos <- 0L
        for (i in seq_len(nMu)) for (j in seq_len(nD)) {
          eval_point(i, j, idx[3], idx[4])
          pos <- pos + 1L
          if (verbose) status(sweep, "temporal", pos, nMu * nD)
        }
        face <- matrix(ebic_grid[, , idx[3], idx[4]], nMu, nD)
        w <- which(face == min(face), arr.ind = TRUE)[1, ]
        idx[1] <- w[1]; idx[2] <- w[2]

        # precision sweep at fixed (lambda_mu, lambda_delta)
        pos <- 0L
        for (a in seq_len(nMg)) for (b in seq_len(nEg)) {
          eval_point(idx[1], idx[2], a, b)
          pos <- pos + 1L
          if (verbose) status(sweep, "precision", pos, nMg * nEg)
        }
        face <- matrix(ebic_grid[idx[1], idx[2], , ], nMg, nEg)
        w <- which(face == min(face), arr.ind = TRUE)[1, ]
        idx[3] <- w[1]; idx[4] <- w[2]

        if (verbose) {
          clear_status()
          cat(sprintf("  sweep %d: lambda (mu, delta, M, E) indices = (%d, %d, %d, %d), EBIC = %.2f\n",
                      sweep, idx[1], idx[2], idx[3], idx[4],
                      ebic_grid[idx[1], idx[2], idx[3], idx[4]]))
        }
        if (all(idx == idx_old)) break
      }
      if (verbose) {
        cat(sprintf("  %s after %d sweep(s), %d model fits (%d lattice points)\n",
                    if (all(idx == idx_old)) "converged" else "stopped at max sweeps",
                    sweep, n_fit, sum(is.finite(ebic_grid))))
      }
    }

    if (is.null(best)) {
      stop(paste0("no valid model found on the penalty grid (all EBIC values ",
                  "non-finite); the latent covariances are likely not jointly ",
                  "positive definite - rebuild with timecop_multitask(..., ",
                  "pd_approx = TRUE)"), call. = FALSE)
    }

    best_pen <- c(lambda_mu    = lambda_mu_seq[best_idx[1]],
                  lambda_delta = lambda_delta_seq[best_idx[2]],
                  lambda_M     = lambda_M_seq[best_idx[3]],
                  lambda_E     = lambda_E_seq[best_idx[4]])

    # a selection at either end of an AUTO grid suggests the grid does not
    # bracket the optimum (see multitask_precision_grid_diagnosis.md, Sec. 10)
    edge <- character(0)
    if (is.null(lambda_mu) && nMu > 1 &&
        best_pen[["lambda_mu"]] %in% lambda_mu_seq[c(1L, nMu)]) edge <- c(edge, "lambda_mu")
    if (is.null(lambda_delta) && nD > 1 &&
        best_pen[["lambda_delta"]] %in% lambda_delta_seq[c(1L, nD)]) edge <- c(edge, "lambda_delta")
    if (is.null(lambda_M) && nMg > 1 &&
        best_pen[["lambda_M"]] %in% lambda_M_seq[c(1L, nMg)]) edge <- c(edge, "lambda_M")
    if (is.null(lambda_E) && nEg > 1 &&
        best_pen[["lambda_E"]] %in% lambda_E_seq[c(1L, nEg)]) edge <- c(edge, "lambda_E")
    if (length(edge)) {
      message("[fit_multitask] EBIC selected ", paste(edge, collapse = ", "),
              " at the boundary of the automatic grid; the grid may not bracket ",
              "the optimum (consider a wider or finer grid).")
    }

    gfit <- best$gfit
    eb   <- best$eb
    results <- list(
      mu_hat           = gfit$mu,
      delta_hat        = gfit$delta,
      B_hat            = gfit$B,
      Omega_hat        = gfit$Omega,
      M_hat            = gfit$M,
      E_hat            = gfit$E,
      lambda_mu        = unname(best_pen["lambda_mu"]),
      lambda_delta     = unname(best_pen["lambda_delta"]),
      lambda_M         = unname(best_pen["lambda_M"]),
      lambda_E         = unname(best_pen["lambda_E"]),
      ebic             = best_ebic,
      ebic_grid        = ebic_grid,
      lambda_mu_seq    = lambda_mu_seq,
      lambda_delta_seq = lambda_delta_seq,
      lambda_M_seq     = lambda_M_seq,
      lambda_E_seq     = lambda_E_seq,
      df               = c(mu = eb$df_mu, delta = eb$df_delta,
                           M = eb$df_M, E = eb$df_E),
      gamma_ebic       = gamma_ebic,
      penalty          = penalty,
      search           = search,
      outer_iter       = gfit$outer_iter,
      obj              = gfit$obj,
      object           = object
    )
    class(results) <- c("timecop_multitask_gvar_fit", "list")
    return(results)
  }

  # ---- fit at a single (lambda_mu, lambda_delta) -----------------------
  fit_point <- function(lmu, ld, mu_ws, delta_ws) {
    if (penalty == "scad") {
      # local linear approximation: update SCAD weights from current estimate
      mu_c    <- mu_ws
      delta_c <- delta_ws
      fit     <- NULL
      for (lla in seq_len(10L)) {
        Wmu <- gvar_weights_A(mu_c, lambda = lmu, penalty = "scad", scad_a = scad_a)
        Wd  <- lapply(delta_c, function(Dk)
          gvar_weights_A(Dk, lambda = ld, penalty = "scad", scad_a = scad_a))
        fit <- multitask_pgd(G, H, lmu, ld, W_mu = Wmu, W_delta = Wd,
                             mu_init = mu_c, delta_init = delta_c,
                             max_iter = max_iter, tol = tol)
        change <- max(abs(fit$mu - mu_c),
                      max(vapply(seq_len(K),
                                 function(k) max(abs(fit$delta[[k]] - delta_c[[k]])),
                                 numeric(1))))
        mu_c    <- fit$mu
        delta_c <- fit$delta
        if (change < tol * 10) break
      }
      fit
    } else {
      multitask_pgd(G, H, lmu, ld, W_mu = W_mu, W_delta = W_delta,
                    mu_init = mu_ws, delta_init = delta_ws,
                    max_iter = max_iter, tol = tol)
    }
  }

  # ---- grid search with warm starts ------------------------------------
  nM <- length(lambda_mu_seq)
  nD <- length(lambda_delta_seq)
  ebic_grid <- matrix(NA_real_, nM, nD)

  best_ebic  <- Inf
  best       <- NULL
  best_lmu   <- NA_real_
  best_ld    <- NA_real_

  # warm-start seeds (cold at the sparse top-left of the grid)
  mu_seed    <- matrix(0, d, d)
  delta_seed <- replicate(K, matrix(0, d, d), simplify = FALSE)

  for (i in seq_len(nM)) {

    lmu       <- lambda_mu_seq[i]
    mu_ws     <- mu_seed
    delta_ws  <- delta_seed
    row_mu    <- NULL
    row_delta <- NULL

    for (j in seq_len(nD)) {

      ld <- lambda_delta_seq[j]
      if (verbose) {
        cat(sprintf("  [%2d, %2d]  lambda_mu = %.4f  lambda_delta = %.4f\n",
                    i, j, lmu, ld))
      }

      fit <- fit_point(lmu, ld, mu_ws, delta_ws)

      eb  <- multitask_ebic(fit$mu, fit$delta, fit$B, S0, S10, N,
                            gamma_ebic = gamma_ebic)
      ebic_grid[i, j] <- eb$ebic

      if (eb$ebic < best_ebic) {
        best_ebic <- eb$ebic
        best      <- list(mu = fit$mu, delta = fit$delta, B = fit$B,
                          sigma = eb$Sigma)
        best_lmu  <- lmu
        best_ld   <- ld
      }

      # warm start along the lambda_delta path
      mu_ws    <- fit$mu
      delta_ws <- fit$delta
      if (j == 1L) { row_mu <- fit$mu; row_delta <- fit$delta }
    }

    # seed the next lambda_mu row from this row's sparsest (first) fit
    mu_seed    <- row_mu
    delta_seed <- row_delta
  }

  # boundary-of-auto-grid diagnostic (as in the gvar branch)
  edge <- character(0)
  if (is.null(lambda_mu) && nM > 1 &&
      best_lmu %in% lambda_mu_seq[c(1L, nM)]) edge <- c(edge, "lambda_mu")
  if (is.null(lambda_delta) && nD > 1 &&
      best_ld %in% lambda_delta_seq[c(1L, nD)]) edge <- c(edge, "lambda_delta")
  if (length(edge)) {
    message("[fit_multitask] EBIC selected ", paste(edge, collapse = ", "),
            " at the boundary of the automatic grid; the grid may not bracket ",
            "the optimum (consider a wider or finer grid).")
  }

  results <- list(
    mu_hat           = best$mu,
    delta_hat        = best$delta,
    B_hat            = best$B,
    lambda_mu        = best_lmu,
    lambda_delta     = best_ld,
    ebic             = best_ebic,
    ebic_grid        = ebic_grid,
    lambda_mu_seq    = lambda_mu_seq,
    lambda_delta_seq = lambda_delta_seq,
    sigma_hat        = best$sigma,
    gamma_ebic       = gamma_ebic,
    penalty          = penalty,
    obj              = object
  )
  class(results) <- c("timecop_multitask_fit", "list")

  results
})
