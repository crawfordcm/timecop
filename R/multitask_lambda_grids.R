#' KKT anchors for the multitask precision penalties
#'
#' Computes the smallest penalties at which the fully sparse precision solution
#' (diagonal Omega_k, i.e. all off-diagonals of M and the E_k zero) satisfies
#' the first-order (KKT) conditions of the weighted multitask graphical-lasso
#' objective
#'   sum_k weights_k (-log det Omega_k + tr(S_eps_k Omega_k))
#'     + lambda_M ||M||_1off + lambda_E sum_k ||E_k||_1off.
#' At a diagonal Omega_k the smooth gradient for edge (i, j) of E_k is
#' weights_k * (S_eps_k)_ij, and for the shared edge M_ij it sums over tasks,
#' so zero is optimal iff
#'   |weights_k * offdiag(S_eps_k)| <= lambda_E   and
#'   |offdiag(sum_k weights_k * S_eps_k)| <= lambda_M.
#' The returned maxima are therefore the natural TOPS of the penalty grids
#' (the exact analogue of max|X'y| in lasso regression). See
#' multitask_precision_grid_diagnosis.md for the full derivation.
#'
#' Under ADAPTIVE penalties the condition for edge (i, j) becomes
#' |gradient_ij| <= lambda * W_ij, i.e. lambda >= gradient_ij / W_ij, so each
#' gradient is divided elementwise by its penalty weight before taking the
#' maximum. With W = NULL this reduces to the uniform-lasso anchors.
#'
#' @param S_eps List. Length-K list of d x d innovation covariance matrices.
#' @param weights Numeric. Length-K likelihood weights (N_k / 2).
#' @param W_M Matrix or NULL. d x d adaptive weights for the shared network M
#'   (off-diagonal). Default \code{NULL} (uniform penalty).
#' @param W_E List or NULL. Length-K list of d x d adaptive weights for the
#'   E_k (off-diagonal). Default \code{NULL} (uniform penalty).
#' @return A list with \code{lam_M_max} and \code{lam_E_max}.
#' @keywords internal

multitask_prec_anchors <- function(S_eps, weights, W_M = NULL, W_E = NULL) {
  # max over off-diagonals of |A| / W (W = 1 when NULL)
  offmax_w <- function(A, W) {
    off <- row(A) != col(A)
    if (is.null(W)) max(abs(A[off])) else max(abs(A[off]) / W[off])
  }
  K <- length(S_eps)
  lam_E_max <- max(vapply(seq_len(K),
                          function(k) offmax_w(weights[k] * S_eps[[k]],
                                               if (is.null(W_E)) NULL else W_E[[k]]),
                          numeric(1)))
  lam_M_max <- offmax_w(Reduce(`+`, Map(`*`, weights, S_eps)), W_M)
  list(lam_M_max = lam_M_max, lam_E_max = lam_E_max)
}

#' Penalty grids for the multitask (graphical) VAR
#'
#' Builds the data-driven log-spaced penalty sequences for all multitask
#' penalty axes, each anchored at its KKT zeroing threshold (the smallest
#' penalty at which the corresponding block is fully sparse), descending three
#' decades:
#'   lambda_mu    : max |sum_k H_k|          (gradient of mu at the origin)
#'   lambda_delta : max_k max |H_k|          (gradient of Delta_k at the origin)
#'   lambda_M     : max |offdiag(sum_k w_k S_eps_k)|   (w_k = N_k / 2)
#'   lambda_E     : max_k w_k max |offdiag(S_eps_k)|
#' The precision anchors are computed from the Yule-Walker innovation
#' covariances via [multitask_prec_anchors()], and only when needed. A supplied
#' penalty (scalar or vector) bypasses its auto grid and is used as the
#' (decreasing) candidate sequence for that axis.
#'
#' Under ADAPTIVE penalties (any of \code{W_mu}, \code{W_delta}, \code{W_M},
#' \code{W_E} non-NULL) the corresponding gradients are divided elementwise by
#' the weights before taking the maxima, because the weighted KKT condition for
#' entry (i, j) is |gradient_ij| <= lambda * W_ij. Without this correction the
#' grid brackets the uniform-lasso problem rather than the (down-shifted,
#' stretched) adaptive one, and EBIC selections pile up at the grid boundary.
#'
#' @param H List. Length-K N_k-weighted lag-1 cross-product blocks
#'   (H_k = N_k t(S10_k)), as used by the temporal solver.
#' @param S0,S10 Lists. Length-K normalized latent lag-0 and lag-1 covariance
#'   blocks. Only used when \code{gvar = TRUE} and a precision penalty is NULL.
#' @param N Numeric. Length-K vector of usable observation counts.
#' @param n_lambda_mu,n_lambda_delta,n_lambda_prec Integers. Grid sizes per
#'   axis (the two precision axes share \code{n_lambda_prec}).
#' @param lambda_mu,lambda_delta,lambda_M,lambda_E Numeric or NULL. Supplied
#'   penalty values (used as-is, sorted decreasing) or NULL for the auto grid.
#' @param W_mu,W_delta,W_M,W_E Adaptive penalty weights (matrix / length-K
#'   lists) or NULL for uniform penalties. Used only to calibrate the auto
#'   grids; supplied penalties are never modified.
#' @param gvar Logical. Build the precision grids too. Default FALSE.
#' @return A list with \code{lambda_mu_seq}, \code{lambda_delta_seq}, and (when
#'   \code{gvar = TRUE}) \code{lambda_M_seq}, \code{lambda_E_seq} (else NULL).
#' @keywords internal

multitask_lambda_grids <- function(H, S0, S10, N,
                                   n_lambda_mu, n_lambda_delta, n_lambda_prec,
                                   lambda_mu = NULL, lambda_delta = NULL,
                                   lambda_M = NULL, lambda_E = NULL,
                                   W_mu = NULL, W_delta = NULL,
                                   W_M = NULL, W_E = NULL,
                                   gvar = FALSE) {

  logseq <- function(top, n) exp(seq(log(top), log(top * 1e-3), length.out = n))
  # max of |A| / W over all entries (W = 1 when NULL)
  allmax_w <- function(A, W) if (is.null(W)) max(abs(A)) else max(abs(A) / W)

  # ---- temporal grids ----------------------------------------------------
  # gradient of mu at the origin is -sum_k H_k; of Delta_k, -H_k. Adaptive
  # weights divide entrywise (weighted KKT: |g_ij| <= lambda * W_ij). H lives
  # in the solver's transposed space while the W_* are natural-orientation
  # (multitask_pgd applies t(W) internally), so the division uses t(W).
  lambda_mu_seq <- if (is.null(lambda_mu)) {
    logseq(allmax_w(Reduce(`+`, H), if (is.null(W_mu)) NULL else t(W_mu)),
           n_lambda_mu)
  } else sort(lambda_mu, decreasing = TRUE)

  lambda_delta_seq <- if (is.null(lambda_delta)) {
    top <- max(vapply(seq_along(H), function(k)
      allmax_w(H[[k]], if (is.null(W_delta)) NULL else t(W_delta[[k]])),
      numeric(1)))
    logseq(top, n_lambda_delta)
  } else sort(lambda_delta, decreasing = TRUE)

  # ---- precision grids (gvar only) ---------------------------------------
  lambda_M_seq <- NULL
  lambda_E_seq <- NULL
  if (gvar) {
    if (is.null(lambda_M) || is.null(lambda_E)) {
      K     <- length(S0)
      B_yw  <- lapply(seq_len(K), function(k) S10[[k]] %*% solve(S0[[k]]))
      S_eps <- lapply(seq_len(K), function(k)
        gvar_resid_cov(B_yw[[k]], S0[[k]], S0[[k]], t(S10[[k]]), S10[[k]]))
      anch  <- multitask_prec_anchors(S_eps, N / 2, W_M = W_M, W_E = W_E)
    }
    lambda_M_seq <- if (is.null(lambda_M)) {
      logseq(anch$lam_M_max, n_lambda_prec)
    } else sort(lambda_M, decreasing = TRUE)
    lambda_E_seq <- if (is.null(lambda_E)) {
      logseq(anch$lam_E_max, n_lambda_prec)
    } else sort(lambda_E, decreasing = TRUE)
  }

  list(lambda_mu_seq = lambda_mu_seq, lambda_delta_seq = lambda_delta_seq,
       lambda_M_seq = lambda_M_seq, lambda_E_seq = lambda_E_seq)
}
