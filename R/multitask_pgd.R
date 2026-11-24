#' Proximal-gradient solver for the multitask (common-plus-individual) VAR
#'
#' Minimizes the stacked Fisher-style objective
#'   sum_k 0.5 * tr(C_k' G_k C_k) - tr(C_k' H_k)
#'     + lambda_mu * sum |W_mu * mu| + lambda_delta * sum_k |W_delta_k * Delta_k|
#' where C_k = mu + Delta_k, using FISTA (accelerated proximal gradient) with
#' monotone restart. Working in the transposed parameter blocks
#' Bt_0 = mu', Bt_k = Delta_k' (transposes do not matter for an entrywise L1
#' penalty), the smooth gradient blocks are
#'   g_k = G_k (Bt_0 + Bt_k) - H_k,    g_0 = sum_k g_k.
#'
#' When per-task innovation precision matrices \code{Omega} are supplied (the
#' multitask graphical VAR case), each gradient block is right-multiplied by
#' Omega_k, giving the precision-weighted gradient
#'   g_k = (G_k (Bt_0 + Bt_k) - H_k) Omega_k.
#' With \code{Omega = NULL} (default) this reduces to the ordinary (Omega_k = I)
#' least-squares gradient.
#'
#' The blocks \code{G} and \code{H} are passed in already weighted (e.g.
#' G_k = N_k times the lag-0 latent covariance and H_k = N_k times the lag-1
#' latent covariance, for subject-specific weighting); this routine is agnostic
#' to how they were scaled.
#'
#' The stacked Gram matrix Q = Z'Z is never assembled: the gradient uses the
#' block formula above, and the Lipschitz constant L = lambda_max(Q) is
#' estimated by power iteration on the operator Bt -> Q Bt.
#'
#' @param G List. Length-K list of d x d predictor Gram blocks (G_k role).
#' @param H List. Length-K list of d x d predictor-outcome blocks (H_k role).
#' @param lambda_mu Numeric. Penalty on the common matrix mu.
#' @param lambda_delta Numeric. Penalty on the individual deviations Delta_k.
#' @param W_mu Matrix or NULL. d x d entry-specific weights for mu (adaptive
#'   lasso), in natural mu orientation. Default \code{NULL} (uniform weights).
#' @param W_delta List or NULL. Length-K list of d x d entry-specific weight
#'   matrices for the Delta_k, in natural Delta_k orientation. Default
#'   \code{NULL} (uniform weights).
#' @param mu_init Matrix or NULL. d x d warm start for mu, in natural mu
#'   orientation. Default zeros.
#' @param delta_init List or NULL. Length-K list of d x d warm starts for the
#'   Delta_k, in natural Delta_k orientation. Default zeros.
#' @param Omega List or NULL. Length-K list of d x d innovation precision
#'   matrices for precision-weighted (multitask graphical VAR) fitting. When
#'   \code{NULL} (default) the gradient is unweighted (Omega_k = I).
#' @param max_iter Integer. Maximum FISTA iterations. Default 1000.
#' @param tol Numeric. Relative objective convergence tolerance. Default 1e-7.
#' @return A list with elements \code{mu} (d x d common matrix), \code{delta}
#'   (length-K list of d x d deviations), \code{B} (length-K list of
#'   person-specific transition matrices mu + Delta_k), \code{iter}, \code{obj}
#'   (final penalized objective), and \code{L} (Lipschitz constant used).
#' @keywords internal

multitask_pgd <- function(G, H,
                          lambda_mu, lambda_delta,
                          W_mu = NULL, W_delta = NULL,
                          mu_init = NULL, delta_init = NULL,
                          Omega = NULL,
                          max_iter = 1000L, tol = 1e-7) {

  K <- length(G)
  d <- nrow(G[[1]]) #everyone should have same d

  weighted <- !is.null(Omega)   # precision-weighted (graphical VAR) gradient?

  soft_thresh <- function(x, tau) sign(x) * pmax(abs(x) - tau, 0)

  # get indices for stacked matrix
  idx <- function(k) (k * d + 1L):((k + 1L) * d)

  # per-block thresholds
  tau <- matrix(0, (K + 1L) * d, d) #initialize tau

  tau[idx(0L), ] <- if (is.null(W_mu)) lambda_mu else lambda_mu * t(W_mu) #common thresholds

  for (k in seq_len(K)) {
    tau[idx(k), ] <- if (is.null(W_delta)) lambda_delta else lambda_delta * t(W_delta[[k]]) #ind thresholds
  }

  # smooth gradient blocks (Q Bt - R), returned stacked
  grad_stacked <- function(Bt) {
    Bt0 <- Bt[idx(0L), , drop = FALSE]
    out <- matrix(0, (K + 1L) * d, d)
    g0  <- matrix(0, d, d)
    for (k in seq_len(K)) {
      Ck <- Bt0 + Bt[idx(k), , drop = FALSE]
      gk <- G[[k]] %*% Ck - H[[k]]
      if (weighted) gk <- gk %*% Omega[[k]]   # precision weighting
      out[idx(k), ] <- gk
      g0 <- g0 + gk
    }
    out[idx(0L), ] <- g0
    out
  }

  # linear operator Bt -> Q Bt (gradient with R = 0); used for power iter
  apply_Q <- function(V) {
    V0  <- V[idx(0L), , drop = FALSE]
    out <- matrix(0, (K + 1L) * d, d)
    q0  <- matrix(0, d, d)
    for (k in seq_len(K)) {
      qk <- G[[k]] %*% (V0 + V[idx(k), , drop = FALSE])
      if (weighted) qk <- qk %*% Omega[[k]]   # precision weighting
      out[idx(k), ] <- qk
      q0 <- q0 + qk
    }
    out[idx(0L), ] <- q0
    out
  }

  # smooth objective: sum_k 0.5 tr(Ck' G_k Ck Omega_k) - tr(Ck' H_k Omega_k)
  f_smooth <- function(Bt) {
    Bt0 <- Bt[idx(0L), , drop = FALSE]
    val <- 0
    for (k in seq_len(K)) {
      Ck  <- Bt0 + Bt[idx(k), , drop = FALSE] #total matrix per k
      if (weighted) {
        val <- val + 0.5 * sum((G[[k]] %*% Ck) * (Ck %*% Omega[[k]])) -
          sum(Ck * (H[[k]] %*% Omega[[k]]))
      } else {
        val <- val + 0.5 * sum((G[[k]] %*% Ck) * Ck) - sum(H[[k]] * Ck)
      }
    }
    val
  }

  obj_full <- function(Bt) f_smooth(Bt) + sum(tau * abs(Bt))

  # ---- Lipschitz constant L = lambda_max(Q) via power iteration ----------
  V <- matrix(stats::rnorm((K + 1L) * d * d), (K + 1L) * d, d)
  nv <- sqrt(sum(V^2)); if (nv > 0) V <- V / nv
  L  <- 1e-8
  for (i in seq_len(100L)) {
    W   <- apply_Q(V)
    nw  <- sqrt(sum(W^2))
    if (nw < .Machine$double.eps) break
    V_new <- W / nw
    L_new <- sum(V_new * apply_Q(V_new))
    if (abs(L_new - L) <= 1e-6 * (abs(L) + 1e-12)) { L <- L_new; V <- V_new; break }
    L <- L_new
    V <- V_new
  }
  if (L < .Machine$double.eps) L <- 1e-8

  # ---- initialize stacked parameter -------------------------------------
  # warm starts
  Bt <- matrix(0, (K + 1L) * d, d)
  if (!is.null(mu_init))    Bt[idx(0L), ] <- t(mu_init)
  if (!is.null(delta_init)) {
    for (k in seq_len(K)) Bt[idx(k), ] <- t(delta_init[[k]])
  }

  # ---- FISTA with monotone restart --------------------------------------
  Y <- Bt
  q <- 1
  obj_old <- obj_full(Bt)

  for (iter in seq_len(max_iter)) {

    grad_Y  <- grad_stacked(Y) #get gradient
    Bt_new  <- soft_thresh(Y - grad_Y / L, tau / L) #gradient step plus soft-thresholding
    obj_new <- obj_full(Bt_new) #update objective

    if (obj_new > obj_old) {
      # monotone restart: ISTA step from current iterate
      grad_B  <- grad_stacked(Bt)
      Bt_new  <- soft_thresh(Bt - grad_B / L, tau / L)
      obj_new <- obj_full(Bt_new)
      q       <- 1
      Y       <- Bt_new
    } else {
      q_new <- (1 + sqrt(1 + 4 * q^2)) / 2
      Y     <- Bt_new + ((q - 1) / q_new) * (Bt_new - Bt)
      q     <- q_new
    }

    rel_change <- abs(obj_old - obj_new) / (abs(obj_old) + 1e-12)
    Bt         <- Bt_new
    obj_old    <- obj_new

    if (rel_change < tol) break
  }

  # ---- unpack (Bt_0 = mu', Bt_k = Delta_k') -----------------------------
  mu    <- t(Bt[idx(0L), , drop = FALSE]) #common
  delta <- lapply(seq_len(K), function(k) t(Bt[idx(k), , drop = FALSE])) #ind
  B     <- lapply(delta, function(Dk) mu + Dk) #total

  list(mu = mu, delta = delta, B = B, iter = iter, obj = obj_old, L = L)
}
