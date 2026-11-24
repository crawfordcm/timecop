#' Per-task positive-definite precision update (multitask graphical lasso)
#'
#' Solves, for one task, over positive-definite Omega,
#'   weight * (-log det Omega + tr(S Omega)) + rho/2 * ||Omega - A||_F^2,
#' which has a closed form via eigendecomposition. Writing
#' C = A - (weight/rho) S = Q diag(c) Q', the solution shares the eigenvectors Q
#' and has eigenvalues (c_j + sqrt(c_j^2 + 4*weight/rho)) / 2, all strictly
#' positive, so the returned matrix is symmetric positive definite.
#'
#' @param S Matrix. d x d innovation covariance for the task.
#' @param A Matrix. d x d ADMM target (M + E_k - U_k).
#' @param weight Numeric. Likelihood weight (e.g. N_k / 2).
#' @param rho Numeric. ADMM penalty parameter.
#' @return A d x d symmetric positive-definite precision matrix.
#' @keywords internal

multitask_omega_update <- function(S, A, weight, rho) {
  C   <- (A + t(A)) / 2 - (weight / rho) * S
  C   <- (C + t(C)) / 2
  eig <- eigen(C, symmetric = TRUE)
  vals <- (eig$values + sqrt(eig$values^2 + 4 * weight / rho)) / 2
  Omega <- eig$vectors %*% (vals * t(eig$vectors))   # Q diag(vals) Q'
  (Omega + t(Omega)) / 2
}

#' Edgewise common-plus-unique precision update (M, E_k)
#'
#' Solves the decomposition step of the multitask precision ADMM. With
#' diag(M) = 0, each task's diagonal is absorbed into E_k, and each
#' off-diagonal edge (i, j) is an independent scalar common-plus-unique lasso
#'   lambda_M |m| + lambda_E sum_k |e_k| + rho/2 sum_k (a_k - m - e_k)^2,
#' solved by coordinate descent (soft is elementwise soft-thresholding):
#'   e_k <- soft(a_k - m, lambda_E/rho),
#'   m   <- soft( mean_k(a_k - e_k), lambda_M/(rho K) ).
#'
#' @param R_list List. Length-K list of d x d ADMM targets (Omega_k + U_k).
#' @param M_start Matrix. d x d warm start for the common component.
#' @param E_start List. Length-K list of d x d warm starts for the unique
#'   components.
#' @param lambda_M Numeric. Penalty on the shared off-diagonal network M.
#' @param lambda_E Numeric. Penalty on the task-specific off-diagonal deviations.
#' @param rho Numeric. ADMM penalty parameter.
#' @param W_M Matrix or NULL. d x d adaptive weights for M (off-diagonal). When
#'   \code{NULL} (default) the M threshold is uniform.
#' @param W_E List or NULL. Length-K list of d x d adaptive weights for the E_k
#'   (off-diagonal). When \code{NULL} (default) the E threshold is uniform.
#' @param inner_max Integer. Maximum per-edge coordinate-descent iterations.
#'   Default 100.
#' @param inner_tol Numeric. Per-edge convergence tolerance. Default 1e-8.
#' @return A list with elements \code{M} (d x d, diag 0) and \code{E} (length-K
#'   list of d x d, carrying the per-task diagonals).
#' @keywords internal

multitask_common_unique_update <- function(R_list, M_start, E_start,
                                           lambda_M, lambda_E, rho,
                                           W_M = NULL, W_E = NULL,
                                           inner_max = 100L, inner_tol = 1e-8) {

  K <- length(R_list)
  d <- nrow(R_list[[1]])
  soft <- function(x, tau) sign(x) * pmax(abs(x) - tau, 0)

  # diag(M) = 0; each E_k absorbs the diagonal of its target R_k
  M_new <- matrix(0, d, d)
  E_new <- lapply(seq_len(K), function(k) {
    Ek <- matrix(0, d, d); diag(Ek) <- diag(R_list[[k]]); Ek
  })

  if (d >= 2) {
    for (i in seq_len(d - 1L)) {
      for (j in (i + 1L):d) {
        a <- vapply(R_list,  function(R) R[i, j], numeric(1))
        m <- M_start[i, j]
        e <- vapply(E_start, function(E) E[i, j], numeric(1))

        # entry-specific thresholds (adaptive) or uniform (W = NULL)
        tau_e <- if (is.null(W_E)) rep(lambda_E / rho, K) else
          vapply(seq_len(K), function(k) lambda_E * W_E[[k]][i, j] / rho, numeric(1))
        tau_m <- (if (is.null(W_M)) lambda_M else lambda_M * W_M[i, j]) / (rho * K)

        for (inner in seq_len(inner_max)) {
          e_new  <- soft(a - m, tau_e)
          m_new  <- soft(mean(a - e_new), tau_m)
          change <- max(abs(c(m_new - m, e_new - e)))
          m <- m_new
          e <- e_new
          if (change < inner_tol) break
        }

        M_new[i, j] <- M_new[j, i] <- m
        for (k in seq_len(K)) E_new[[k]][i, j] <- E_new[[k]][j, i] <- e[k]
      }
    }
  }

  list(M = M_new, E = E_new)
}

#' Multitask common-plus-unique graphical lasso via ADMM
#'
#' Estimates per-task innovation precision matrices Omega_k that share a common
#' off-diagonal network and have task-specific deviations:
#'   Omega_k = M + E_k,    diag(M) = 0,    Omega_k positive definite.
#' Minimizes
#'   sum_k weight_k (-log det Omega_k + tr(S_k Omega_k))
#'     + lambda_M * offdiagL1(M) + lambda_E * sum_k offdiagL1(E_k)
#' (offdiagL1 = sum of absolute off-diagonal entries)
#' by ADMM: a per-task positive-definite precision update
#' ([multitask_omega_update()]), an edgewise common-plus-unique decomposition
#' ([multitask_common_unique_update()]), and a scaled dual update.
#'
#' @param S_list List. Length-K list of d x d innovation covariance matrices.
#' @param lambda_M Numeric. Penalty on the shared off-diagonal network M.
#' @param lambda_E Numeric. Penalty on the task-specific off-diagonal deviations.
#' @param weights Numeric or NULL. Length-K positive likelihood weights (e.g.
#'   N_k / 2). Default \code{NULL} (all ones).
#' @param W_M Matrix or NULL. d x d adaptive weights for the shared network M.
#'   Default \code{NULL} (uniform penalty).
#' @param W_E List or NULL. Length-K list of d x d adaptive weights for the
#'   unique components E_k. Default \code{NULL} (uniform penalty).
#' @param Omega_init,M_init,E_init,U_init Warm starts for the ADMM state (the
#'   precision matrices, shared network, unique components, and scaled duals).
#'   When \code{NULL} (default) the corresponding variable is cold-initialized.
#'   Supplying the state from a nearby solve (e.g. the previous outer iteration)
#'   greatly reduces the number of ADMM iterations.
#' @param rho Numeric. ADMM penalty parameter. Default 1.
#' @param max_iter Integer. Maximum ADMM iterations. Default 1000.
#' @param abs_tol Numeric. Absolute tolerance for the stopping rule. Default
#'   1e-5.
#' @param rel_tol Numeric. Relative tolerance for the stopping rule. Default
#'   1e-4.
#' @param inner_max Integer. Maximum per-edge coordinate-descent iterations in
#'   the decomposition step. Default 100.
#' @param inner_tol Numeric. Per-edge convergence tolerance. Default 1e-8.
#' @return A list with elements \code{Omega} (length-K list of precision
#'   matrices), \code{M} (shared off-diagonal network), \code{E} (length-K list
#'   of unique components), \code{U} (length-K list of scaled duals, for
#'   warm-starting), \code{converged}, \code{iterations}, \code{lambda_M},
#'   \code{lambda_E}, \code{weights}, and \code{rho}.
#' @keywords internal

multitask_glasso_admm <- function(S_list, lambda_M, lambda_E,
                                  weights = NULL, W_M = NULL, W_E = NULL,
                                  Omega_init = NULL, M_init = NULL,
                                  E_init = NULL, U_init = NULL,
                                  rho = 1,
                                  max_iter = 1000L, abs_tol = 1e-5,
                                  rel_tol = 1e-4, inner_max = 100L,
                                  inner_tol = 1e-8) {

  if (!is.list(S_list) || length(S_list) < 1) {
    stop("'S_list' must be a non-empty list of covariance matrices", call. = FALSE)
  }
  if (lambda_M < 0 || lambda_E < 0) {
    stop("penalty parameters must be non-negative", call. = FALSE)
  }
  if (rho <= 0) stop("'rho' must be positive", call. = FALSE)

  sym <- function(A) (A + t(A)) / 2
  K <- length(S_list)
  d <- nrow(S_list[[1]])
  S_list <- lapply(S_list, sym)

  if (is.null(weights)) weights <- rep(1, K)
  if (length(weights) != K || any(weights <= 0)) {
    stop("'weights' must be positive and have length K", call. = FALSE)
  }

  # initialization (warm start when supplied, else cold)
  M     <- if (is.null(M_init)) matrix(0, d, d) else M_init
  E     <- if (is.null(E_init)) lapply(S_list, function(S) diag(1 / pmax(diag(S), 1e-8), d)) else E_init
  Omega <- if (is.null(Omega_init)) E else Omega_init
  U     <- if (is.null(U_init)) lapply(seq_len(K), function(k) matrix(0, d, d)) else U_init

  converged <- FALSE
  for (iter in seq_len(max_iter)) {

    M_old <- M
    E_old <- E

    # 1. per-task precision updates
    for (k in seq_len(K)) {
      A_k        <- M + E[[k]] - U[[k]]
      Omega[[k]] <- multitask_omega_update(S_list[[k]], A_k, weights[k], rho)
    }

    # 2. common-plus-unique decomposition
    R_list <- lapply(seq_len(K), function(k) Omega[[k]] + U[[k]])
    dec <- multitask_common_unique_update(R_list, M, E, lambda_M, lambda_E, rho,
                                          W_M = W_M, W_E = W_E,
                                          inner_max = inner_max, inner_tol = inner_tol)
    M <- dec$M
    E <- dec$E

    # 3. scaled dual updates
    for (k in seq_len(K)) U[[k]] <- U[[k]] + Omega[[k]] - M - E[[k]]

    # convergence (primal/dual residuals)
    primal <- sqrt(sum(vapply(seq_len(K),
      function(k) sum((Omega[[k]] - M - E[[k]])^2), numeric(1))))
    dual <- rho * sqrt(sum(vapply(seq_len(K),
      function(k) sum(((M + E[[k]]) - (M_old + E_old[[k]]))^2), numeric(1))))

    omega_norm <- sqrt(sum(vapply(Omega, function(O) sum(O^2), numeric(1))))
    decomp_norm <- sqrt(sum(vapply(seq_len(K),
      function(k) sum((M + E[[k]])^2), numeric(1))))
    u_norm <- sqrt(sum(vapply(U, function(Uk) sum(Uk^2), numeric(1))))

    eps_primal <- sqrt(K * d^2) * abs_tol + rel_tol * max(omega_norm, decomp_norm)
    eps_dual   <- sqrt(K * d^2) * abs_tol + rel_tol * rho * u_norm

    if (primal <= eps_primal && dual <= eps_dual) {
      converged <- TRUE
      break
    }
  }

  list(Omega = Omega, M = M, E = E, U = U,
       converged = converged, iterations = iter,
       lambda_M = lambda_M, lambda_E = lambda_E,
       weights = weights, rho = rho)
}
