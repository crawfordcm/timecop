#' Repair the joint latent covariance to be positive definite
#'
#' Projects the joint lag-0/lag-1 latent covariance matrix onto the
#' intersection of (i) the positive-definite cone and (ii) the stationary
#' block-Toeplitz structure
#'   J = rbind(cbind(S0, t(S10)), cbind(S10, S0))
#' (equal diagonal blocks, transpose-consistent off-diagonal blocks). Joint
#' positive definiteness of J guarantees a positive semidefinite innovation
#' covariance S_eps(B) = (I, -B) J (I, -B)' for EVERY transition matrix B, so a
#' single repair at construction covers all downstream estimation.
#'
#' A single nearPD() call is not sufficient: it projects onto the PSD cone
#' WITHOUT the structure constraint, so its output generally has two different
#' diagonal blocks. Extracting one of them and reassembling the structured
#' matrix can land outside the PSD cone again (this was the bug that let
#' pd_approx = TRUE return non-PD joint matrices). The alternation below
#' projects onto the PSD cone (nearPD), then onto the structure set (average
#' the diagonal blocks; average the off-diagonal blocks with their
#' transposes), and repeats until the STRUCTURED matrix itself is positive
#' definite; both sets are convex. If the alternation stalls, a fallback
#' shrinks the lag-1 block toward zero, whose block-diagonal limit is positive
#' definite once S0 is.
#'
#' With \code{correlation = TRUE} (the latent scale is unit-variance by
#' construction) the diagonal of S0 is set to exactly 1 before any check, so
#' even the "no repair needed" path standardizes the diagonal. Assumes VAR
#' order p = 1 (a two-block Toeplitz structure).
#'
#' @param S10 Matrix. d x d lag-1 latent covariance, Cov(z_t, z_{t-1}).
#' @param S0 Matrix. d x d lag-0 latent covariance (symmetrized on entry).
#' @param d Numeric. The number of variables.
#' @param eig.tol Numeric. Eigenvalue tolerance passed to nearPD. Default 1e-6.
#' @param pd.tol Numeric. Minimum eigenvalue above which a matrix counts as
#'   positive definite. Default 1e-8.
#' @param conv.tol Numeric. Relative-change tolerance for detecting a stalled
#'   alternation. Default 1e-8.
#' @param maxit Integer. Maximum alternating-projection iterations. Default 50.
#' @param correlation Logical. Treat the matrices as correlation-scale (unit
#'   diagonal enforced, nearPD run with corr = TRUE). Default TRUE.
#' @return A list with the repaired \code{S10} and \code{S0} (in that order),
#'   plus diagnostics: \code{repaired}, \code{alternating_converged} (only when
#'   repaired), \code{iterations}, and \code{min_eigenvalue} of the final
#'   joint matrix.
#' @keywords internal

check_pd <- function(S10,
                     S0,
                     d,
                     eig.tol = 1e-6,
                     pd.tol = 1e-8,
                     conv.tol = 1e-8,
                     maxit = 50L,
                     correlation = TRUE) {

  sym <- function(A) {
    (A + t(A)) / 2
  }

  build <- function(S0, S10) {
    rbind(
      cbind(S0, t(S10)),
      cbind(S10, S0)
    )
  }

  min_eig <- function(M) {
    min(eigen(
      sym(M),
      symmetric = TRUE,
      only.values = TRUE
    )$values)
  }

  is_pd <- function(M) {
    min_eig(M) > pd.tol
  }

  project_structure <- function(J) {
    i1 <- seq_len(d)
    i2 <- d + seq_len(d)

    J11 <- J[i1, i1, drop = FALSE]
    J12 <- J[i1, i2, drop = FALSE]
    J21 <- J[i2, i1, drop = FALSE]
    J22 <- J[i2, i2, drop = FALSE]

    S0_new <- sym((J11 + J22) / 2)
    S10_new <- (J21 + t(J12)) / 2

    if (correlation) {
      diag(S0_new) <- 1
    }

    list(
      S0 = S0_new,
      S10 = S10_new,
      J = build(S0_new, S10_new)
    )
  }

  if (!is.matrix(S0) || !all(dim(S0) == c(d, d))) {
    stop("S0 must be a d x d matrix.")
  }

  if (!is.matrix(S10) || !all(dim(S10) == c(d, d))) {
    stop("S10 must be a d x d matrix.")
  }

  if (any(!is.finite(S0)) || any(!is.finite(S10))) {
    stop("S0 and S10 must contain only finite values.")
  }

  S0 <- sym(S0)

  if (correlation) {
    diag(S0) <- 1
  }

  J <- build(S0, S10)

  # Already jointly PD: return unchanged apart from the symmetrization and
  # (when correlation = TRUE) the unit diagonal enforced above.
  if (is_pd(J)) {
    return(list(
      S10 = S10,
      S0 = S0,
      repaired = FALSE,
      iterations = 0L,
      min_eigenvalue = min_eig(J)
    ))
  }

  converged <- FALSE
  iterations <- 0L

  # Alternation between:
  # 1. positive-definite matrix with the desired diagonal;
  # 2. stationary block-Toeplitz structure.
  for (iter in seq_len(maxit)) {

    J_old <- J

    Jp <- as.matrix(
      Matrix::nearPD(
        J,
        corr = correlation,
        keepDiag = !correlation,
        eig.tol = eig.tol,
        posd.tol = pd.tol,
        doDykstra = TRUE,
        do2eigen = TRUE
      )$mat
    )

    projected <- project_structure(Jp)

    S0 <- projected$S0
    S10 <- projected$S10
    J <- projected$J

    iterations <- iter

    if (is_pd(J)) {
      converged <- TRUE
      break
    }

    relative_change <-
      norm(J - J_old, type = "F") /
      max(1, norm(J_old, type = "F"))

    if (relative_change < conv.tol) {
      break
    }
  }

  # Fallback: first ensure S0 itself is PD.
  if (!is_pd(S0)) {
    S0 <- sym(as.matrix(
      Matrix::nearPD(
        S0,
        corr = correlation,
        keepDiag = !correlation,
        eig.tol = eig.tol,
        posd.tol = pd.tol,
        doDykstra = TRUE,
        do2eigen = TRUE
      )$mat
    ))

    if (correlation) {
      diag(S0) <- 1
    }
  }

  # Shrink lag covariance toward zero if still needed.
  if (!is_pd(build(S0, S10))) {
    shrink <- 1

    repeat {
      J_candidate <- build(S0, shrink * S10)

      if (is_pd(J_candidate)) {
        S10 <- shrink * S10
        J <- J_candidate
        break
      }

      shrink <- shrink * 0.95

      if (shrink <= 1e-10) {
        S10 <- matrix(0, d, d)
        J <- build(S0, S10)
        break
      }
    }
  }

  # Final certification of exactly what downstream code will use.
  final_min_eigenvalue <- min_eig(J)

  if (!is_pd(J)) {
    stop(
      sprintf(
        paste0(
          "Failed to construct a positive-definite stationary ",
          "joint matrix; minimum eigenvalue = %.3e."
        ),
        final_min_eigenvalue
      )
    )
  }

  list(
    S10 = S10,
    S0 = S0,
    repaired = TRUE,
    alternating_converged = converged,
    iterations = iterations,
    min_eigenvalue = final_min_eigenvalue
  )
}
