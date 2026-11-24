#' SCAD penalty value
#'
#' Elementwise SCAD penalty p_lambda(theta) of Fan & Li (2001), evaluated at
#' theta = |x|. Returns a value/matrix of the same shape as \code{theta}.
#'   p(t) = lambda * t                                 if t <= lambda
#'   p(t) = (2*a*lambda*t - t^2 - lambda^2)/(2*(a-1)) if lambda < t <= a*lambda
#'   p(t) = lambda^2 * (a + 1) / 2                     if t > a*lambda
#'
#' @param theta Numeric or matrix. Non-negative magnitudes.
#' @param lambda Numeric. Penalty parameter.
#' @param a Numeric. SCAD shape parameter (default 3.7).
#' @return A value/matrix of penalty contributions, same shape as \code{theta}.
#' @keywords internal
gvar_scad_penalty <- function(theta, lambda, a = 3.7) {
  if (lambda == 0) return(theta * 0)
  ifelse(theta <= lambda,
          lambda * theta,
          ifelse(theta <= a * lambda,
                 (2 * a * lambda * theta - theta^2 - lambda^2) / (2 * (a - 1)),
                 lambda^2 * (a + 1) / 2))
}

#' MRCE objective for graphical VAR
#'
#' Computes the MRCE objective
#'   tr(S_eps(A) * Omega) - log|Omega| + P_Omega(Omega) + 2 * P_A(A)
#' where the penalty terms depend on \code{penalty}:
#' \describe{
#'   \item{lasso}{P_Omega = lambda_Omega * sum over i,j off-diagonal of |omega_ij|,
#'     P_A = lambda_A * sum |a_ij|.}
#'   \item{adaptive}{weighted versions using \code{W_Omega} and \code{W_A};
#'     each weight matrix multiplies its entries elementwise. \code{W_Omega}
#'     is assumed to have a zero diagonal.}
#'   \item{scad}{P_Omega = sum over i,j off-diagonal of p_scad(|omega_ij|),
#'     P_A = sum p_scad(|a_ij|), with the non-convex SCAD penalty.}
#' }
#'
#' @param A Matrix. d x d temporal coefficient matrix.
#' @param Omega Matrix. d x d innovation precision matrix.
#' @param S0_plus Matrix. Lag-0 current-side covariance.
#' @param S0_minus Matrix. Lag-0 past-side covariance.
#' @param S01 Matrix. Lag-1 past-current covariance.
#' @param S10 Matrix. Lag-1 current-past covariance.
#' @param lambda_A Numeric. Temporal penalty parameter.
#' @param lambda_Omega Numeric. Precision matrix penalty parameter.
#' @param penalty Character. One of \code{"lasso"}, \code{"adaptive"},
#'   \code{"scad"}. Default \code{"lasso"}.
#' @param W_A Matrix or NULL. Adaptive weights for A. Only used when
#'   \code{penalty = "adaptive"}.
#' @param W_Omega Matrix or NULL. Adaptive weights for Omega (zero diagonal).
#'   Only used when \code{penalty = "adaptive"}.
#' @param scad_a Numeric. SCAD shape parameter. Only used when
#'   \code{penalty = "scad"}. Default 3.7.
#' @return Numeric scalar. Value of the MRCE objective.
#' @keywords internal

gvar_objective <- function(A, Omega, S0_plus, S0_minus, S01, S10,
                           lambda_A, lambda_Omega,
                           penalty = "lasso", W_A = NULL, W_Omega = NULL,
                           scad_a = 3.7) {

  S_eps   <- gvar_resid_cov(A, S0_plus, S0_minus, S01, S10)
  log_det <- as.numeric(determinant(Omega, logarithm = TRUE)$modulus)

  # tr(S_eps %*% Omega) — valid since Omega is symmetric
  tr_term <- sum(S_eps * Omega)

  if (penalty == "lasso") {
    off_sum   <- sum(abs(Omega)) - sum(abs(diag(Omega)))
    pen_Omega <- lambda_Omega * off_sum
    pen_A     <- 2 * lambda_A * sum(abs(A))

  } else if (penalty == "adaptive") {
    # W_Omega carries a zero diagonal, so the diagonal is excluded automatically
    pen_Omega <- lambda_Omega * sum(W_Omega * abs(Omega))
    pen_A     <- 2 * lambda_A * sum(W_A * abs(A))

  } else if (penalty == "scad") {
    Omega_off <- Omega
    diag(Omega_off) <- 0   # p_scad(0) = 0, so diagonal contributes nothing
    pen_Omega <- sum(gvar_scad_penalty(abs(Omega_off), lambda_Omega, scad_a))
    pen_A     <- 2 * sum(gvar_scad_penalty(abs(A), lambda_A, scad_a))

  } else {
    stop("'penalty' must be 'lasso', 'adaptive', or 'scad'", call. = FALSE)
  }

  tr_term - log_det + pen_Omega + pen_A
}

# Smooth part only (used inside FISTA for convergence checks)
gvar_objective_smooth <- function(A, Omega, S0_plus, S0_minus, S01, S10) {
  S_eps <- gvar_resid_cov(A, S0_plus, S0_minus, S01, S10)
  0.5 * sum(S_eps * Omega)
}
