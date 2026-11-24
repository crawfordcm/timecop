#' Compute penalty weights for the precision matrix
#'
#' Returns a d x d symmetric matrix of entry-specific penalty weights for use
#' in the graphical lasso precision update. The diagonal is always zero since
#' it is never penalised.
#'
#' For \code{penalty = "adaptive"}, the weight for off-diagonal entry (i, j) is
#'   W(i,j) = 1 / (|Omega_current(i,j)| + epsilon)^gamma
#'
#' For \code{penalty = "scad"}, the weight is the SCAD derivative divided by
#' lambda (LLA linearisation):
#'   W(i,j) = 1                                           if |w| <= lambda
#'   W(i,j) = (a*lambda - |w|) / ((a-1) * lambda)        if lambda < |w| <= a*lambda
#'   W(i,j) = 0                                           if |w| > a*lambda
#' where a = scad_a (default 3.7, Fan & Li 2001).
#'
#' @param Omega_current Matrix. d x d current estimate of the precision matrix.
#' @param lambda Numeric. Current penalty parameter. Required for
#'   \code{penalty = "scad"}; ignored for \code{penalty = "adaptive"}.
#' @param penalty Character. Either \code{"adaptive"} or \code{"scad"}.
#'   Default \code{"adaptive"}.
#' @param gamma Numeric. Exponent for adaptive glasso weights. Only used when
#'   \code{penalty = "adaptive"}. Default 1.
#' @param scad_a Numeric. Shape parameter for the SCAD penalty. Default 3.7.
#' @return A d x d symmetric matrix of non-negative weights with zero diagonal.
#' @keywords internal

gvar_weights_Omega <- function(Omega_current, lambda = NULL,
                               penalty = "adaptive", gamma = 1, scad_a = 3.7) {

  if (penalty == "adaptive") {
    eps <- 1e-3 * max(abs(Omega_current), 1e-8)
    W   <- 1 / (abs(Omega_current) + eps)^gamma
    diag(W) <- 0
    return(W)
  }

  if (penalty == "scad") {
    if (is.null(lambda) || lambda < 0) {
      stop("'lambda' must be a non-negative number for SCAD weights", call. = FALSE)
    }
    if (lambda == 0) {
      W <- matrix(1, nrow(Omega_current), ncol(Omega_current))
      diag(W) <- 0
      return(W)
    }
    t <- abs(Omega_current)
    W <- ifelse(t <= lambda,
                1,
                ifelse(t <= scad_a * lambda,
                       (scad_a * lambda - t) / ((scad_a - 1) * lambda),
                       0))
    diag(W) <- 0
    return(W)
  }

  stop("'penalty' must be 'adaptive' or 'scad'", call. = FALSE)
}
