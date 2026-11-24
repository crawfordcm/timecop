#' Compute penalty weights for the temporal coefficient matrix
#'
#' Returns a d x d matrix of entry-specific penalty weights for use in the
#' FISTA temporal update.
#'
#' For \code{penalty = "adaptive"}, the weight for entry (i, j) is
#'   W(i,j) = 1 / (|A_current(i,j)| + epsilon)^gamma
#'
#' For \code{penalty = "scad"}, the weight is the SCAD derivative divided by
#' lambda (LLA linearisation):
#'   W(i,j) = 1                                       if |a| <= lambda
#'   W(i,j) = (a*lambda - |a|) / ((a-1) * lambda)    if lambda < |a| <= a*lambda
#'   W(i,j) = 0                                       if |a| > a*lambda
#' where a = scad_a (default 3.7, Fan & Li 2001).
#'
#' @param A_current Matrix. d x d current estimate of the temporal coefficient
#'   matrix.
#' @param lambda Numeric. Current penalty parameter. Required for
#'   \code{penalty = "scad"}; ignored for \code{penalty = "adaptive"}.
#' @param penalty Character. Either \code{"adaptive"} or \code{"scad"}.
#'   Default \code{"adaptive"}.
#' @param gamma Numeric. Exponent for adaptive lasso weights. Only used when
#'   \code{penalty = "adaptive"}. Default 1.
#' @param scad_a Numeric. Shape parameter for the SCAD penalty. Default 3.7.
#' @return A d x d matrix of non-negative weights.
#' @keywords internal

gvar_weights_A <- function(A_current, lambda = NULL,
                           penalty = "adaptive", gamma = 1, scad_a = 3.7) {

  if (penalty == "adaptive") {
    eps <- 1e-3 * max(abs(A_current), 1e-8)
    return(1 / (abs(A_current) + eps)^gamma)
  }

  if (penalty == "scad") {
    if (is.null(lambda) || lambda < 0) {
      stop("'lambda' must be a non-negative number for SCAD weights", call. = FALSE)
    }
    if (lambda == 0) return(matrix(1, nrow(A_current), ncol(A_current)))
    t <- abs(A_current)
    W <- ifelse(t <= lambda,
                1,
                ifelse(t <= scad_a * lambda,
                       (scad_a * lambda - t) / ((scad_a - 1) * lambda),
                       0))
    return(W)
  }

  stop("'penalty' must be 'adaptive' or 'scad'", call. = FALSE)
}
