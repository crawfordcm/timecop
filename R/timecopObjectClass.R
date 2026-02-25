# check timecop object
check_timecop <- function(object) {

  errors <- character()

  if (any(is.na(object@data))) {
    msg    <- c("timecop error: remove NA values before proceeding")
    errors <- c(errors, msg)
  }

  if (length(errors) == 0) TRUE else errors

}

#' timecop object class
#'
#' An object to be used with fit_timecop
#'
#' @slot data Matrix. A d (variables) by n (time points) multivariate time series matrix
#' @slot cov_x_hat Array. An array of observed covariances
#' @slot cov_z_hat Array. An array of latent covariances
#' @slot gamma_hat Matrix. A matrix of stacked covariances (lag 1 to p)
#' @slot Gamma_hat Matrix. A Toeplitz matrix of covariances
#' @slot d Numeric. The number of variables
#' @slot p Numeric. The VAR order. Default is 1
#' @slot n Numeric. Time series length
#' @slot marg_num Numeric. Total number of marginal parameters
#' @slot family List. A list of marginal distributions
#' @slot pd_approx Logical. A logical indicating if adjustments were made to ensure positive definiteness
#' @slot corr Logical. A logical indicating if observed correlations or covariances were used
#' @export


setClass(
  Class = "timecop",
  slots = list(
    data = "matrix",
    cov_x_hat = "array",
    cov_z_hat = "array",
    gamma_hat = "matrix",
    Gamma_hat = "matrix",
    d = "numeric",
    p = "numeric",
    n = "numeric",
    marg_num = "numeric",
    family = "list",
    pd_approx = "logical",
    corr = "logical"
  ), validity = check_timecop
)

#' Construct an object of class timecop
#'
#' Takes an n (time points) by d (variables) data matrix and a list of marginal
#' distribution families, computes observed and latent covariances, and builds
#' the Yule-Walker matrices needed for VAR estimation.
#'
#' @param data Matrix. An n (time points) by d (variables) multivariate time series matrix
#' @param family List. A list of length d with the names of each distribution.
#'   Supported values: `"Bernoulli"`, `"Poisson"`, `"Gaussian"`.
#' @param p Numeric. The VAR order. Default is 1. Only p=1 is currently supported.
#' @param corr Logical. Use correlations instead of covariances. Default is `FALSE`.
#' @param pd_approx Logical. Check if latent covariance matrices are PD and
#'   apply nearest positive-definite approximation if not. Default is `FALSE`.
#' @param eig.tol Numeric. Eigenvalue tolerance for the positive-definite
#'   approximation. Only relevant if `pd_approx = TRUE`. Default is `1e-4`.
#' @return An S4 object of class [timecop-class].
#' @examples
#' \donttest{
#' sim <- latent_var_sim(
#'   d = 2, n = 200, p = 1,
#'   param = list(0.5, 0.5),
#'   phi_lv = matrix(c(0.4, 0.2, 0.2, 0.4), 2, 2),
#'   family = list("Bernoulli", "Bernoulli")
#' )
#' obj <- timecop(data = t(sim$X_t), family = list("Bernoulli", "Bernoulli"))
#' }
#' @export
timecop <- function(data = NULL,
                    family = NULL,
                    p = 1,
                    corr = FALSE,
                    pd_approx = FALSE,
                    eig.tol = 1e-4) {

  if (!is.logical(corr) || length(corr) != 1) {
    stop("'corr' must be a single logical value", call. = FALSE)
  }
  if (!is.logical(pd_approx) || length(pd_approx) != 1) {
    stop("'pd_approx' must be a single logical value", call. = FALSE)
  }
  if (!is.numeric(p) || length(p) != 1 || p != as.integer(p) || p < 1) {
    stop("'p' must be a positive integer", call. = FALSE)
  }
  if (p != 1) {
    stop("Only lag 1 is currently supported", call. = FALSE)
  }
  if (!is.numeric(eig.tol) || length(eig.tol) != 1 || eig.tol <= 0) {
    stop("'eig.tol' must be a positive numeric value", call. = FALSE)
  }

  # setup data
  data <- setup_data(data)

  d <- nrow(data)   # number of variables
  n <- ncol(data)   # time series length

  if (n < 2 * p + 2) {
    stop(sprintf("Time series length n=%d is too short (need at least %d)", n, 2 * p + 2),
         call. = FALSE)
  }

  # check if marginals are correctly specified
  check_marginals(family, d)

  # get total number of marginal parameters
  if ("Gaussian" %in% family) {
    marg_num <- d + length(which(family == "Gaussian"))
  } else {
    marg_num <- d
  }

  k <- 100

  # compute observed and latent covariances
  cov_x_hat  <- observed_var_cov(data, d, p, n, corr)
  ell_ij_hat <- latent_var_link(data, d, n, k, family, corr)
  cov_z_hat  <- latent_var_invlink(cov_x_hat, d, p, ell_ij_hat)

  # construct covariance matrices for Yule-Walker
  # stacked covariances
  gamma_hat <- array(NA,dim=c(p*d,d))
  for (h in 1:p) {
    gamma_hat[((h-1)*d+1):(h*d),] <- cov_z_hat[,,(p+1-h)]
  }

  # Toeplitz covariance matrix
  Gamma_hat <- array(NA,dim=c(p*d,p*d))
  for (i in 1:p) {
    for (j in 1:p) {
      if (i < j) {
        Gamma_hat[((i-1)*d+1):(i*d),((j-1)*d+1):(j*d)] <- cov_z_hat[,,(p+1+abs(i-j))]
      }else if (i == j) {
        Gamma_hat[((i-1)*d+1):(i*d),((j-1)*d+1):(j*d)] <- cov_z_hat[,,(p+1)]
      } else { # i > j
        Gamma_hat[((i-1)*d+1):(i*d),((j-1)*d+1):(j*d)] <- cov_z_hat[,,(p+1-abs(i-j))]
      }
    }
  }

  # ensure that the interpolation didn't compute impossible values
  gamma_hat[gamma_hat > 1] <- 1
  gamma_hat[gamma_hat < -1] <- -1
  Gamma_hat[Gamma_hat > 1] <- 1
  Gamma_hat[Gamma_hat < -1] <- -1

  # check if covariances matrices are PD (only works for p = 1 currently)
  if (pd_approx) {
    cov_mats  <- check_pd(gamma_hat, Gamma_hat, d, eig.tol)
    gamma_hat <- cov_mats[[1]]
    Gamma_hat <- cov_mats[[2]]
  }

  obj <- new(
    "timecop",
    data = data,
    cov_x_hat = cov_x_hat,
    cov_z_hat = cov_z_hat,
    gamma_hat = gamma_hat,
    Gamma_hat = Gamma_hat,
    d = d,
    p = p,
    n = n,
    marg_num = marg_num,
    family = family,
    pd_approx = pd_approx,
    corr = corr
  )

  return(obj)

}

#' Fit a latent Gaussian VAR model
#'
#' Estimates the VAR(1) transition matrix and standard errors from a
#' [timecop-class] object using the Yule-Walker method.
#'
#' @usage fit_timecop(object)
#' @param object A [timecop-class] object built using [timecop()].
#' @return An object of class \code{"timecop_fit"} (a named list) with elements:
#'   \item{estimates}{A d x d matrix of estimated VAR(1) transition coefficients.}
#'   \item{se}{A d x d matrix of standard errors for the estimates.}
#'   \item{obj}{The original \code{\linkS4class{timecop}} object.}
#'
#' Use \code{summary()} on the result for a formatted table with z-values and
#' p-values.
#' @examples
#' \donttest{
#' sim <- latent_var_sim(
#'   d = 2, n = 200, p = 1,
#'   param = list(0.5, 0.5),
#'   phi_lv = matrix(c(0.4, 0.2, 0.2, 0.4), 2, 2),
#'   family = list("Bernoulli", "Bernoulli")
#' )
#' obj <- timecop(data = t(sim$X_t), family = list("Bernoulli", "Bernoulli"))
#' results <- fit_timecop(obj)
#' results$estimates       # coefficient matrix
#' results$se              # standard errors
#' summary(results)        # formatted summary table
#' }
#' @export
#' @aliases fit_timecop,timecop-method

setGeneric(name = "fit_timecop", def = function(object){standardGeneric("fit_timecop")})
setMethod(f = "fit_timecop", signature = "timecop", definition = function(object) {

  est <- fit_var(
    object@gamma_hat,
    object@Gamma_hat,
    object@d
  )

  se <- se_var(
    object@data,
    object@gamma_hat,
    object@Gamma_hat,
    object@cov_x_hat,
    object@d,
    object@p,
    object@n,
    object@family,
    object@marg_num,
    object@corr
  )

  results <- list(
    estimates = est,
    se        = se,
    obj       = object
  )
  class(results) <- c("timecop_fit", "list")

  return(results)

})
