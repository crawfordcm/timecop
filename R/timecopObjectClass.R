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
#' @export


setClass(
  Class = "timecop",
  slots = list(
    data = "matrix",
    cov_x_hat = "array",
    cov_z_hat = "array",
    gamma_hat = "matrix", # lag 1 matrix (if p = 1)
    Gamma_hat = "matrix", # lag 0 matrix (if p = 1)
    d = "numeric",
    p = "numeric",
    n = "numeric",
    marg_num = "numeric",
    family = "list",
    pd_approx = "logical"
  ), validity = check_timecop
)

#' Construct an object of class timecop
#'
#' @param data Matrix. An n (time points) by d (variables) multivariate time series matrix
#' @param family List. A list of length d with the names of each distribution
#' @param p Numeric. The VAR order. Default is 1
#' @export
timecop <- function(data = NULL,
                    family = NULL,
                    p = 1) {

  if (p != 1) {
    stop("Only lag 1 is currently supported", call. = FALSE)
  }

  # setup data
  data <- setup_data(data)

  d <- nrow(data)   # number of variables
  n <- ncol(data)   # time series length

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
  cov_x_hat  <- observed_var_cov(data, d, p, n)
  ell_ij_hat <- latent_var_link(data, d, n, k, family)
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

  # check positive definiteness of Gamma_hat
  # this is inelegant: could be better
  pd_approx <- !is.positive.definite(Gamma_hat)

  if (pd_approx) {
    eig.tol  <- 1e-6
    posd.tol <- 1e-6
    pd_mat <- nearPD(Gamma_hat, corr = TRUE, keepDiag = TRUE, eig.tol = eig.tol, posd.tol = posd.tol)
    Gamma_hat <- as.matrix(pd_mat$mat)

    while (any(abs(solve(Gamma_hat)) > 5)) {
      eig.tol <- eig.tol * 10
      posd.tol <- posd.tol * 10
      pd_mat <- nearPD(Gamma_hat, corr = TRUE, keepDiag = TRUE, eig.tol = eig.tol, posd.tol = posd.tol)
      Gamma_hat <- as.matrix(pd_mat$mat)
    }
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
    pd_approx = pd_approx
  )

  return(obj)

}

#' Fit function for timecop
#'
#' @usage fit_timecop(object)
#' @param object timecop object built using \code{timecop}

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
    object@marg_num
  )

  return(list(est, se))

})





