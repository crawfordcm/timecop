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


setClass(
  Class = "timecop",
  slots = list(
    cov_z_hat = "array",
    #gamma_hat = "matrix", # lag 1 matrix (if p = 1)
    #Gamma_hat = "matrix", # lag 0 matrix (if p = 1)
    cov_x_hat = "array",
    d = "numeric",
    p = "numeric",
    n = "numeric",
    family = "list"
  ), validity = check_timecop
)

#' Construct an object of class timecop
#'
#' @param data Matrix. An n (time points) by d (variables) multivariate time series matrix
#' @param p Numeric. The VAR order. Defualt is 1
#' @param family List. A list of length d with the names of each distribution
#' @param grid Numeric. A vector of grid values used for interpolation. Default is NULL
#' @param k Numeric. The value at which Hermite coefficient infinite sums terminate. Default is 100
#' @export
timecop <- function(data = NULL,
                    p = 1,
                    family = NULL,
                    grid = NULL,
                    k = 100) {

  if (p != 1) {
    stop("Only lag 1 is currently supported", call. = FALSE)
  }

  # setup data
  data <- setup_data(data)

  d <- nrow(data)   # number of variables
  n <- ncol(data)   # time series length

  # check if marginals are correctly specified
  check_marginals(family, d)

  # setup grid for interpolation
  if (is.null(grid)) {
    grid_u <- c(seq(-1,-0.71,0.01), seq(-0.7,0,0.2), seq(0.1,0.7,0.2), seq(0.71,1,0.01))
  } else {
    grid_u <- grid
  }

  # compute observed and latent covariances
  cov_x_hat  <- observed_var_cov(data, d, p, n)
  ell_ij_hat <- latent_var_link(data, d, n, k, family)
  cov_z_hat  <- latent_var_invlink(cov_x_hat, d, p, ell_ij_hat, grid_u)

  # construct implied latent covariances
  gamma_hat <- array(NA,dim=c(p*d,d))
  for (h in 1:p) {
    gamma_hat[((h-1)*d+1):(h*d),] <- cov_z_hat[,,(p+1-h)]
  }

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

  obj <- new("timecop",
             cov_z_hat = cov_z_hat,
             cov_x_hat = cov_x_hat,
             d = d,
             p = p,
             n = n,
             family = family)

  return(obj)

}

#' Fit function for timecop
#'
#'






