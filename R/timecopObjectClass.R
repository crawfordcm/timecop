# check timecop object
check_timecop <- function(object) {

  errors <- character()

  if (any(is.na(object@data))) {
    msg    <- c("timecop error: remove NA values before proceeding")
    errors <- c(errors, msg)
  }

  if (length(family) != d) {
    msg    <- c("timecop error: number of families must match number of variables")
    errors <- c(errors, msg)
  }

  if (p != 1) {
    msg    <- c("timecrop error: only lag 1 is currently supported")
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
    gamma_hat = "matrix", # lag 1 matrix (if p = 1)
    Gamma_hat = "matrix", # lag 0 matrix (if p = 1)
    Cov_X_hat = "matrx",
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

  cov_x_hat  <- observed_var_cov(data, d, p, n)
  ell_ij_hat <- latent_var_link(data, d, n, k, family)
  cov_z_hat  <- latent_var_invlink(cov_x_hat, d, p, ell_ij_hat, grid_u)

}


