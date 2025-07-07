# check timecop object
check_timecop <- function(object){

  errors <- character()

  if (any(is.na(object@data))){
    msg    <- c("timecop error: remove NA values before proceeding")
    errors <- c(errors, msg)
  }

  if (length(family) != d){
    msg    <- c("timecop error: number of families must match number of variables")
    errors <- c(errors, msg)
  }

  if (p != 1){
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
#' @param data Matrix. An n by d multivariate time series matrix
#' @param p Numeric. The VAR order. Defualt is 1
#' @param d Numeric. The number of variables
#' @param n Numeric. Time series length
#' @param family List. A list of length d with the names of each distribution
#' @param grid Numeric. A vector of grid values used for interpolation. Default is NULL
#' @param k Numeric. The value at which Hermite coefficient infinite sums terminate. Default is 100
#' @export
timecop <- function(data = NULL,
                    p = 1,
                    d = NULL,
                    n = NULL,
                    family = NULL,
                    grid = NULL,
                    k = 100){

  if (p != 1){
    stop("timecrop error: only lag 1 is currently supported")
  }

  if (length(family) != d){
    stop("timecop error: number of families must match number of variables")
  }

  if (any)

  data <- t(data)

  if (is.null(grid)){
    grid_u <- c(seq(-1,-0.71,0.01),seq(-0.7,0,0.2),seq(0.1,0.7,0.2),seq(0.71,1,0.01))
  } else{
    grid_u <- grid
  }

  ell_ij_hat <- latent_var_link(d = d, n, data, )

}


check_distributions <- function(dist_list) {
  # Accepted distributions — stored as a list
  accepted <- list("Gaussian", "Poisson", "Bernoulli")

  ## ---- case-insensitive membership test -------------------------------
  dist_lc     <- tolower(unlist(dist_list))   # flatten user list
  accepted_lc <- tolower(unlist(accepted))    # flatten accepted list

  bad_idx <- which(!dist_lc %in% accepted_lc) # indices of problematic entries

  if (length(bad_idx)) {
    bad <- unique(unlist(dist_list[bad_idx]))

    stop(
      sprintf(
        "Distribution%s not supported: %s.\nAllowed: %s",
        if (length(bad) > 1) "s" else "",          # plural “s” if needed
        paste(bad, collapse = ", "),
        paste(unlist(accepted), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)    # silent success
}

## ── examples ────────────────────────────────────────────────
check_distributions(list("Gaussian", "Gaussian", "Poisson"))
# (returns silently)

check_distributions(list("Gaussian", "Gamma"))
# Error: Distribution not supported: Gamma.
#        Allowed: Gaussian, Poisson, Bernoulli

check_distributions(list("Gamma", "Beta", "Beta"))
# Error: Distributions not supported: Gamma, Beta.
#        Allowed: Gaussian, Poisson, Bernoulli



