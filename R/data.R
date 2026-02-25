#' Example bivariate Bernoulli time series
#'
#' A simulated 200 x 2 matrix of bivariate Bernoulli time series data generated
#' from a latent Gaussian VAR(1) process. Generated using [latent_var_sim()]
#' with `set.seed(42)`, success probabilities of 0.5 for both variables, and
#' a transition matrix with 0.4 on the diagonal and 0.2 on the off-diagonal.
#'
#' @format A 200 x 2 numeric matrix. Each column is a binary (0/1) time series
#'   of length 200.
#' @source Simulated data.
"timecop_example"
