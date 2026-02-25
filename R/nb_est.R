#' Estimate negative binomial parameters
#'
#' Method-of-moments estimator for the negative binomial distribution.
#'
#' TODO: Integrate with main estimation pipeline for negative binomial marginals.
#'
#' @param data Numeric. A vector of count data.
#' @param n Numeric. The length of the data vector.
#' @return A list with components `r` (number of successes) and `p` (probability
#'   of success).
#' @keywords internal
nb_est <- function(data, n) {

  mu <- mean(data)
  s2 <- var(data)
  r  <- mu^2 / (s2 - mu)
  p  <- r / (r + mu)

  return(list(r = r, p = p))
}
