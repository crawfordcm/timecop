#' Compute Hermite coefficients
#'
#' @param param Numeric. Estimate of marginal parameter
#' @param k Numeric. The value at which Hermite coefficient infinite sums terminate. Default is 100
#' @param family Character. Name of marginal distribution
#' @keywords internal

hermite_coefs <- function(param, k, family){
  prob <- param

  if (family == "Bernoulli"){
    q <- qnorm(1-prob)

    hk <- lapply(seq_len(k), function(i) {
      her <- as.function(Polys[[i]])
      her(q)
    })

    g <- unlist(lapply(seq_len(k), function(i) {
      exp(-q^2/2) * hk[[i]] / (sqrt(2*pi)*factorial(i))
    }))

  } else if (family == "Poisson"){
    lambda <- param

    g <- unlist(lapply(seq_len(k), function(i) {do.call("sum", lapply(0:50, function(j) {
      # get Q
      c <- ppois(j, lambda)
      q <- qnorm(c)

      # Hermite polynomials/coefficients
      if (c == 1 | c == 0) { # make sure Her isn't Inf
        coef <- 0
      } else {
        her <- as.function(Polys[[i]])
        hk <- her(q)
        coef <- exp(-q^2/2) * hk / (sqrt(2*pi)*factorial(i))
      }
      return(coef)

    }))}))

  } else if (family == "Gaussian"){
    g <- numeric(k)
    g[1] <- 1
  }

  return(g)

}
