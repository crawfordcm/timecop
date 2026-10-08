#' Compute Hermite coefficients
#'
#' @param param Numeric. Estimate of marginal parameter
#' @param k Numeric. The value at which Hermite coefficient infinite sums terminate. Default is 100
#' @param family Character. Name of marginal distribution
#' @param ordinal_levels Numeric. A vector of levels for Ordinal family and NULL otherwise
#' @return A numeric vector of k Hermite coefficients.
#' @keywords internal

hermite_coefs <- function(param, k, family, ordinal_levels) {

  if (family == "Bernoulli") {

    prob <- param
    q <- qnorm(1-prob)

    hk <- lapply(seq_len(k), function(i) {
      her <- as.function(Polys[[i]])
      her(q)
    })

    g <- unlist(lapply(seq_len(k), function(i) {
      exp(-q^2/2) * hk[[i]] / (sqrt(2*pi)*factorial(i))
    }))

  } else if (family == "Poisson") {

    lambda <- param

    g <- unlist(lapply(seq_len(k), function(i) {do.call("sum", lapply(0:50, function(j) {

      # get Q
      c <- ppois(j, lambda)
      q <- qnorm(c)

      # Hermite polynomials/coefficients
      if (c == 1 | c == 0) { # make sure Her isn't Inf
        coef <- 0 # see jia 2021
      } else {
        her <- as.function(Polys[[i]])
        hk <- her(q)
        coef <- exp(-q^2/2) * hk / (sqrt(2*pi)*factorial(i))
      }
      return(coef)

    }))}))

  } else if (family == "Ordinal"){

    probs <- param

    if (length(probs) != length(ordinal_levels)) {
      stop(
        paste0(
          "Hermite coefficient calculation error: the number of probabilities ",
          "does not match the number of ordinal levels."
        ),
        call. = FALSE
      )
    }

    n_levels <- length(probs)

    # need number of levels - 1 here
    q <- qnorm(cumsum(probs)[seq_len(n_levels - 1L)])

    g <- unlist(lapply(seq_len(k), function(i) {

      her <- as.function(Polys[[i]])

      coef <- sum(sapply(q, function(q0) {

        if (!is.finite(q0)) {
          return(0)
        }

        hk <- her(q0)

        exp(-q0^2/2) * hk/(sqrt(2*pi)*factorial(i))
      }))

      return(coef)

    }))


  } else if (family == "Gaussian") {
    g <- numeric(k)
    g[1] <- 1
  }

  return(g)

}
