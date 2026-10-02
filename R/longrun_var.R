#' Compute long-run variance
#'
#' @param data Matrix. A d by n multivariate time series matrix
#' @param d Numeric. The number of variables
#' @param n Numeric. Time series length
#' @param family List. A list of marginal distributions
#' @param ordinal_levels List. A list of length d containing the numeric levels
#'   for each Ordinal variable and `NULL` for each non-Ordinal variable.
#' @return A long-run variance-covariance matrix.
#' @importFrom cointReg getLongRunVar
#' @keywords internal

longrun_var <- function(data, d, n, family, ordinal_levels) {

  # center data
  dat_c <- data - rowMeans(data)

  F_list <- list()

  # estimate marginal probabilities for Ordinal case
  if ("Ordinal" %in% family) {
    param_hat <- estimate_marginal_params(data, family, ordinal_levels, d)
  }

  for(i in 2:(n-1)){

    xt <- matrix(dat_c[,i], d, 1)

    a1 <- xt %*% t(xt)
    a2 <- matrix(dat_c[,i-1], d, 1) %*% t(xt)
    a  <- vec(rbind(a1,a2))

    b1 <- matrix(dat_c[,i+1], d, 1) %*% t(xt)
    b2 <- xt %*% t(xt)
    b  <- vec(rbind(b1, b2))

    # added Ordinal-Ordinal case
    rt <- numeric(0)

    for (j in seq_along(family)) {

      if (family[[j]] == "Ordinal") {

        levels_j <- ordinal_levels[[j]]
        probs_j <- param_hat[[j]]
        x_raw <- data[j, i - 1L]

        # Only K - 1 probabilities are free parameters.
        ordinal_scores <- vapply(
          seq_len(length(levels_j) - 1L),
          function(h) {
            as.numeric(x_raw == levels_j[h]) - probs_j[h]
          },
          numeric(1)
        )

        rt <- c(rt, ordinal_scores)

      } else if (family[[j]] == "Gaussian") {

        # TODO: Confirm the correct time index for the marginal estimating-function
        # components. The previous Gaussian implementation used dat_c[j, i] through
        # xt[j], whereas the current implementation uses dat_c[j, i - 1L], consistent
        # with the Ordinal, Bernoulli, and Poisson components. Verify this choice
        # against the long-run variance derivation and Gaussian benchmark tests.

        x_centered <- dat_c[j, i - 1L]

        # Gaussian has two marginal parameters.
        rt <- c(
          rt,
          x_centered,
          x_centered^2
        )

      } else {

        # Bernoulli and Poisson each have one marginal parameter.
        rt <- c(rt, dat_c[j, i - 1L])
      }
    }

    rt <- matrix(rt, ncol = 1L)
    F_list[[i - 1L]] <- rbind(rt, a, b)
  }

  F_mat <- do.call(cbind, F_list)

  longrun <- getLongRunVar(
    t(F_mat),
    bandwidth = "nw",
    kernel = "ba",
    demeaning = TRUE,
    check = FALSE
  )

  Sigma <- longrun$Omega

  return(Sigma)
}

