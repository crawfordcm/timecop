#' Compute numerical derivatives
#'
#' @param data Matrix. A d by n multivariate time series matrix
#' @param cov_x_hat Array. An array of observed covariances
#' @param d Numeric. The number of variables
#' @param p Numeric. The VAR order
#' @param n Numeric. Time series length
#' @param family List. A list of marginal distributions
#' @param corr Logical. Correlations or covariances
#' @return A Jacobian matrix of numerical derivatives.
#' @importFrom numDeriv jacobian
#' @keywords internal

numderiv <- function(data, cov_x_hat, d, p, n, family, corr) {

  x1 <- numeric()
  for (i in seq_along(family)) {
    if (family[[i]] %in% c("Bernoulli", "Poisson")) {
      x1 <- c(x1,
              mean(data[i,], na.rm = TRUE))
    } else if (family[[i]] == "Gaussian") {
      x1 <- c(x1,
              mean(data[i,], na.rm = TRUE),
              var(data[i,],  na.rm = TRUE))
    }
  }

  marg_num <- length(x1)

  x2 <- c(vec(cov_x_hat[,,2]), vec(cov_x_hat[,,1]), vec(cov_x_hat[,,3]))
  x  <- c(x1, x2)

  func <- function(x, d, marg_num, p, n, data, family, corr) {

    param_hat <- vector("list", d)
    j <- 1
    for (i in seq_along(family)) {
      if (family[[i]] %in% c("Bernoulli", "Poisson")) {
        param_hat[[i]] <- x[j]
        j <- j + 1
      } else if (family[[i]] == "Gaussian") {
        param_hat[[i]] <- c(x[j], x[j+1])
        j <- j + 2
      }
    }

    cov_x_hat <- array(0, dim=c(d,d,3))
    cov_x_hat[,,2] <- revVec(x[(marg_num + 1):(marg_num + d^2)])
    cov_x_hat[,,1] <- revVec(x[(marg_num + d^2 + 1):(marg_num + 2*d^2)])
    cov_x_hat[,,3] <- revVec(x[(marg_num + 2*d^2 + 1):(marg_num + 3*d^2)])

    ell_ij_hat <- latent_var_link_numderiv(d, n, param_hat, family, corr)
    cov_z_hat <- latent_var_invlink(cov_x_hat, d, p, ell_ij_hat)
    res <- c(vec(cov_z_hat[,,1]), vec(cov_z_hat[,,2]), vec(cov_z_hat[,,3]))

    return(res)
  }

  jacob <- jacobian(func, x, d = d, marg_num = marg_num, p = p, n = n, data = data, family = family, corr = corr)

  return(jacob)

}
