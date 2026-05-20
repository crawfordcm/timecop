#' Latent link function
#'
#' @param data Matrix. A d by n multivariate time series matrix
#' @param d Numeric. The number of variables
#' @param n Numeric. Time series length
#' @param k Numeric. The value at which Hermite coefficient infinite sums terminate. Default is 100
#' @param family List. A list of length d with the names of each distribution
#' @param corr Logical. Correlations or covariances
#' @return A k x d x d array of link function coefficients.
#' @keywords internal

latent_var_link <- function(data, d, n, k, family, corr) {
  param_hat <- estimate_marginal_params(data, family, d)
  .link_coef_array(d, k, param_hat, family, corr)
}

#' Latent link function for numderiv
#'
#' @param d Numeric. The number of variables
#' @param n Numeric. Time series length
#' @param param_hat List. List of estimated marginal parameters
#' @param family List. List of marginal distributions
#' @param corr Logical. Correlations or covariances
#' @return A k x d x d array of link function coefficients.
#' @keywords internal

latent_var_link_numderiv <- function(d, n, param_hat, family, corr) {
  .link_coef_array(d, k = 100, param_hat, family, corr)
}

.link_coef_array <- function(d, k, param_hat, family, corr) {
  ell_ij_hat <- array(NA, dim = c(k, d, d))
  for (i in seq_len(d)) {
    param_list  <- list(param_hat[[i]])
    family_list <- list(family[[i]])
    for (j in seq_len(d)) {
      if (j >= i) {
        param_list[[2]]  <- param_hat[[j]]
        family_list[[2]] <- family[[j]]
        ell_ij_hat[, i, j] <- link_coefs(param_list, k, family_list, corr)
      } else {
        ell_ij_hat[, i, j] <- ell_ij_hat[, j, i]
      }
    }
  }
  ell_ij_hat
}
