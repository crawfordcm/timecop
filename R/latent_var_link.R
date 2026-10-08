#' Latent link function
#'
#' @param data Matrix. A d by n multivariate time series matrix
#' @param d Numeric. The number of variables
#' @param n Numeric. Time series length
#' @param k Numeric. The value at which Hermite coefficient infinite sums terminate. Default is 100
#' @param family List. A list of length d with the names of each distribution
#' @param ordinal_levels List. A list of length d containing the numeric levels
#'   for each Ordinal variable and `NULL` for each non-Ordinal variable.
#' @param corr Logical. Correlations or covariances
#' @return A k x d x d array of link function coefficients.
#' @keywords internal

latent_var_link <- function(data, d, n, k, family, ordinal_levels, corr) {
  param_hat <- estimate_marginal_params(data, family, ordinal_levels, d)
  .link_coef_array(d, k, param_hat, family, ordinal_levels, corr)
}

#' Latent link function for numderiv
#'
#' @param d Numeric. The number of variables
#' @param n Numeric. Time series length
#' @param param_hat List. List of estimated marginal parameters
#' @param family List. List of marginal distributions
#' @param ordinal_levels List. A list of length d containing the numeric levels
#'   for each Ordinal variable and `NULL` for each non-Ordinal variable.
#' @param corr Logical. Correlations or covariances
#' @return A k x d x d array of link function coefficients.
#' @keywords internal

latent_var_link_numderiv <- function(d, n, param_hat, family, ordinal_levels, corr) {
  .link_coef_array(d, k = 100, param_hat, family, ordinal_levels, corr)
}

.link_coef_array <- function(d, k, param_hat, family, ordinal_levels, corr) {

  ell_ij_hat <- array(NA, dim = c(k, d, d))

  for (i in seq_len(d)) {
    for (j in seq_len(d)) {
      if (j >= i) {
        param_list <- param_hat[c(i, j)]
        family_list <- family[c(i, j)]
        levels_list <- ordinal_levels[c(i, j)]

        ell_ij_hat[, i, j] <- link_coefs(
          param_list = param_list,
          k = k,
          family_list = family_list,
          ordinal_levels = levels_list,
          corr = corr
        )
      } else {
        ell_ij_hat[, i, j] <- ell_ij_hat[, j, i]
      }
    }
  }

  ell_ij_hat
}
