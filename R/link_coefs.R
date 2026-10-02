#' Compute Hermite coefficients and link function constants
#'
#' @param param_list List. A two-entry list of marginal parameters
#' @param k Numeric. The value at which Hermite coefficient infinite sums terminate. Default is 100
#' @param family_list List. A two-entry list of marginal distributions
#' @param ordinal_levels List. A two-entry list containing the ordinal levels
#'   for each variable and `NULL` for non-Ordinal variables.
#' @param corr Logical. Correlations or covariances
#' @return A numeric vector of k link function coefficients.
#' @keywords internal

link_coefs <- function(param_list, k, family_list, ordinal_levels, corr){

  # get Hermite coefficients
  g_i <- hermite_coefs(param_list[[1]], k, family_list[[1]], ordinal_levels[[1]])
  g_j <- hermite_coefs(param_list[[2]], k, family_list[[2]], ordinal_levels[[2]])

  # Polynomial of link function
  if (corr) {

    # get SDs
    sd_marg <- vapply(seq_along(family_list), function(i) {
      marginal_sd(
        param = param_list[[i]],
        family = family_list[[i]],
        ordinal_levels = ordinal_levels[[i]]
      )}, numeric(1))

    l_ij <- (g_i * g_j * factorial(seq_len(k)) / prod(sd_marg))

  } else {

    l_ij <- g_i * g_j * factorial(seq_len(k))

  }

  return(l_ij)
}

