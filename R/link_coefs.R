#' Compute Hermite coefficients and link function constants
#'
#' @param param_list List. A list of marginal parameters
#' @param k Numeric. The value at which Hermite coefficient infinite sums terminate. Default is 100
#' @param family_list List. List of marginal distributions
#' @keywords internal

link_coefs <- function(param_list, k, family_list){

  # get hermite coefficients
  g_i <- hermite_coefs(param_list[[1]], k, family_list[[1]])
  g_j <- hermite_coefs(param_list[[2]], k, family_list[[2]])

  # Polynomial of link function
  l_ij <- g_i*g_j*factorial(1:k)

  return(l_ij)
}
