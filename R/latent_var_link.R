#' Latent link function
#'
#' @param data Matrix. A d by n multivariate time series matrix
#' @param d Numeric. The number of variables
#' @param n Numeric. Time series length
#' @param k Numeric. The value at which Hermite coefficient infinite sums terminate. Default is 100
#' @param family List. A list of length d with the names of each distribution
#' @keywords internal

latent_var_link <- function(data, d, n, k, family) {

  # estimate marginal parameters
  # need to change for more parameterized marginals
  param_hat <- list()
  for (i in seq_len(d)) {
    param_hat[[i]] <- sum(data[i,])/n
  }

  # Estimate coefficients of link function
  ell_ij_hat <- array(NA, dim = c(k,d,d))
  for (i in seq_len(d)) {
    #dist_pair <- vector("numeric",length=2)
    param_list <- list()
    family_list <- list()

    param_list[[1]] <- param_hat[[i]]
    family_list[[1]] <- family[[i]]

    for (j in seq_len(d)) {
      if (j >= i) {
        param_list[[2]] <- param_hat[[j]]
        family_list[[2]] <- family[[j]]
        ell_ij_hat[,i,j] <- link_coefs(param_list, k, family_list)
      } else {
        ell_ij_hat[,i,j] <- ell_ij_hat[,j,i]
      }
    }
  }
  return(ell_ij_hat)
}
