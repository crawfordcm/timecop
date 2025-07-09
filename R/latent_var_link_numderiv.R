#' Latent link function for numderiv
#'
#' @param d Numeric. The number of variables
#' @param n Numeric. Time series length
#' @param param_hat List. List of estimated marginal parameters
#' @param family List. List of marginal distributions
#' @keywords internal

latent_var_link_numderiv <- function(d, n, param_hat, family){

  k <- 100

  ell_ij_hat <- array(NA,dim=c(k,d,d))

  for (i in 1:d) {
    param_list <- list()
    family_list <- list()

    param_list[[1]] <- param_hat[[i]]
    family_list[[1]] <- family[[i]]

    for (j in 1:d) {
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
