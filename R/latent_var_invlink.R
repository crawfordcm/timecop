#' Compute latent Gaussian covariances
#'
#' @param cov_x_hat Array. An array of observed covariance matrices
#' @param d Numeric. Number of variables
#' @param p Numeric. The VAR order. Default is 1
#' @param ell_ij_hat Array. An array of link function coefficients
#' @keywords internal

latent_var_invlink <- function(cov_x_hat, d, p, ell_ij_hat){

  grid_u <- c(seq(-1,-0.71,0.01), seq(-0.7,0,0.2), seq(0.1,0.7,0.2), seq(0.71,1,0.01))

  cov_z_hat <- array(NA, dim(cov_x_hat))
  for (i in seq_len(d)) {
    for (j in seq_len(d)){
      cov_z_hat[i,j,] <- interpolation(ell_ij_hat[,i,j], grid_u, cov_x_hat[i,j,])
    }
  }

  return(cov_z_hat)

}
