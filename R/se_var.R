#' Compute standard errors
#'
#' @param data Matrix. A d by n multivariate time series matrix
#' @param gamma_hat Matrix. Stacked covariances matrix
#' @param Gamma_hat Matrix. Toeplitz covariance matrix
#' @param cov_x_hat Array. An array of observed covariances
#' @param d Numeric. The number of variables
#' @param p Numeric. The VAR order
#' @param n Numeric. Time series length
#' @param family List. A list of marginal distributions
#' @param marg_num Numeric. The total number of marginal parameters
#' @keywords internal

se_var <- function(data,
                   gamma_hat,
                   Gamma_hat,
                   cov_x_hat,
                   d,
                   p,
                   n,
                   family,
                   marg_num) {

  jacob <- numderiv(data, cov_x_hat, d, p, n, family)

  # derivs wrt marginal params
  f1 <- jacob[1:d^2, 1:marg_num]
  f2 <- jacob[(1 + d^2):(2*d^2), 1:marg_num]

  # derivs wrt to count covs
  if (d == 1) {
    # case for d == 1
    Gamma_deriv <- jacob[(1 + d^2):(2*d^2), (marg_num + 1):(marg_num + d^2)]
    gamma_deriv <- jacob[1:d^2, (d^2 + d + 1):(2*d^2 + d)]

  } else {
    # case for d > 1
    Gamma_deriv <- revVec( diag( jacob[(1 + d^2):(2*d^2), (marg_num + 1):(marg_num + d^2)] ) )
    gamma_deriv <- revVec( diag( jacob[1:d^2, (d^2 + marg_num + 1):(2*d^2 + marg_num)] ) )

    diag(Gamma_deriv) <- 0
  }

  # set up needed matrices
  S1_t <- cbind(diag(d), matrix(0,d,d))
  S2_t <- cbind(diag(d)[,d:1], matrix(0,d,d))
  S3_t <- cbind(matrix(0,d,d), diag(d))

  A1 <- (diag(d) %x% solve(Gamma_hat)) %*% (1 %x% f1)
  A2 <- -(t(gamma_hat) %x% diag(d)) %*% (t(solve(Gamma_hat)) %x% solve(Gamma_hat)) %*% (1 %x% f2)
  B  <- diag(d) %x% solve(Gamma_hat)
  C  <- (S3_t %x% S1_t)
  D <- -(t(gamma_hat) %x% diag(d)) %*% (t(solve(Gamma_hat)) %x% solve(Gamma_hat))
  E  <- S1_t %x% S1_t
  J  <- matrix(1,marg_num,1)

  Q1 <- rbind(vec(gamma_deriv),vec(Gamma_deriv))
  Q2 <- rbind(C,E)

  Sigma <- longrun_var(data, d, n)

  sigma11 <- Sigma[1:marg_num,1:marg_num]
  sigma21 <- Sigma[(marg_num+1):(marg_num+(p+1)^2*d^2),1:marg_num]
  sigma12 <- Sigma[1:marg_num,(marg_num+1):(marg_num+(p+1)^2*d^2)]
  sigma22 <- Sigma[(marg_num+1):(marg_num+(p+1)^2*d^2),(marg_num+1):(marg_num+(p+1)^2*d^2)]

  # put everything together
  variance <- (A1 + A2) %*% sigma11 %*% t(A1 + A2) +
    (A1 + A2) %*% ((J %*% t(Q1)) * (sigma12 %*% t(Q2))) %*% t(cbind(B,D)) +
    cbind(B,D) %*% ((Q1 %*% t(J)) * (Q2 %*% sigma21)) %*% t(A1 + A2) +
    cbind(B,D) %*% ((Q1 %*% t(Q1)) * (Q2 %*% sigma22 %*% t(Q2))) %*% t(cbind(B,D))

  se <- sqrt(diag(variance)/n)

  return(matrix(se, nrow = d))

}
