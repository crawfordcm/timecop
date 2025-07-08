#' Simulate multivariate time series
#'
#' @param d Numeric. The number of variables
#' @param n Numeric. Time series length
#' @param p Numeric. The VAR order
#' @param param List. A list of marginal parameters
#' @param phi_lv Matrix. A d by d transition matrix
#' @param family List. A list of marginal distributions
#' @export

latent_var_sim <- function(d, n, p, param, phi_lv, family){

  Sigma_star <- diag(1,nrow=d)
  A_star <- array(phi_lv,dim=c(d,d,p))

  decay.rate <- 1
  test.mat <- companion_form_phi(A_star,d,p)
  test.eigen <- max(abs(eigen(test.mat)$values))
  if( test.eigen > 0.95 ){
    decay.rate <- 0.9
    while ( test.eigen > 0.95 ){
      for (h in 1:p){
        A_h_star <- A_star[,,h]
        A_h_star[upper.tri(A_h_star)] <- decay.rate*A_h_star[upper.tri(A_h_star)]
        A_h_star[lower.tri(A_h_star)] <- decay.rate*A_h_star[lower.tri(A_h_star)]
        A_star[,,h] <- A_h_star
      }
      test.mat <- companion_form_phi(A_star,d,p)
      decay.rate <- decay.rate*0.9
      test.eigen <- max(abs(eigen(test.mat)$values))
    }
  }
  A_star <- A_star*decay.rate

  if (p == 1){
    Sigma_U_star <- as(Sigma_star,"sparseMatrix")
    C_star <- as(A_star[,,1],"sparseMatrix")
  }else{
    Sigma_U_star <- matrix(0,p*d,p*d)
    Sigma_U_star[1:d,1:d] <- Sigma_star
    Sigma_U_star <- as(Sigma_U_star,"sparseMatrix")
    C_star <- as(companion_form_phi(A_star,d,p),"sparseMatrix")
  }

  diff <- 1
  order <- 0
  Sigma_Y_star <- Sigma_U_star
  while( diff > 1e-5){
    order <- order + 1
    if ( d == 1){
      #update_mat <- C_star^order * Sigma_U_star * (t(C_star))^order
      update_mat <- C_star^order * Sigma_U_star * as(C_star, "TsparseMatrix")^order
    }else{
      #update_mat <- C_star^order %*% Sigma_U_star %*% (t(C_star))^order
      update_mat <- C_star^order %*% Sigma_U_star %*% as(C_star, "TsparseMatrix")^order
    }
    Sigma_Y_star <- Sigma_Y_star + update_mat
    if ( d == 1){
      diff <- update_mat[1:d,1:d]
    }else{
      diff <- norm(as.matrix(update_mat[1:d,1:d]),type="o")
    }
  }
  if ( d == 1){
    Sigma_X_star <- Sigma_Y_star
  }else{
    Sigma_X_star <- diag(diag(as.matrix(Sigma_Y_star[1:d,1:d])))
  }

  A <- array(NA,dim=c(d,d,p))
  for ( h in 1:p){
    if ( d == 1){
      A[,,h] <- A_star[,,h]/as.numeric(Sigma_X_star)
    }else{
      A[,,h] <- as.matrix(solve(sqrt(Sigma_X_star))%*%A_star[,,h]%*%sqrt(Sigma_X_star))
    }
  }
  Sigma <- as.matrix(solve(sqrt(Sigma_X_star))%*%Sigma_U_star[1:d,1:d]%*%t(solve(sqrt(Sigma_X_star))))

  B <- matrix(NA,p*d,d)
  for (h in 1:p){
    B[((h-1)*d+1):(h*d),] <- t(A[,,h])
  }
  B <- as(B,"sparseMatrix")
  beta <- as.vector(B)


  Burn <- 500
  Z_t <- array(NA, dim=c(d,n+Burn))
  eps_t <- t(mvtnorm::rmvnorm(n+Burn, rep(0,d), Sigma, method="eigen"))

  for (t in 1:(n+Burn) ){
    if (t <= p){
      Z_t[,t] <- eps_t[,t]
    }else{
      if ( d == 1){
        Z_t[,t] <- sum(sapply(c(1:p),function(x){sum(A[,,x]%*%Z_t[,(t-x), drop = FALSE])})) + eps_t[,t]
      }else{
        Z_t[,t] <- rowSums(sapply(c(1:p),function(x){rowSums(A[,,x]%*%Z_t[,(t-x)])})) +eps_t[,t]
      }
    }
  }
  Z_t <- Z_t[,-c(1:Burn), drop = FALSE]

  # get count data
  X_t <- array(NA, dim=c(d,n))

  for (i in 1:d){
    # allow for different marginals
    if (family[[i]] == "Bernoulli"){
      Thrs <- qnorm(1 - param[[i]])
      X_t[i,] <- ( Z_t[i,] > Thrs )*1

    } else if (family[[i]] == "Poisson"){
      X_t[i,] <- qpois(pnorm(Z_t[i,]), param[[i]])

    } else if (family[[i]] == "Gaussian"){
      X_t[i,] <- Z_t[i,]

    }
  }

  Y_cal_X <- t(X_t[,(p+1):n])
  Y_cal_Z <- t(Z_t[,(p+1):n])
  if ( p == 1 ){
    X_cal_X <- t(X_t[,p:(n-1)])
    X_cal_Z <- t(Z_t[,p:(n-1)])
  }else{
    tmp_X <- t(X_t[,p:(n-1)])
    tmp_Z <- t(Z_t[,p:(n-1)])
    for (h in 2:p){
      tmp_X <- cbind(tmp_X,t(X_t[,(p-h+1):(n-h)]))
      tmp_Z <- cbind(tmp_Z,t(Z_t[,(p-h+1):(n-h)]))
    }
    X_cal_X <- tmp_X
    X_cal_Z <- tmp_Z
  }
  vec_bf_Y <- as.vector(Y_cal_X)
  mat_cal_Z <- bdiag(replicate(d,X_cal_X,simplify=FALSE))
  vec_bf_Y_test <- as.vector(Y_cal_Z)
  mat_cal_Z_test <- bdiag(replicate(d,X_cal_Z,simplify=FALSE))


  output <- list( X_t = X_t, Z_t = Z_t, A_true = A, Sigma_true = Sigma, beta_true = beta,
                  Y_vec_test = vec_bf_Y_test, Z_bmat_test = mat_cal_Z_test,
                  Y_mat_test = Y_cal_Z, X_mat_test = X_cal_Z, param = param,
                  Y_vec = vec_bf_Y, Z_bmat = mat_cal_Z,
                  Y_mat = Y_cal_X, X_mat = X_cal_X)
  return(output)
}
