#' Simulate multivariate time series from a latent Gaussian VAR
#'
#' Generates a multivariate discrete-valued time series from a latent Gaussian
#' VAR process. The latent Gaussian series is transformed to observed marginals
#' (Bernoulli, Poisson, or Gaussian) via the probability integral transform.
#'
#' The latent process is generated as a stationary Gaussian VAR process and then
#' rescaled so that each latent component has marginal variance one. By default,
#' the latent innovations are independent. Users may instead supply either an
#' innovation covariance matrix through \code{sigma_lv} or an innovation
#' precision matrix through \code{omega_lv}. Supplying \code{omega_lv} is useful
#' for simulating a latent graphical VAR, where zeros in the precision matrix
#' encode contemporaneous conditional independencies among the innovations.
#'
#' @param d Numeric. The number of variables.
#' @param n Numeric. Time series length.
#' @param p Numeric. The VAR order. Currently only \code{p = 1} is supported.
#' @param param List. A list of length \code{d} containing marginal parameters.
#'   For Bernoulli variables, the success probability; for Poisson variables,
#'   the rate parameter. Gaussian variables do not currently use this value.
#' @param phi_lv Matrix. A \code{d} by \code{d} transition matrix for the latent
#'   VAR process.
#' @param family List. A list of length \code{d} with marginal distribution
#'   names. Supported options are \code{"Bernoulli"}, \code{"Poisson"}, and
#'   \code{"Gaussian"}.
#' @param sigma_lv Optional matrix. A \code{d} by \code{d} innovation covariance
#'   matrix for the latent VAR process. Default is \code{NULL}. Only one of
#'   \code{sigma_lv} and \code{omega_lv} may be supplied.
#' @param omega_lv Optional matrix. A \code{d} by \code{d} innovation precision
#'   matrix for the latent VAR process. Default is \code{NULL}. If supplied,
#'   the innovation covariance matrix is computed as \code{solve(omega_lv)}.
#'   Only one of \code{sigma_lv} and \code{omega_lv} may be supplied.
#'
#' @return A list with components:
#'   \item{X_t}{A \code{d} by \code{n} matrix of observed time series.}
#'   \item{Z_t}{A \code{d} by \code{n} matrix of latent Gaussian time series.}
#'   \item{A_true}{A \code{d} by \code{d} by \code{p} array of true scaled VAR
#'     coefficient matrices.}
#'   \item{Sigma_true}{A \code{d} by \code{d} scaled latent innovation covariance
#'     matrix used to simulate the latent innovations.}
#'   \item{Omega_true}{A \code{d} by \code{d} scaled latent innovation precision
#'     matrix, equal to \code{solve(Sigma_true)} when \code{omega_lv} is supplied;
#'     otherwise \code{NULL}.}
#'   \item{beta_true}{Vectorized true VAR coefficients.}
#'   \item{innovation_type}{Character string indicating how the latent innovation
#'     covariance was specified: \code{"identity"}, \code{"covariance"}, or
#'     \code{"precision"}.}
#'   \item{Y_vec}{Vectorized response for observed data.}
#'   \item{Z_bmat}{Block diagonal design matrix for observed data.}
#'   \item{Y_mat}{Response matrix for observed data.}
#'   \item{X_mat}{Predictor matrix for observed data.}
#'   \item{Y_vec_test}{Vectorized response for latent data.}
#'   \item{Z_bmat_test}{Block diagonal design matrix for latent data.}
#'   \item{Y_mat_test}{Response matrix for latent data.}
#'   \item{X_mat_test}{Predictor matrix for latent data.}
#'   \item{param}{Input marginal parameters.}
#'
#' @examples
#' sim <- latent_var_sim(
#'   d = 2, n = 100, p = 1,
#'   param = list(0.5, 0.3),
#'   phi_lv = matrix(c(0.3, 0.1, 0.1, 0.3), 2, 2),
#'   family = list("Bernoulli", "Bernoulli")
#' )
#' dim(sim$X_t)  # 2 x 100
#'
#' # With user-specified latent innovation covariance
#' sigma_lv <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
#' sim_cov <- latent_var_sim(
#'   d = 2, n = 100, p = 1,
#'   param = list(0.5, 0.3),
#'   phi_lv = matrix(c(0.3, 0.1, 0.1, 0.3), 2, 2),
#'   family = list("Bernoulli", "Bernoulli"),
#'   sigma_lv = sigma_lv
#' )
#'
#' # With user-specified latent innovation precision matrix
#' # for a latent graphical VAR
#' phi_lv <- matrix(
#'   c(
#'     0.30,  0.00,  0.10,
#'     0.20,  0.40,  0.00,
#'     0.00, -0.15,  0.25
#'   ),
#'   nrow = 3,
#'   byrow = TRUE
#' )
#'
#' omega_lv <- matrix(
#'   c(
#'     1.00, -0.30,  0.00,
#'    -0.30,  1.00,  0.25,
#'     0.00,  0.25,  1.00
#'   ),
#'   nrow = 3,
#'   byrow = TRUE
#' )
#'
#' sim_gvar <- latent_var_sim(
#'   d = 3, n = 100, p = 1,
#'   param = list(0.5, 3, 0.3),
#'   phi_lv = phi_lv,
#'   family = list("Bernoulli", "Poisson", "Bernoulli"),
#'   omega_lv = omega_lv
#' )
#'
#' @importFrom mvtnorm rmvnorm
#' @importFrom expm %^%
#' @importFrom Matrix bdiag
#' @export

latent_var_sim <- function(d, n, p, param, phi_lv, family,
                           sigma_lv = NULL, omega_lv = NULL){

  if (!is.numeric(d) || length(d) != 1 || d < 1 || d != as.integer(d)) {
    stop("'d' must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1 || n < 1 || n != as.integer(n)) {
    stop("'n' must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(p) || length(p) != 1 || p != 1) {
    stop("Only p=1 is currently supported", call. = FALSE)
  }
  if (!is.list(param) || length(param) != d) {
    stop("'param' must be a list of length d", call. = FALSE)
  }
  if (!is.matrix(phi_lv) || nrow(phi_lv) != d || ncol(phi_lv) != d) {
    stop("'phi_lv' must be a d x d matrix", call. = FALSE)
  }
  if (!is.list(family) || length(family) != d) {
    stop("'family' must be a list of length d", call. = FALSE)
  }

  if (!is.null(sigma_lv) && !is.null(omega_lv)) {
    stop("Provide only one of 'sigma_lv' or 'omega_lv', not both.", call. = FALSE)
  }

  if (is.null(sigma_lv) && is.null(omega_lv)) {

    Sigma_star <- diag(1, nrow = d)
    innovation_type <- "identity"

  } else if (!is.null(sigma_lv)) {

    if (!is.matrix(sigma_lv) || nrow(sigma_lv) != d || ncol(sigma_lv) != d) {
      stop("'sigma_lv' must be a d x d covariance matrix.", call. = FALSE)
    }

    if (max(abs(sigma_lv - t(sigma_lv))) > 1e-8) {
      stop("'sigma_lv' must be symmetric.", call. = FALSE)
    }

    eig_sigma <- eigen(sigma_lv, symmetric = TRUE, only.values = TRUE)$values
    if (min(eig_sigma) <= 0) {
      stop("'sigma_lv' must be positive definite.", call. = FALSE)
    }

    Sigma_star <- sigma_lv
    innovation_type <- "covariance"

  } else {

    if (!is.matrix(omega_lv) || nrow(omega_lv) != d || ncol(omega_lv) != d) {
      stop("'omega_lv' must be a d x d precision matrix.", call. = FALSE)
    }

    if (max(abs(omega_lv - t(omega_lv))) > 1e-8) {
      stop("'omega_lv' must be symmetric.", call. = FALSE)
    }

    eig_omega <- eigen(omega_lv, symmetric = TRUE, only.values = TRUE)$values
    if (min(eig_omega) <= 0) {
      stop("'omega_lv' must be positive definite.", call. = FALSE)
    }

    Sigma_star <- solve(omega_lv)
    innovation_type <- "precision"
  }


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
  while( diff > 1e-8){
    order <- order + 1
    if ( d == 1){
      update_mat <- C_star^order * Sigma_U_star * as(C_star, "TsparseMatrix")^order
    }else{
      update_mat <- as.matrix(C_star) %^% order %*% as.matrix(Sigma_U_star) %*% t((as.matrix(C_star)) %^% order)
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


  output <- list(
    X_t = X_t,
    Z_t = Z_t,

    A_true = A,
    Sigma_true = Sigma,
    beta_true = beta,

    Omega_true = if (!is.null(omega_lv)) solve(Sigma) else NULL,
    innovation_type = innovation_type,

    Y_vec_test = vec_bf_Y_test,
    Z_bmat_test = mat_cal_Z_test,
    Y_mat_test = Y_cal_Z,
    X_mat_test = X_cal_Z,

    param = param,

    Y_vec = vec_bf_Y,
    Z_bmat = mat_cal_Z,
    Y_mat = Y_cal_X,
    X_mat = X_cal_X,

    family = family
  )

  return(output)

}
