#' Compute splines
#'
#' @param coef Numeric. A vector of link function coefficients
#' @param u Numeric. A vector of grid values
#' @return A list with three components: spline coefficients matrix, second
#'   derivatives vector, and interval widths vector.
#' @keywords internal

nat_spline <- function(coef,u){

  v <- vector("numeric",length=length(u))
  for (i in 1:length(u)){
    pow <- 1:length(coef)
    v[i] <- coef[1:length(coef)]%*%(u[i]^pow)[1:length(coef)]
  }

  L_inv_v <- u;
  n <- length(v);

  h <- vector("numeric",length=n);
  LHS <- matrix(0,nrow=n,ncol=n);
  RHS <- vector("numeric",length=n);

  h[1] <- 0; h[n] <- v[n]-v[(n-1)];
  LHS[1,1] <- 1; LHS[n,n] <- 1;
  RHS[1] <- 0; RHS[n] <- 0;

  for (i in 2:(n-1)){
    h[i] <- v[i] - v[(i-1)];
    LHS[i,(i-1)] <- (v[i] - v[(i-1)])/6;
    LHS[i,i] <- ( (v[i] - v[(i-1)]) + (v[(i+1)] - v[i]) )/3;
    LHS[i,(i+1)] <- (v[(i+1)] - v[i])/6;
    RHS[i] <- ( L_inv_v[(i+1)]-L_inv_v[i] ) / (v[(i+1)] - v[i]) - ( L_inv_v[i]-L_inv_v[(i-1)] ) / (v[i] - v[(i-1)]);
  }

  c <- matrix(0,nrow=2,ncol=(n-1));
  d <- solve(LHS)%*%RHS;

  for (i in 1:(n-1)){
    c[1,i] <- L_inv_v[(i+1)]/h[(i+1)]-d[(i+1)]*h[(i+1)]/6;
    c[2,i] <- L_inv_v[i]/h[(i+1)]-d[i]*h[(i+1)]/6;
  }

  return(list(c,d,h));
}
