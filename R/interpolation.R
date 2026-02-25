#' Compute inverse link function to get latent covariances
#'
#' @param coef Numeric. A vector of link function coefficients
#' @param u Numeric. A vector of grid values
#' @param v Numeric. A vector of observed covariances
#' @return A numeric vector of inverse link function values.
#' @keywords internal

interpolation <- function(coef,u,v){
  Spline_result <- nat_spline(coef,u);
  c <- Spline_result[[1]];
  d <- Spline_result[[2]];
  h <- Spline_result[[3]];

  knot <- vector("numeric",length=length(u))
  for (i in 1:length(u)){
    pow <- 1:length(coef)
    knot[i] <- coef[1:length(coef)]%*%(u[i]^pow)[1:length(coef)] #power series expansion
  }
  n <- length(knot);

  L_inv_v <- vector("numeric",length(v));
  for (i in 1:length(v)){
    # Cutoff L(u) < -1 or L(u) > 1:
    if ( v[i] < knot[1] |  v[i] > knot[n]){
      if ( v[i] < knot[1] ){
        # Cutoff L(u) < -1
        L_inv_v[i] <- -1;
        next;
      }
      else {
        # Cutoff L(u) > 1
        L_inv_v[i] <- 1;
        next;
      }
    }
    else{
      idx <- 1;
      while(v[i] > knot[idx]){
        idx <- idx+1;
      }
      idx <- idx-1;
      First_term <- d[(idx+1)]*(v[i]-knot[idx])^3/(6*h[(idx+1)])
      Second_term <- d[idx]*(knot[(idx+1)]-v[i])^3/(6*h[(idx+1)])
      Third_term <- c[1,idx]*(v[i]-knot[idx])
      Fourth_term <- c[2,idx]*(knot[(idx+1)]-v[i])

      L_inv_v[i] <- First_term + Second_term + Third_term + Fourth_term
      next;
    }
  }
  return(L_inv_v)
}
