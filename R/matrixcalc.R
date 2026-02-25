#' Column vectorize a matrix
#' @param x A matrix.
#' @return A column vector.
#' @keywords internal
vec <- function(x) t(t(as.vector(x)))

#' Half-vectorize a symmetric matrix
#' @param x A symmetric matrix.
#' @return A column vector of lower triangular elements including diagonal.
#' @keywords internal
vech <- function(x) t(t(x[!upper.tri(x)]))

#' Strict lower triangular vectorization
#' @param x A square matrix.
#' @return A column vector of strictly lower triangular elements.
#' @keywords internal
vecp <- function(x) t(t(x[!upper.tri(x, diag = TRUE)]))

#' Diagonal vectorization
#' @param x A square matrix.
#' @return A column vector of diagonal elements.
#' @keywords internal
vecd <- function(x) t(t(diag(x)))

#' Reverse vec operator
#'
#' Reshapes a vector back into a matrix (column-major).
#'
#' @param x A numeric vector.
#' @param ncol Integer. Number of columns (optional for square matrices).
#' @param nrow Integer. Number of rows (optional for square matrices).
#' @return A matrix.
#' @keywords internal
revVec <- function(x, ncol, nrow){

  if (length(x)==1)  { return(x) }

  d <- sqrt(length(x))

  if (missing(ncol) | missing(nrow)){
    ncol <- nrow <- d
    if (round(d) != d){
      stop("timecop: Dimensions needed for non-square matrices.")
    }
  }

  revvecx <- matrix(0, nrow = nrow, ncol = ncol)

  for (j in 1:ncol){
    revvecx[,j] <- x[c(1:nrow) + (j-1)*nrow]
  }
  return(revvecx)
}

#' Reverse vech operator
#' @param x A numeric vector.
#' @return A symmetric matrix.
#' @keywords internal
revVech <- function(x){

  if (length(x)==1) { return(x) }

  d <- (-1 + sqrt(8*length(x) + 1))/2

  if (round(d) != d){
    stop("timecop: Matrix is not square.")
  }

  revvechx <- matrix(0, nrow=d, ncol=d)

  for (j in 1:d){
    revvechx[j:d,j] <- x[1:(d-j+1)+ (j-1)*(d - 1/2*(j-2))]
  }

  revvechx <- revvechx + t(revvechx) - diag(diag(revvechx))

  return(revvechx)
}

#' Reverse vecp operator
#' @param x A numeric vector.
#' @return A symmetric matrix with unit diagonal.
#' @keywords internal
revVecp <- function(x){

  d <- (1 + sqrt(8*length(x) + 1))/2

  if (round(d) != d){
    stop("timecop: Matrix is not square.")
  }
  revvecp <- matrix(0, nrow=d, ncol=d)
  revvecp[lower.tri(revvecp , diag=FALSE)] <- x
  diag(revvecp) <- 1
  revvecp <- revvecp + t(revvecp) - diag(diag(revvecp))

  return(revvecp)
}

#' Create unit basis vector e_i
#' @param i Integer. Position of the 1.
#' @param n Integer. Length of the vector.
#' @return A numeric vector of length n with 1 in position i.
#' @keywords internal
make_e_basis <- function(i, n) {

  replace(numeric(n), i, 1)

}

#' Create symmetric unit basis vector
#' @param i Integer. Row index.
#' @param j Integer. Column index.
#' @param n Integer. Matrix dimension.
#' @return A numeric unit vector of length n(n+1)/2.
#' @keywords internal
make_u_basis <- function(i,j,n) {

  replace(numeric(.5*n*(n+1)), (j-1)*n + i - .5*j*(j-1), 1)

}

#' Create strictly lower triangular unit basis vector
#' @param i Integer. Row index.
#' @param j Integer. Column index.
#' @param n Integer. Matrix dimension.
#' @return A numeric unit vector of length n(n-1)/2.
#' @keywords internal
make_v_basis <- function(i,j,n) {

  replace(numeric(.5*n*(n-1)), (j-1)*n + i - .5*j*(j+1), 1)

}

#' Commutation matrix
#'
#' Constructs the mn x mn commutation matrix.
#'
#' @param m Integer. First dimension.
#' @param n Integer. Second dimension (defaults to m).
#' @return An mn x mn commutation matrix.
#' @keywords internal
comMat <- function(m, n=m){
  K <- matrix(0, nrow = n*m, ncol = n*m)
  for(i in 1:m){
    for(j in 1:n){
      H <- t(t(make_e_basis(i,m))) %*% t(make_e_basis(j,n))
      K <- K + kronecker(H, t(H))
    }
  }
  return(K)
}

#' Symmetrizer matrix N
#' @param n Integer. Matrix dimension.
#' @return An n^2 x n^2 symmetrizer matrix.
#' @keywords internal
nMat <- function(n){
  return(.5 * (diag(n^2) + comMat(n)))
}

#' Elimination matrix
#'
#' @param n Integer. Matrix dimension.
#' @param type Character. Type of elimination: `"s"` (symmetric), `"l"`
#'   (strictly lower triangular), or `"d"` (diagonal).
#' @return An elimination matrix.
#' @keywords internal
eliMat <- function(n, type = "s"){


  if (type == "s"){

    M <- matrix(0, nrow= .5*n*(n+1), ncol = n^2)

    for(i in 1:n){
      for(j in 1:n){
        if( i >= j){

          M <- M + kronecker(kronecker(make_u_basis(i,j,n),t(make_e_basis(j,n))),t(make_e_basis(i,n)))

        }
      }
    }

  } else if (type == "l"){

    M <- matrix(0, nrow= .5*n*(n-1), ncol = n^2)

    for(i in 1:n){
      for(j in 1:n){
        if( i > j){

          M <- M + kronecker(kronecker(make_v_basis(i,j,n),t(make_e_basis(j,n))),t(make_e_basis(i,n)))

        }
      }
    }

  } else if (type == "d"){

    M <- matrix(0, nrow = n, ncol = n^2)

    for(i in 1:n){

      M <- M + kronecker(kronecker(make_e_basis(i,n), t(make_e_basis(i,n))), t(make_e_basis(i,n)))

    }

  } else {

    M <- NULL

  }

  return(M)

}

#' Duplication matrix
#'
#' @param n Integer. Matrix dimension.
#' @param type Character. Type of duplication: `"s"` (symmetric), `"l"`
#'   (strictly lower triangular), or `"d"` (diagonal).
#' @return A duplication matrix.
#' @keywords internal
dupMat <- function(n, type = "s"){

  if (type == "s"){

    Dt <- matrix(0, nrow = .5*n*(n+1), ncol= n^2)

    for(i in 1:n){
      for(j in 1:n){
        if( i >= j){
          Tm <- matrix(0,n,n); Tm[i,j] <- Tm[j,i] <- 1
          Dt <- Dt + make_u_basis(i,j,n) %*% t(vec(Tm))
        }
      }
    }

    return(t(Dt))

  } else if (type == "l"){

    Dt <- matrix(0, nrow = .5*n*(n-1), ncol= n^2)

    for(i in 1:n){
      for(j in 1:n){
        if( i > j){
          T1 <- matrix(0,n,n); T1[i,j] <- 1
          T2 <- matrix(0,n,n); T2[j,i] <- 1
          Tm <- T1-T2
          Dt <- Dt + make_v_basis(i,j,n) %*% t(vec(Tm))
        }
      }
    }

    return(t(Dt))

  } else if( type == "d"){

    Dt <- eliMat(n, type = "d")
    return(t(Dt))


  }

}

#' Commutation matrix (transposed)
#'
#' @param n Integer. First dimension.
#' @param m Integer. Second dimension (defaults to n).
#' @return A transposed commutation matrix.
#' @keywords internal
makeCommutationMat <- function(n, m=n){

  K <- matrix(0, nrow = n*m, ncol = n*m)

  for(i in 1:m){
    for(j in 1:n){

      H <- t(t(make_e_basis(i,m))) %*% t(make_e_basis(j,n))
      K <- K + kronecker(H, t(H))

    }
  }
  return(t(K))
}
