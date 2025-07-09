#' Compute long-run variance
#'
#' @param data Matrix. A d by n multivariate time series matrix
#' @param d Numeric. The number of variables
#' @param n Numeric. Time series length
#' @keywords internal

longrun_var <- function(data, d, n) {

  # center data
  dat_c <- data - rowMeans(data)

  F_list <- list()

  for(i in 2:(n-1)){
    xt <- matrix(dat_c[,i], d, 1)
    a1 <- xt %*% t(xt)
    a2 <- matrix(dat_c[,i-1], d, 1) %*% t(xt)
    a  <- vec(rbind(a1,a2))
    b1 <- matrix(dat_c[,i+1], d, 1) %*% t(xt)
    b2 <- xt %*% t(xt)
    b  <- vec(rbind(b1, b2))

    # for the Gaussian case (see AN_counts)
    if ("Gaussian" %in% family){ # note that this is set up for the case with d = 3 and 3rd component is Gaussian
      rt <- t(cbind(t(xt), matrix(dat_c[3,i]^2, 1, 1)))
      F_list[[i-1]] <- rbind(rt,a,b)
    } else {
      F_list[[i-1]] <- rbind(matrix(dat_c[,i-1], d, 1), a, b) #sol1
    }
    #F_list[[i-1]] <- rbind(xt,a,b) #sol 2
  }

  F_mat <- do.call(cbind, F_list)

  #sol1
  longrun <- getLongRunVar(t(F_mat), bandwidth = "nw", kernel = "ba", demeaning = TRUE, check = FALSE)
  Sigma   <- longrun$Omega

  return(Sigma)

}
