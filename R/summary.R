#' Summarize a fitted timecop model
#'
#' Produces a summary of the latent Gaussian VAR estimation results,
#' including coefficient estimates, standard errors, z-values, and p-values
#' for each equation.
#'
#' @param object A \code{timecop_fit} object returned by \code{\link{fit_timecop}}.
#' @param ... Additional arguments (currently unused).
#' @return An object of class \code{"summary.timecop_fit"} containing:
#'   \item{coefficients}{A list of coefficient matrices, one per equation.
#'     Each matrix has columns: Estimate, Std. Error, z value, Pr(>|z|).}
#'   \item{varnames}{Character vector of variable names.}
#'   \item{d}{Number of endogenous variables.}
#'   \item{n}{Time series length.}
#'   \item{p}{VAR order.}
#'   \item{family}{List of marginal distributions.}
#'   \item{corr}{Logical indicating whether correlations were used.}
#' @export
summary.timecop_fit <- function(object, ...) {
  est <- object$estimates
  se  <- object$se
  obj <- object$obj
  d   <- obj@d

  varnames <- rownames(obj@data)
  if (is.null(varnames)) {
    varnames <- paste0("Y", seq_len(d))
  }

  coefficients <- vector("list", d)
  names(coefficients) <- varnames

  for (i in seq_len(d)) {
    est_row <- est[i, ]
    se_row  <- se[i, ]
    z_row   <- est_row / se_row
    p_row   <- 2 * pnorm(abs(z_row), lower.tail = FALSE)

    coef_table <- cbind(
      Estimate     = est_row,
      `Std. Error` = se_row,
      `z value`    = z_row,
      `Pr(>|z|)`   = p_row
    )
    rownames(coef_table) <- paste0(varnames, ".l1")
    coefficients[[i]] <- coef_table
  }

  result <- list(
    coefficients = coefficients,
    varnames     = varnames,
    d            = d,
    n            = obj@n,
    p            = obj@p,
    family       = obj@family,
    corr         = obj@corr
  )
  class(result) <- "summary.timecop_fit"
  result
}

#' Print a timecop model summary
#'
#' @param x A \code{summary.timecop_fit} object.
#' @param ... Additional arguments passed to \code{\link[stats]{printCoefmat}}.
#' @return Invisibly returns \code{x}.
#' @export
print.summary.timecop_fit <- function(x, ...) {
  cat("\nLatent Gaussian VAR(1) Estimation Results:\n")
  cat("==========================================\n")
  cat("Endogenous variables:", paste(x$varnames, collapse = ", "), "\n")
  cat("Marginal distributions:", paste(x$family, collapse = ", "), "\n")
  cat("Sample size:", x$n, "\n")
  if (x$corr) cat("Estimation based on correlations\n")
  cat("\n")

  for (i in seq_len(x$d)) {
    nm <- x$varnames[i]
    cat(sprintf("Estimation results for equation %s:\n", nm))
    cat(paste(rep("=", nchar(nm) + 36), collapse = ""), "\n")

    predictors <- paste0(x$varnames, ".l1")
    cat(nm, "=", paste(predictors, collapse = " + "), "\n\n")

    printCoefmat(x$coefficients[[i]], P.values = TRUE, has.Pvalue = TRUE,
                 signif.stars = TRUE, ...)
    cat("\n")
  }

  invisible(x)
}

#' Print a fitted timecop model
#'
#' @param x A \code{timecop_fit} object returned by \code{\link{fit_timecop}}.
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns \code{x}.
#' @export
print.timecop_fit <- function(x, ...) {
  cat("Latent Gaussian VAR(1) fit\n\n")
  cat("Coefficient matrix:\n")
  print(x$estimates)
  cat("\nUse summary() for detailed results with standard errors and p-values.\n")
  invisible(x)
}

#' Summarize a fitted graphical VAR model
#'
#' Produces a summary of a latent graphical VAR fit, including the selected
#' tuning parameters, penalty type, and edge counts for the temporal and
#' contemporaneous networks.
#'
#' @param object A \code{timecop_gvar} object returned by
#'   \code{\link{fit_graphical_var}}.
#' @param ... Additional arguments (currently unused).
#' @return An object of class \code{"summary.timecop_gvar"}.
#' @export
summary.timecop_gvar <- function(object, ...) {
  obj <- object$obj
  d   <- obj@d

  varnames <- rownames(obj@data)
  if (is.null(varnames)) varnames <- paste0("Y", seq_len(d))

  A     <- object$A_hat
  Omega <- object$Omega_hat
  Omega_off <- Omega
  diag(Omega_off) <- 0

  result <- list(
    A_hat          = A,
    Omega_hat      = Omega,
    varnames       = varnames,
    d              = d,
    n              = obj@n,
    family         = obj@family,
    penalty        = object$penalty,
    lambda_A       = object$lambda_A,
    lambda_Omega   = object$lambda_Omega,
    ebic           = object$ebic,
    gamma_ebic     = object$gamma_ebic,
    n_edges_A      = sum(abs(A) > 1e-8),
    n_possible_A   = d * d,
    n_edges_Omega  = sum(abs(Omega_off[upper.tri(Omega_off)]) > 1e-8),
    n_possible_Omega = d * (d - 1L) / 2
  )
  class(result) <- "summary.timecop_gvar"
  result
}

#' Print a graphical VAR model summary
#'
#' @param x A \code{summary.timecop_gvar} object.
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns \code{x}.
#' @export
print.summary.timecop_gvar <- function(x, ...) {
  cat("\nLatent Graphical VAR(1) Estimation Results:\n")
  cat("===========================================\n")
  cat("Endogenous variables:  ", paste(x$varnames, collapse = ", "), "\n")
  cat("Marginal distributions:", paste(x$family, collapse = ", "), "\n")
  cat("Sample size:           ", x$n, "\n")
  cat("Penalty:               ", x$penalty, "\n")
  cat(sprintf("Selected lambda_A:      %.5f\n", x$lambda_A))
  cat(sprintf("Selected lambda_Omega:  %.5f\n", x$lambda_Omega))
  cat(sprintf("EBIC (gamma = %.2f):     %.3f\n", x$gamma_ebic, x$ebic))
  cat("\n")
  cat(sprintf("Temporal network (A):          %d / %d edges\n",
              x$n_edges_A, x$n_possible_A))
  cat(sprintf("Contemporaneous network (Omega): %d / %d edges\n",
              x$n_edges_Omega, x$n_possible_Omega))

  cat("\nTemporal coefficient matrix (A):\n")
  print(round(x$A_hat, 4))
  cat("\nInnovation precision matrix (Omega):\n")
  print(round(x$Omega_hat, 4))

  invisible(x)
}

#' Print a fitted graphical VAR model
#'
#' @param x A \code{timecop_gvar} object returned by
#'   \code{\link{fit_graphical_var}}.
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns \code{x}.
#' @export
print.timecop_gvar <- function(x, ...) {
  d <- x$obj@d
  n_edges_A     <- sum(abs(x$A_hat) > 1e-8)
  Omega_off     <- x$Omega_hat
  diag(Omega_off) <- 0
  n_edges_Omega <- sum(abs(Omega_off[upper.tri(Omega_off)]) > 1e-8)

  cat("Latent Graphical VAR(1) fit\n\n")
  cat("Penalty:              ", x$penalty, "\n")
  cat(sprintf("Selected lambda_A:     %.5f\n", x$lambda_A))
  cat(sprintf("Selected lambda_Omega: %.5f\n", x$lambda_Omega))
  cat(sprintf("Temporal edges:        %d / %d\n", n_edges_A, d * d))
  cat(sprintf("Contemporaneous edges: %d / %d\n",
              n_edges_Omega, d * (d - 1L) / 2))
  cat("\nUse summary() for the estimated network matrices.\n")
  invisible(x)
}
