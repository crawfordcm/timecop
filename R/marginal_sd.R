#' Compute a marginal standard deviation
#'
#' @param param Numeric. Marginal parameter or probability vector.
#' @param family Character. Name of the marginal distribution.
#' @param ordinal_levels Numeric or NULL. Supplied values for an Ordinal
#'   variable; `NULL` for other families.
#' @return A numeric marginal standard deviation.
#' @keywords internal

marginal_sd <- function(param,
                        family,
                        ordinal_levels = NULL) {

  if (family == "Bernoulli") {

    sd_value <- sqrt(param * (1 - param))

  } else if (family == "Poisson") {

    sd_value <- sqrt(param)

  } else if (family == "Ordinal") {

    if (length(param) != length(ordinal_levels)) {
      stop(
        sprintf(
          paste0(
            "For the Ordinal family, the number of probabilities (%d) ",
            "must equal the number of supplied levels (%d)."
          ),
          length(param),
          length(ordinal_levels)
        ),
        call. = FALSE
      )
    }

    if (any(!is.finite(param)) ||
        any(param < 0) ||
        abs(sum(param) - 1) > 1e-8) {
      stop(
        paste0(
          "For the Ordinal family, probabilities must be finite, ",
          "nonnegative, and sum to 1."
        ),
        call. = FALSE
      )
    }

    mean_value <- sum(ordinal_levels * param)

    variance <- sum(
      (ordinal_levels - mean_value)^2 * param
    )

    sd_value <- sqrt(variance)

  } else if (family == "Gaussian") {

    sd_value <- 1

  } else {

    stop(
      sprintf("Unsupported marginal family: %s", family),
      call. = FALSE
    )
  }

  if (length(sd_value) != 1L ||
      !is.finite(sd_value) ||
      sd_value <= 0) {
    stop(
      sprintf(
        "Cannot calculate a correlation for a %s variable with zero variance",
        family
      ),
      call. = FALSE
    )
  }

  return(sd_value)
}
