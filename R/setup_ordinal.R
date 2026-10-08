#' Setup ordinal levels
#'
#' @param data Matrix. A d by n multivariate time series matrix
#' @param family List. A list of length d with the names of each distribution
#' @param ordinal_levels List. A list containing one numeric vector of levels for each ordinal variable
#' @param d Numeric. The number of variables
#' @return A list of length d containing levels for each ordinal variable and NULL otherwise.
#' @keywords internal

setup_ordinal_levels <- function(data,
                                 family,
                                 ordinal_levels,
                                 d) {

  ordinal_idx <- which(
    vapply(family, identical, logical(1), "Ordinal")
  )

  n_ordinal <- length(ordinal_idx)

  if (n_ordinal == 0L) {
    if (!is.null(ordinal_levels) && length(ordinal_levels) > 0L) {
      stop(
        "'ordinal_levels' was supplied, but no variables use the Ordinal family",
        call. = FALSE
      )
    }

    return(vector("list", d))
  }

  if (!is.list(ordinal_levels) || length(ordinal_levels) != n_ordinal) {
    stop(
      sprintf(
        "'ordinal_levels' must be a list of length %d, one entry for each Ordinal variable",
        n_ordinal
      ),
      call. = FALSE
    )
  }

  levels_expanded <- vector("list", d)
  levels_expanded[ordinal_idx] <- ordinal_levels

  # check that the levels are correctly specified
  for (i in ordinal_idx) {
    levels_i <- levels_expanded[[i]]

    if (!is.numeric(levels_i) ||
        any(!is.finite(levels_i))) {
      stop(
        sprintf(
          "Ordinal levels for variable %d must be finite numeric values",
          i
        ),
        call. = FALSE
      )
    }

    if (anyDuplicated(levels_i) > 0L ||
        is.unsorted(levels_i, strictly = TRUE)) {
      stop(
        sprintf(
          "Ordinal levels for variable %d must be unique and strictly increasing",
          i
        ),
        call. = FALSE
      )
    }

    if (any(levels_i != floor(levels_i))) {
      stop(
        sprintf(
          "Ordinal levels for variable %d must be integers",
          i
        ),
        call. = FALSE
      )
    }

    if (length(levels_i) < 3L) {
      if (length(levels_i) == 2L) {
        stop(
          sprintf(
            'Variable %d has only two ordinal levels. Use the Bernoulli family instead and encode as 0 and 1.',
            i
          ),
          call. = FALSE
        )
      }

      stop(
        sprintf(
          "Ordinal variable %d must have at least three levels",
          i
        ),
        call. = FALSE
      )
    }

    if (any(diff(levels_i) != 1)) {
      stop(
        sprintf(
          "Ordinal levels for variable %d must be consecutive with no gaps",
          i
        ),
        call. = FALSE
      )
    }

    invalid <- setdiff(unique(data[i, ]), levels_i)

    if (length(invalid) > 0L) {
      stop(
        sprintf(
          "Ordinal variable %d contains values not listed in ordinal_levels: %s",
          i,
          paste(invalid, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  return(levels_expanded)
}
