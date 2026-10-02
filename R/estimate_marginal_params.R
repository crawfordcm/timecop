#' Estimate marginal parameters from data
#'
#' @param data Matrix. A d by n multivariate time series matrix
#' @param family List. A list of length d with the names of each distribution
#' @param ordinal_levels List. A list of length d containing the numeric levels
#'   for each Ordinal variable and `NULL` for each non-Ordinal variable.
#' @param d Numeric. The number of variables
#' @return A list of length d of estimated marginal parameters.
#' @keywords internal

estimate_marginal_params <- function(data, family, ordinal_levels, d) {

  param_hat <- vector("list", d)

  for (i in seq_len(d)) {
    if (family[[i]] %in% c("Bernoulli", "Poisson")) {

      param_hat[[i]] <- mean(data[i, ], na.rm = TRUE)

    } else if (family[[i]] == "Ordinal") {

      counts <- tabulate(match(data[i, ], ordinal_levels[[i]]),
                         nbins = length(ordinal_levels[[i]]))

      if (any(counts == 0L)) {
        unobserved <- ordinal_levels[[i]][counts == 0L]

        stop(
          sprintf(
            paste0(
              "Cannot estimate ordinal marginal probabilities for variable %d. ",
              "The following supplied levels were not observed in the data: %s"
            ),
            i,
            paste(unobserved, collapse = ", ")
          ),
          call. = FALSE
        )
      }

      # maybe add check here for small probabilities
      # can suggest to user to collapse levels
      param_hat[[i]] <- counts / sum(counts)

    } else if (family[[i]] == "Gaussian") {

      param_hat[[i]] <- c(mean(data[i, ], na.rm = TRUE),
                          var(data[i, ],  na.rm = TRUE))
    }
  }
  param_hat
}
