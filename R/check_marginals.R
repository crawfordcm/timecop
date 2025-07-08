#' Check that user-specified marginal distributions are supported
#'
#' @param family List. A list of length d containing the names of marginal distributions
#' @param d Numeric. The number of variables
#' @keywords internal

check_marginals <- function(family, d) {

  if (is.null(family)) {
    stop("User must provide a list of marginal distributions", call. = FALSE)
  }

  if (length(family) != d) {
    stop("Number of families must match number of variables", call. = FALSE)
  }

  # supported marginals
  supported <- list("Bernoulli", "Poisson", "Gaussian")

  # indices of incorrect marginals
  bad_idx <- which(!family %in% supported)

  if (length(bad_idx) > 0) {
    bad <- unique(unlist(family[bad_idx]))

    stop(
      sprintf(
        "Distribution%s not supported: %s\nAllowed: %s",
        if (length(bad) > 1) "s" else "",
        paste(bad, collapse = ", "),
        paste(unlist(supported), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}
