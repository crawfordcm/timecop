#' Check that user-specified marginal distributions are supported
#' @param family List. A list of length d containing the names of marginal distributions
#' @keywords internal

check_marginals <- function(family, d) {

  if (is.null(family)){
    stop("timecop error: user must provide a list of marginal distributions", call. = FALSE)
  }

  if (length(family) != d){
    stop("timecop error: number of families must match number of variables", call. = FALSE)
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
