#' Check and setup data
#'
#' @param data Matrix.  An n (time points) by d (variables) multivariate time series matrix
#' @keywords internal

setup_data <- function(data) {

  if (is.null(data)) {
    stop("Data matrix is not supplied", call. = FALSE)
  }

  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }

  if (!is.numeric(data)) {
    stop("Data must be numeric", call. = FALSE)
  }

  # transpose to conform with estimation code
  data <- t(data)
  return(data)

}
