#' Correct e-value so that it never equals zero
#'
#' Transforms an e-value \code{e} so that it is never zero, by replacing
#' \code{e[i]} with \code{1/i + (1-1/i)*e[i]}.
#'
#' @param e numeric vector of e-values.
#'
#' @return
#' Vector of corrected e-values.
#'
#' @export
evalue_correct <- function(e) {
  inv_index <- 1 / seq_along(e)
  inv_index + (1 - inv_index) * e
}
