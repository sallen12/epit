#' Check forecast lag
#'
#' @param h forecast lag.
#'
#' @return
#' Error message if \code{h} is not a single positive integerl.
#'
#' @keywords internal
check_h <- function(h) {
  if (!(is.numeric(h) && (length(h) == 1) && (floor(h) == h)))
    stop("'h' should be a single positive integer")
}

#' Check PIT
#'
#' @param z PIT or quantile PIT.
#'
#' @return
#' Error message if \code{z} is not a vector in [0,1] without missing values.
#'
#' @keywords internal
check_z <- function(z) {
  if (!(is.vector(z, "numeric") & (min(z) >= 0) & (max(z) <= 1)))
    stop("(quantile) PIT should be numeric vector in [0,1]")
}

#' Check ranks
#'
#' @param r ranks for rank histogram.
#' @param m ensemble size.
#'
#' @return
#' Error message if \code{r} is not a vector of positive integers smaller or
#' equal to \code{m+1}. Warning if the maximum of \code{r} is stritcly less than
#' \code{m+1}.
#'
#' @keywords internal
check_ranks <- function(r, m) {
  if (!(is.numeric(m) && (length(m) == 1) && (floor(m) == m)))
    stop("'m' should be a single positive integer")
  if (!(is.vector(r, "numeric") & (min(r) > 0) & !any(r != floor(r))))
    stop("'r' PIT should be a vector of positive integers")
  mr <- max(r)
  if (mr > m + 1)
    stop("maximum of 'r' should be less or equal to m+1")
  if (mr < m + 1)
    warning("maximum of r is strictly smaller than m+1. Is 'm' correct?")
}

#' Check strategy
#'
#' @param strategy chosen strategy for computing e-values.
#' @param type \code{"pit"}, \code{"quantile_pit"}, or \code{"rank_histogram"}.
#'
#' @return
#' Error message if \code{"strategy"} is not among the available strategies
#'
#' @keywords internal
check_strategy <- function(strategy, type) {
  if (type == "pit") {
    if (!identical(strategy, "beta") & !identical(strategy, "kernel"))
      stop("'strategy' must be 'beta' or 'kernel'")
  } else if (type == "quantile_pit") {
    if (!identical(strategy, "grenander") & !identical(strategy, "bernstein"))
      stop("'strategy' must be 'grenander' or 'bernstein'")
  } else if (type == "rank_histogram") {
    if (!identical(strategy, "empirical") & !identical(strategy, "betabinom"))
      stop("'strategy' must be 'betabinom' or 'empirical'")
  }
}

#' Check forecasts
#'
#' @param y vector of forecasts.
#' @param mat matrix of forecast quantiles or ensemble forecasts.
#'
#' @return
#' Error message if objects are not of specified types or contain NA or sizes
#' do not match.
#'
#' @keywords internal
check_forecasts <- function(y, mat) {
  if (!(is.vector(y, "numeric")) || anyNA(y))
    stop("'y' should be a numeric vector without NA")
  if (!is.matrix(mat) || !is.numeric(mat) || anyNA(mat))
    stop("forecasts should be a numeric matrix without NA")
  if (length(y) != nrow(mat))
    stop("length of forecasts and observations do not match")
}
