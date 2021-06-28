#' Correct e-value so that it never equals zero
#'
#' Transforms an e-value \code{e} so that it is never zero, by replacing
#' \code{e[i]} with \code{1/i+(1-1/i)*e[i]}.
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

#' Combine lag h e-values
#'
#' Combine e-values for forecasts with lag \code{h>1}.
#'
#' @param es list of e-values.
#'
#' @details
#' It is assumed that the first enrty of \code{es} contains all e-values with
#' indices \code{1,1+h,1+2h,...}, and more general the \code{j}th entry contains
#' the e-values with indices \code{j,j+h,j+2h,...}.
#'
#' @return
#' Combined e-value (vector of length \code{sum(lengths(es))}).
evalue_combine_h <- function(es) {
  h <- length(es)
  ns <- lengths(es)
  n <- sum(ns)
  e <- rep(exp(cumsum(log(es[[1]]))), each = h, length.out = n)
  for (j in 2:h) {
    e <- e +
      c(rep(0, j - 1), rep(exp(cumsum(log(es[[j]]))), each = h))[seq_len(n)]
  }
  e[seq_len(h)] <- 1
  e / h
}

#' Merge e-values to perform hypothesis tests
#'
#' Merges a vector of e-values to a single e-value to perform hypothesis tests.
#'
#' @param e_output output of functions \code{\link{e_pit}},
#'     \code{\link{e_quantile_pit}}, or \code{\link{e_rank_histogram}}.
#'
#' @details
#' If the lag in the e-values is 1, e-values are merged by (cumulative) product.
#' The null hypothesis can be rejected at level \code{alpha} if the cumulative
#' product exceeds the level \code{1/alpha} at least once.
#'
#' For lag \code{h>1}, all e-values with time lag \code{h} are merged by
#' product, and the average is taken. The null hypothesis can be rejected if
#' this process exceeds \code{1/alpha} at least once.
#'
#' @return
#' A list containing the maximum of the cumulated e-values (\code{e}), the index
#' where this maximum is attained (\code{max_ind}), and the e-value at the end
#' of the observation period (\code{e_end}).
#'
#' @export
evalue_merge <- function(e_output) {
  if (e_output$h > 1) {
    e <- evalue_combine_h(lapply(e_output$evalues_h, function(x) x$e))
  } else {
    e <- cumprod(e_output$e)
  }
  max_ind <- which.max(e)[1]
  list(e = e[max_ind], max_ind = max_ind, e_end = e[length(e)])
}

#' Quantile PIT
#'
#' Computes the quantile PIT for forecasts specified by a collection of
#' quantiles.
#'
#' @param y observations (numeric vector).
#' @param quantiles quantile forecasts. A numeric matrix, the number of rows
#'     should equal the length of \code{y}. It is assumed that the \code{j}th
#'     column contains the \code{j/(K+1)}th quantile, where
#'     \code{K=ncol(quantiles)}.
#' @param seed seed for random number generation (the PITs are randomized). Can
#'     be omitted.
#' @param check check for correct format of input parameters.
#'
#' @return
#' A list containing the upper and lower quantile PIT (\code{zu} and \code{zl}).
#'
#' @importFrom stats runif
#'
#' @export
quantile_pit <- function(
    y,
    quantiles,
    seed = NULL,
    check = FALSE
  ) {
  if (check) {
    check_forecasts(y, quantiles)
  }
  n <- length(y)
  quantiles <- cbind(-Inf, quantiles)
  pos <- mapply(
    function(obs, qs) {
      c(
        ind = findInterval(x = obs, vec = qs),
        ind_minus = findInterval(x = obs, vec = qs, left.open = TRUE)
      )
    },
    obs = y,
    qs = asplit(quantiles, 1)
  )
  alpha <- 1 / ncol(quantiles)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  V <- stats::runif(n)
  pos <- pos - 1
  zu <- pos["ind_minus", ]
  zu <- (zu + V * (pos["ind", ] - zu)) * alpha
  zl <- zu + alpha
  list(zu = zu, zl = zl)
}

#' Ranks for rank histogram
#'
#' Computes the rank of an observation among an ensembles of forecasts.
#'
#' @param y observations (vector).
#' @param ensemble ensemble forecasts (matrix with \code{length(y)} rows).
#' @param check check for correct format of input parameters.
#'
#' @details
#' In case of ties, the ranks are randomized.
#'
#' @return
#' Vector of ranks.
#'
#' @export
ensemble_rank <- function(y, ensemble, check = FALSE) {
  if (check) {
    check_forecasts(y, ensemble)
  }
  apply(
    cbind(y, ensemble), 1, function(row) rank(row, ties.method = "random")[1]
  )
}
