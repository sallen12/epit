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

#' Quantile PIT
#'
#' Computes the quantile PIT for forecasts specified by a collection of
#' quantiles.
#'
#' @param y observations (numeric vector).
#' @param quantiles quantile forecasts. A numeric matrix, the number of rows
#'     should equal the length of \code{y}.
#' @param alphas levels of the quantiles. Either a numeric matrix of the same
#'     dimension as \code{quantiles}, giving for each quantile the level
#'     (for example 0.2 for the 20% quantile). If the quantiles are equispaced
#'     (for example, always the 0.1, 0.2, ..., 0.9 or 0.05, 0.10, ..., 0.90,
#'     0.95 quantiles for all forecasts), then \code{alphas} can be specified
#'     as a single number giving the grid spacing of the quantiles (0.1 and
#'     0.05, respectively, in the examples before).
#' @param seed seed for random number generation (the PIT are randomized). Can
#'     be omitted.
#'
#' @return
#' A list containing the upper and lower quantile PIT (\code{zu} and \code{zl}).
#'
#' @importFrom stats runif
#'
#' @export
quantile_pit <- function(y, quantiles, alphas, seed = NULL) {
  if (!is.vector(y, "numeric"))
    stop("'y' must be a numeric vector")
  if (!is.matrix(quantiles) && nrow(quantiles) == length(y))
    stop("'quantiles' must be a numeric matrix nrow(quantiles) == length(y)")
  if (!(is.vector(alphas, "numeric") && length(alphas) == 1)) {
    if (!(is.matrix(alphas) & identical(dim(alphas), dim(quantiles))))
      stop("'alphas' must be a single number or of the same dimension as 'quantiles'")
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

  if (!is.null(seed)) {
    set.seed(seed)
  }
  V <- stats::runif(n)

  if (is.matrix(alphas)) {
    alphas <- cbind(0, alphas, 1)
    nseq <- seq_len(n)
    F_u <- alphas[cbind(nseq, pos["ind", ])]
    F_u_minus <- alphas[cbind(nseq, pos["ind_minus", ])]
    pos <- pos + 1
    F_l <- alphas[cbind(nseq, pos["ind", ])]
    F_l_minus <- alphas[cbind(nseq, pos["ind_minus", ])]
    zu <- F_u_minus + V * (F_u - F_u_minus)
    zl <- F_l_minus + V * (F_l - F_l_minus)
  } else {
    pos <- pos - 1
    zu <- pos["ind_minus", ]
    zu <- (zu + V * (pos["ind", ] - zu)) * alphas
    zl <- zu + alphas
  }
  list(zu = zu, zl = zl)
}

#' Ranks for rank histogram
#'
#' Computes the rank of an observation among an ensembles of forecasts.
#'
#' @param y observations (vector).
#' @param ensemble ensemble forecasts (matrix with \code{length(y)} rows).
#'
#' @details
#' In case of ties, the ranks are randomized.
#'
#' @return
#' Vector of ranks.
#'
#' @export
ensemble_rank <- function(y, ensemble) {
  apply(
    cbind(y, ensemble), 1, function(row) rank(row, ties.method = "random")[1]
  )
}

#' Combine lag h e-values
#'
#' Combine e-values for forecasts with lag \code{h > 1}.
#'
#' @param es list of e-values.
#'
#' @details
#' It is assumed that the first enrty of \code{es} contains all e-values with
#' indices \code{1,1+h,1+2h,...}, and more general the \code{j}th entry contains
#' the e-values with indices \code{j,j+h,j+2h,...}.
#'
#' @return
#' Combine e-value.
combine_e_h <- function(es) {
  h <- length(es)
  ns <- lengths(es)
  n <- sum(ns)
  e <- numeric(n)
  pos <- (seq_len(n) %% h + 1L)
  for (j in seq_len(h)) {
    tmp <- cumprod(es[[j]])
    e[pos] <- c(tmp[1], diff(tmp))
  }
  e <- cumsum(e) / h
}

#' Merge e-values and perform hypothesis test
#'
#' Merges a vector of e-values to a single e-value and performs a hypothesis
#' test (if desired).
#'
#' @param e_output output of functions \code{\link{e_pit}},
#'     \code{\link{e_quantile_pit}}, or \code{\link{e_rank_histogram}}.
#'
#' @details
#' If the lag in the e-values is 1, e-values are merged by (cumulative) product.
#' The null hypothesis can be rejected at level \code{\alpha}if the cumulative
#' product exceeds the level \code{1/alpha} at least once.
#'
#' For lag h, all e-values with time lag h are merged by product, and the
#' average is taken. The null hypothesis can be rejected if this process exceeds
#' \code{1/alpha} at least once.
#'
#' @return
#' A list containing the maximum of the cumulated e-values (\code{e}), the index
#' where this maximum is attained (\code{max_ind}), and the e-value at the end
#' of the observation period (\code{e_end}).
#'
#' @export
merge_e <- function(e_output) {
  if (e_output$h > 1) {
    e <- combine_e_h(lapply(e_output$evalues_h, function(x) x$e))
  } else {
    e <- cumprod(e_output$e)
  }
  max_ind <- which.max(e)[1]
  list(e = e[max_ind], max_ind = max_ind, end = e[length(e)])
}
