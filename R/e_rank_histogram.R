#' E-values for rank histograms
#'
#' Tests null hypothesis that rank histogram follows a uniform distribution.
#'
#' @param r observed ranks (numbers in \code{1,...,m}).
#' @param h forecast lag.
#' @param m size of the ensemble PLUS ONE (positive integer, the \code{m} above).
#' @param strategy strategy for evaluating calibration. Available are
#'     \code{"empirical"} for the empirical distribution, and \code{"betabinom"}
#'     for the beta-binomial distribution (estimated with maximum likelihood).
#' @param options options for the given strategy (see \code{\link{betabinom_e}}
#'     and \code{\link{empirical_e}} for the available parameters).
#' @param check check for correct format of input parameters.
#'
#' @details
#' It is recommended to use \code{strategy="empirical"} only for large sample
#' sizes (>1000).
#'
#' If r is a matrix, then the rows will be treated as different time points, and
#' the columns as different variables to aggregate over. This accounts for cases
#' where more than one forecast and observation is available at each point in time,
#' e.g. from multiple points in space.
#'
#' @return
#' If \code{h} equals 1: Returns a list containing the vector of e-values
#' (\code{e}) and the forecast lag \code{h}.
#'
#' If \code{h} is greater than 1: Instead of a vector of e-values, the list
#' contains for each \code{j=1,2,...,h} a list with the e-values for all
#' observations with indices \code{j,j+h,j+2h,...}.
#'
#' @seealso
#' \code{\link{ensemble_rank}} to compute \code{r} from ensemble forecasts and
#' observations.
#'
#' @export
#'
#' @examples
#' n <- 360
#' r <- simulate_pit(n, "rank_histogram", K = 20, bias = 0.2, dispersion = 0)$r
#'
#' # If z is a lag 1 forecast:
#' e <- e_rank_histogram(r, h = 1, m = 20, strategy = "betabinom")
#' evalue_merge(e)
#' prod(e$e)
#' max(cumprod(e$e))
#'
#' # Lag 2:
#' e <- e_rank_histogram(r, h = 2, m = 20, strategy = "betabinom")
#' str(e)
#' evalue_merge(e)
#'
#'
#' # Multivariate
#' d <- 3
#' r <- simulate_pit(n*d, "rank_histogram", K = 20, bias = 0.2, dispersion = 0)$r
#' r <- matrix(r, ncol = d)
#' e <- e_rank_histogram(r, h = 2, m = 20, strategy = "betabinom")
#' str(e)
#' evalue_merge(e)
e_rank_histogram <- function(
    r,
    h,
    m,
    strategy = "betabinom",
    options = list(),
    check = FALSE
  ) {
  if (check) {
    check_ranks(r, m)
    check_h(h)
    check_strategy(strategy, "rank_histogram")
  }
  if (is.matrix(r)) {
    n <- nrow(r)
    e_func <- get(paste0(strategy, "_e_agg"))
  } else {
    n <- length(r)
    e_func <- get(paste0(strategy, "_e"))
  }
  if (h == 1) {
    e <- rep(1, n)
    if (is.matrix(r)) {
      not_na <- rowSums(is.na(r)) == 0
      evalues <- do.call(e_func, c(list(r = r[not_na, ], m = m), options))
    } else {
      not_na <- !is.na(r)
      evalues <- do.call(e_func, c(list(r = r[not_na], m = m), options))
    }
    e[not_na] <- evalues$e
    evalues$e <- e
    c(evalues, list(na = which(!not_na), h = 1))
  } else {
    evalues <- vector("list", h)
    f <- rep(seq_len(h), ceiling(n / h))[seq_len(n)]
    if (is.matrix(r)) {
      r_split <- lapply(split(x = data.frame(r), f), function(x) unname(as.matrix(x)))
    } else {
      r_split <- unname(split(x = r, f))
    }
    for (j in seq_len(h)) {
      tmp <- e_rank_histogram(
        r = r_split[[j]],
        h = 1,
        m = m,
        options = options,
        strategy = strategy
      )
      tmp[[length(tmp)]] <- NULL
      evalues[[j]] <- tmp
    }
    list(evalues_h = evalues, h = h)
  }
}

#' Test discrete uniform distribution against empirical distribution
#'
#' Tests the null hypothesis that \code{r} are generated from a discrete
#' uniform distribution, against the alternative defined by their empirical
#' distribution.
#'
#' @param r observations.
#' @param m size of discrete uniform distribution.
#' @param n0 minimum number of observations for starting. All e-values until
#'     \code{n0} (included) are equal to 1.
#' @param ... currently not used, only for calling the function with
#'     \code{do.call}.
#'
#' @return
#' List containing the e-values.
#'
#' @export
#'
#' @examples
#' r <- simulate_pit(1000, "rank_histogram", K = 20, bias = 0.2, dispersion = 0)
#' e <- empirical_e(r = r$r, m = 20)
#' max(cumprod(e$e))
#'
#' @name empirical_e
empirical_e <- function(r, m, n0 = 10, ...) {
  sequential_ranks(r = r, m = m, n0 = n0)
}

#' @rdname empirical_e
#' @export
empirical_e_agg <- function(r, m, avg = TRUE, n0 = 10, ...) {
  e <- sequential_ranks_agg(r = r, m = m, n0 = n0)
  if (avg) {
    e$e <- rowMeans(e$e)
  }
  return(e)
}

#' Test discrete uniform distribution against betabinomial
#'
#' Tests the null hypothesis that \code{r} are generated from a discrete
#' uniform distribution, against a betabinomial distribution (estimated by
#' maximum likelihood.)
#'
#' @param r vector or matrix of observations.
#' @param m size of discrete uniform distribution.
#' @param n0 minimum number of observations for starting. All e-values until
#'     \code{n0} (included) are equal to 1.
#' @param tol tolerance for likelihood maximization. The maximization algorithm
#'     is stopped as soon as the difference between subsequent parameter
#'     estimates is smaller than \code{tol}.
#' @param max_it maximum number of iterations for likelihood maximization.
#' @param ... currently not used, only for calling the function with
#'     \code{do.call}.
#' @param avg logical specify whether e-values should be averaged (across rows)
#'     when \code{r} is a matrix.
#'
#' @details
#' If \code{n0} is too small, the estimated parameters may diverge for small
#' samples.
#'
#' If r is a matrix, estimation is performed separately for each row - the data
#' is aggregated across the columns. In this case, \code{n0} is the number of
#' rows of observations that are used for starting. This is done using the
#' \code{betabinom_e_agg()} function.
#'
#' @return
#' A list containing the e-values (\code{e}) and the estimated parameters
#' (matrix, \code{pars}).
#'
#' @export
#'
#' @examples
#' n <- 360
#'
#' # vector
#' r <- simulate_pit(n, "rank_histogram", K = 20, bias = 0.2, dispersion = 0)
#' e <- betabinom_e(r = r$r, m = 20)
#' max(cumprod(e$e))
#'
#' # matrix
#' d <- 3
#' r <- simulate_pit(n*d, "rank_histogram", K = 20, bias = 0.2, dispersion = 0)$r
#' r <- matrix(r, ncol = d)
#' e <- betabinom_e_agg(r = r, m = 20)
#'
#' @name betabinom
betabinom_e <- function(r, m, n0 = 20, tol = 1e-7, max_it = 20, ...) {
  betabinom_e_cpp(r = r, N = m, tol = tol, max_it = max_it, n0 = n0)
}

#' @rdname betabinom
#' @export
betabinom_e_agg <- function(r, m, avg = TRUE, n0 = 20, tol = 1e-7, max_it = 20, ...) {
  e <- betabinom_e_agg_cpp(r = r, N = m, tol = tol, max_it = max_it, n0 = n0)
  if (avg) {
    e$e <- rowMeans(e$e)
  }
  return(e)
}
