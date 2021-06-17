#' E-values for rank histograms
#'
#' Tests null hypothesis that rank histogram follows a uniform distribution.
#'
#' @param r observed ranks (numbers in 1,...,m+1).
#' @param h forecast lag.
#' @param m size of the ensemble (positive integer).
#' @param strategy strategy for evaluating calibration. Available are
#'     \code{"empirical"} for the empirical distribution, and \code{"betabinom"}
#'     for the beta-binomial distribution (estimated with maximum likelihood).
#' @param options options for the given strategy (see \code{\link{betabinom_e}}
#'     and \code{\link{empirical_e}} for the available parameters).
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
e_rank_histogram <- function(r, h, m, strategy, options = list()) {
  e_func <- get(paste0(strategy, "_e"))
  n <- length(r)
  if (h == 1) {
    c(do.call(e_func, c(list(r = r, m = m + 1), options)), list(h = 1))
  } else {
    evalues <- vector("list", h)
    f <- seq_along(r) %% h
    r_split <- unname(split(x = r, f))
    for (j in seq_len(h)) {
      tmp <- e_rank_histogram(
        r = r_split[[j]],
        h = 1,
        m = m + 1,
        options = options()
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
#' @param m size of discrete uniform distribution PLUS ONE (positive integer,
#'     greater or equal to \code{max(r)+1}).
#' @param n0 minimum number of observations for starting. All e-values until
#'     \code{n0} (included) are equal to 1.
#' @param ... currently not used, only for calling the function with
#'     \code{do.call}.
#'
#' @return
#' List containing the e-values.
#'
#' @export
empirical_e <- function(r, m, n0 = 20, ...) {
  sequential_ranks(r = r, m = m + 1, n0 = n0)
}

#' Test discrete uniform distribution against betabinomial
#'
#' Tests the null hypothesis that \code{r} are generated from a discrete
#' uniform distribution, against a betabinomial distribution (estimated by
#' maximum likelihood.)
#'
#' @param r observations.
#' @param m size of discrete uniform distribution PLUS ONE (positive integer,
#'     greater or equal to \code{max(r)+1}).
#' @param n0 minimum number of observations for starting. All e-values until
#'     \code{n0} (included) are equal to 1.
#' @param tol tolerance for likelihood maximization. The maximization algorithm
#'     is stopped as soon as the difference between the log-likelihood of two
#'     subsequent iterations is at most \code{tol}.
#' @param max_it maximum number of iterations for likelihood maximization.
#' @param ... currently not used, only for calling the function with
#'     \code{do.call}.
#'
#' @details
#' If \code{n0} is too small, the estimated parameters may diverge for small
#' samples.
#'
#' @return
#' A list containing the e-values (\code{e}) and the estimated parameters
#' (matrix, \code{pars}).
#'
#' @export
betabinom_e <- function(r, m, n0 = 20, tol = 1e-7, max_it = 20, ...) {
  betabinom_e_cpp(r = r, N = m, tol = tol, max_it = max_it, n0 = n0)
}
