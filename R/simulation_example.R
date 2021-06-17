#' Simulation example for testing the functions
#'
#' Generates observations \code{y} following a standard normal distribution,
#' and the probabilistic forecast is normal with mean \code{0+bias} and
#' variance \code{1+(dispersion error)}.
#'
#' @param n number of simulations.
#' @param type \code{"pit"} for the usual PIT, \code{"quantile_pit"} for the
#'     quantile PIT, \code{"rank_histogram"} for the ensemble rank histograms.
#' @param bias the bias parameter (single number).
#' @param dispersion the dispersion error (single number greater than -1).
#' @param K number of equispaced quantiles for \code{type="quantile_pit"}, or
#'     ensemble size for \code{type="rank_histogram"}. Has no effect for the
#'     usual PIT.
#' @param seed seed for random number generation (can be omitted).
#'
#' @details
#' For the quantile PIT, the \code{K} equispaced quantiles are extracted from
#' the forecast distribution (for example, the 0.1,0.2,...0.9 quantiles for
#' \code{K=9}). For the rank histograms, \code{K} ensemble forecasts are
#' generated randomly from the forecast distribution.
#'
#' @return
#' A list containing the pit \code{z} for the \code{type="pit"}, the upper and
#' lower quantile pit \code{zu,zl} for \code{type="quantile_pit"} or the
#' ranks \code{r} for \code{type="rank_histogram"}.
#'
#' @export
#'
#' @importFrom stats rnorm
#' @importFrom stats pnorm
#' @importFrom stats qnorm
simulate_pit <- function(n, type, bias, dispersion, K = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  y <- stats::rnorm(n)
  if (identical(type, "pit")) {
    list(z = stats::pnorm(y, bias, sqrt(1 + dispersion)))
  } else if (identical(type, "quantile_pit")) {
    alpha <- 1 / (K + 1)
    qseq <- seq(alpha, 1 - alpha, alpha)
    quantiles <- matrix(
      stats::qnorm(qseq, bias, sqrt(1 + dispersion)),
      nrow = n,
      ncol = K,
      byrow = TRUE
    )
    quantile_pit(y, quantiles)
  } else if (identical(type, "rank_histogram")) {
    y <- rnorm(n)
    ens <- matrix(nrow = n, ncol = K, stats::rnorm(n * K))
    list(r = ensemble_rank(y, ens))
  } else {
    stop("'type' must be one of 'pit', 'quantile_pit', or 'rank_histogram'")
  }
}
