#' E-values for probabilistic calibration
#'
#' Tests if the probability integral transform (PIT) of forecasts is uniformly
#' distributed.
#'
#' @param z probability integral transform (numbers in [0,1]).
#' @param h forecast lag (positive integer). For example, daily forecasts for
#'     the next day have lag 1, daily forecasts for an event two days ahead have
#'     lag 2.
#' @param strategy strategy for evaluating calibration. Available are
#'     \code{"beta"} for betting against a beta distribution (estimated by
#'     maximum likelihood) or \code{"kernel"} for an e-value constructed with
#'     kernel density estimation.
#' @param options options for the given method (see \code{\link{beta_e}} and
#'     \code{\link{kernel_e}} for the available parameters).
#'
#' @return
#' If \code{h} equals 1: Returns a list containing the vector of e-values
#' (\code{e}), parameters used for the computation of the e-value (depending
#' on \code{strategy}), and the forecast lag \code{h}.
#'
#' If \code{h} is greater than 1: Instead of a vector of e-values, the list
#' contains for each \code{j=1,2,...,h} a list with the e-values for all
#' observations with indices \code{j,j+h,j+2h,...} and the corresponding
#' parameters.
#'
#' @seealso
#' \code{\link{beta_e}}, \code{\link{kernel_e}} for directly computing
#' e-values for lag 1 forecasts with the corresponding method.
#'
#' @export
e_pit <- function(z, h, strategy, options = list()) {
  e_func <- get(paste0(strategy, "_e"))
  n <- length(z)
  if (h == 1) {
    not_zero_one <- (z > 0 & z < 1)
    e <- rep(1, n)
    evalues <- do.call(e_func, c(list(z = z[not_zero_one]), options))
    e[not_zero_one] <- evalue_correct(evalues$e)
    evalues$e <- e
    c(evalues, list(zero_one = which(!not_zero_one), h = h))
  } else {
    evalues <- vector("list", h)
    f <- seq_along(z) %% h
    z_split <- unname(split(x = z, f))
    for (j in seq_len(h)) {
      evalues[[j]] <- e_pit(
        z = z_split[[j]],
        h = 1,
        strategy = strategy,
        options = options
      )
      evalues[[j]] <- evalues[[j]][length(evalues[[j]])] <- NULL
    }
    list(evalues_h = evalues, h = h)
  }
}

#' Test uniform against beta distribution
#'
#' Tests the null hypothesis that \code{z} is uniformly distributed against
#' the alternative that it follows a beta distribution with parameters not
#' equal to 1.
#'
#' @param z vector of numbers in [0,1].
#' @param n0 minimum number of observations for starting. All e-values until
#'     \code{n0} (included) are equal to 1, and the density is estimated
#'     starting from the \code{n0+1}th observation.
#' @param tol tolerance for likelihood maximization. The maximization algorithm
#'     is stopped as soon as the difference between the log-likelihood of two
#'     subsequent iterations is at most \code{tol}.
#' @param max_it maximum number of iterations for likelihood maximization.
#' @param ... currently not used, only for calling the function with
#'     \code{do.call}.
#'
#' @return
#' A list containing the vector of e-values (\code{e}) and a matrix of the
#' estimated parameters at all time points (\code{pars}).
#'
#' @export
beta_e <- function(
  z,
  n0 = 10,
  tol = 1e-6,
  max_it = 20,
  ...
) {
  beta_e_cpp(z = z, n0 = n0, tol = tol, max_it = max_it)
}

#' Test for uniformity with kernel density estimator
#'
#' Tests the null hypothesis that \code{z} is uniformly distributed against
#' any other non-uniform distribution, estimated with a kernel density
#' estimator.
#'
#' @param z vector of numbers in [0,1].
#' @param n0 minimum number of observations for starting. All e-values until
#'     \code{n0} (included) are equal to 1, and the density is estimated
#'     starting from the \code{n0+1}th observation.
#' @param scalest scale estimator in plug-in bandwidth selection (see
#'     \code{\link[KernSmooth]{dpik}}).
#' @param level number of levels of functional estimation in plug-in rule
#'     (see \code{\link[KernSmooth]{dpik}}).
#' @param gridsize see \code{\link[KernSmooth]{dpik}}.
#' @param ... currently not used, only for calling the function with
#'     \code{do.call}.#'
#' @details
#' The density is estimated with
#' \code{\link[bde]{jonesCorrectionMuller94BoundaryKernel}} and bandwidth
#' selection uses \code{\link[KernSmooth]{dpik}}. In some cases, this
#' estimator is not a density (it does not integrate to one), and it is
#' therefore approximated on the grid 0,0.01,...,0.99,1 and appropriately
#' rescaled.
#'
#' @return
#' A list containing the vector of e-values (\code{e}) and a vector of the
#' bandwidths at the corresponding indices (\code{bws}).
#'
#' @importFrom KernSmooth dpik
#' @importFrom bde jonesCorrectionMuller94BoundaryKernel
#' @importFrom bde densityCache
#' @importFrom bde distribution
#'
#' @export
kernel_e <- function(
  z,
  n0 = 10,
  scalest = "minim",
  level = 2L,
  gridsize = 401L,
  ...
) {
  n <- length(z)
  e <- numeric(n)
  bws <- numeric(n)
  bws[seq_len(n0)] <- NA
  e[seq_len(n0)] <- 1
  pos <- findInterval(z, seq(0, 1, 0.01))
  for (i in (n0 + 1):n) {
    bws[i] <- KernSmooth:::dpik(
      x = z[seq_len(i - 1)],
      scalest = scalest,
      level = level,
      kernel = "epanech",
      gridsize = gridsize,
      range.x = c(0, 1)
    )
    kernel_estim <- bde::jonesCorrectionMuller94BoundaryKernel(
      dataPoints = z[seq_len(i - 1)],
      b = min(0.5, bws[i]),
      mu = 1
    )
    e[i] <- kernel_estim@densityCache[pos[i]] /
      bde::distribution(kernel_estim, x = 1, FALSE)
  }
  list(e = e, bws = bws)
}
