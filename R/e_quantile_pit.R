#' E-values for calibration of quantile PIT
#'
#' Tests if quantile forecasts are probabilistically calibrated.
#'
#' @param zu,zl quantile PIT (can be computed with \code{\link{quantile_pit}}).
#' @param h forecast lag (positive integer). For example, daily forecasts for
#'     the next day have lag 1, daily forecasts for an event two days ahead have
#'     lag 2.
#' @param strategy strategy for evaluating calibration. Available are
#'     \code{"grenander"} for the Grenander estimator, and \code{"bernstein"}
#'     for e-values based on Bernstein polynomials.
#' @param options options for the given strategy (see \code{\link{grenander_e}}
#'     and \code{\link{bernstein_e}} for the available parameters).
#' @param check check for correct format of input parameters.
#'
#' @details
#' For continuously distributed observations (when \code{zu,zl} may only
#' attain finitely many different values), it is recommended to use
#' \code{strategy="grenander"}.
#'
#' @return
#' If \code{h} equals 1: Returns a list containing the vector of e-values
#' (\code{e}), separate e-values for the upper and lower quantile PIT
#' (\code{eu}, \code{el}), and the forecast lag \code{h}.
#'
#' If \code{h} is greater than 1: the list contains for each \code{j=1,2,...,h}
#' a list with the e-values for all observations with indices
#' \code{j,j+h,j+2h,...}.
#'
#' @seealso
#' \code{\link{quantile_pit}} for computing the quantile PIT from forecasts.
#'
#' \code{\link{grenander_e}}, \code{\link{bernstein_e}} for related tests of
#' stochastic dominance.
#'
#' @export
#'
#' @examples
#' n <- 360
#' sim <- simulate_pit(n, "quantile_pit", K = 19, bias = 0.2, dispersion = 0)
#'
#' # If z is a lag 1 forecast:
#' e <- e_quantile_pit(sim$zu, sim$zl, h = 1, strategy = "grenander")
#' evalue_merge(e)
#' prod(e$e)
#' max(cumprod(e$e))
#'
#' # Lag 2:
#' e <- e_quantile_pit(sim$zu, sim$zl, h = 2, strategy = "grenander")
#' str(e)
#' evalue_merge(e)
e_quantile_pit <- function(
    zu,
    zl,
    h,
    strategy = "grenander",
    options = list(),
    check = FALSE
  ) {
  if (check) {
    check_z(zu)
    check_z(zl)
    check_h(h)
    check_strategy(strategy, "quantile_pit")
  }
  e_func <- get(paste0(strategy, "_e"))
  n <- length(zu)
  if (h == 1) {
    zu <- 1 - zu
    not_na_u <- !is.na(zu)
    evalues_u <- rep(1, n)
    evalues_u[not_na_u] <- do.call(e_func, c(list(z = zu[not_na_u]), options))
    eu <- evalues_u$e

    not_na_l <- !is.na(zl)
    evalues_l <- rep(1, n)
    evalues_l[not_na_l] <- do.call(e_func, c(list(z = zl[not_na_l]), options))
    el <- evalues_l$e

    e <- 0.5 * (eu + el)
    list(e = e, eu = eu, el = el, na = which(na_u | na_l),  h = 1)
  } else {
    evalues <- vector("list", h)
    f <- seq_along(zu) %% h
    zu_split <- unname(split(x = zu, f))
    zl_split <- unname(split(x = zl, f))
    for (j in seq_len(h)) {
      tmp <- e_quantile_pit(
        zu = zu_split[[j]],
        zl = zl_split[[j]],
        h = 1,
        strategy = strategy,
        optiions = options
      )
      tmp[[length(tmp)]] <- NULL
      evalues[[j]] <- tmp
    }
    list(evalues_h = evalues, h = h)
  }
}

#' Test stochastic dominance with Grenander estimator.
#'
#' Tests the null hypothesis that \code{z} are generated from a distribution
#' stochastically greater than the uniform distribution.
#'
#' @param z numbers in (0, 1] (exact zeros are not allowed).
#' @param n0 minimum number of observations for starting. All e-values until
#'     \code{n0} (included) are equal to 1, and the e-value is computed
#'     starting from the \code{n0+1}th observation.
#' @param ... currently not used, only for calling the function with
#'     \code{do.call}.
#'
#' @details
#' For the reverse hypothesis (smaller than uniform distribution), replace
#' \code{z} by \code{1-z}. Apply \code{\link{evalue_correct}} to prevent e-values
#' of exactly zero.
#'
#' @return
#' A list containing the e-values (\code{e}).
#'
#' @export
#'
#' @examples
#' # Test hypothesis "greater than uniform distribution in stochastic dominance"
#'
#' ## Hypothesis is violated
#' z <- rbeta(300, 0.75, 1)
#' e <- grenander_e(z = z)
#' max(cumprod(e$e))
#'
#' ## Hypothesis is correct
#' z <- rbeta(300, 1, 0.75)
#' e <- grenander_e(z = z)
#' max(cumprod(e$e))
grenander_e <- function(z, n0 = 10, ...) {
  df <- stats::aggregate(
    data.frame(ind = seq_along(z)),
    by = list(z = z),
    FUN = list
  )
  pos_Z <- integer(length(z))
  pos_Z[unlist(df$ind)] <- rep(seq_len(nrow(df)), times = lengths(df$ind))

  e <- sequential_grenander(z = c(0, df$z), pos_Z = pos_Z)
  e[seq_len(n0)] <- 1
  list(e = e)
}


#' Test stochastic dominance with Bernstein polynomials
#'
#' Tests the null hypothesis that \code{z} are generated from a distribution
#' stochastically greater than the uniform distribution.
#'
#' @param z numbers in [0, 1].
#' @param n0 minimum number of observations for starting. All e-values until
#'     \code{n0} (included) are equal to 1, and the e-value is computed
#'     starting from the \code{n0+1}th observation.
#' @param m maximum degree of the Bernstein polynomials in the mixture. If
#'     \code{NULL}, an optimal \code{m} is estimated with the AIC, BIC or CN
#'     criterion specified.
#' @param crit the type of criterion to use for selecting the number of weights.
#'     Defaults to \code{"CN"}, other options are \code{"AIC"} and \code{"BIC"}.
#' @param settings options for \code{\link[osqp]{osqp}}.
#' @param ... currently not used, only for calling the function with
#'     \code{do.call}.
#'
#' @details
#' For the reverse hypothesis (smaller than uniform distribution), replace
#' \code{z} by \code{1-z}. Apply \code{\link{evalue_correct}} to prevent e-values
#' of exactly zero.
#'
#' @return
#' A list containing the e-values (\code{e}).
#'
#' @export
#'
#' @examples
#' # Test hypothesis "greater than uniform distribution in stochastic dominance"
#'
#' ## Hypothesis is violated
#' z <- rbeta(300, 0.75, 1)
#' e <- bernstein_e(z = z)
#' max(cumprod(e$e))
#'
#' ## Hypothesis is correct
#' z <- rbeta(300, 1, 0.75)
#' e <- bernstein_e(z = z)
#' max(cumprod(e$e))
bernstein_e <- function(
  z,
  n0 = 10,
  m = 20,
  crit = "CN",
  settings =
    list(eps_abs = 1e-5, eps_rel = 1e-5, max_iter = 4000L, verbose = FALSE),
  ...
) {
  n <- length(z)
  e <- rep(1, n)
  for (i in max(3, n0):n){
    e[i] <- umd(
      data = z[seq_len(i - 1)],
      monotone = "Decr",
      m = m,
      settings = settings
    )$dumd(z[i])
  }

  list(e = e)
}

