#------------------------------------------------------------------------------#
# Functions used in simulation examples
#------------------------------------------------------------------------------#

#' Kolmogorov-Smirnov test with uniform distribution
#' 
#' @param z usual probability transform or \code{NULL} if \code{zu,zl} are not
#'     \code{NULL}
#' @param zu,zl upper and lower quantile PIT, or \code{NULL} if \code{z} is not
#'     \code{NULL}
#' 
#' @return 
#' P-value of the test.
ks.test.unif <- function(z = NULL, zu = NULL, zl = NULL) {
  if (is.null(zu)) {
    ks.test(x = z, y = "punif", alternative = "two.sided")$p.value
  } else {
    pu <- ks.test(x = zu, y = "punif", alternative = "less")$p.value
    pl <- ks.test(x = zl, y = "punif", alternative = "greater")$p.value
    2 * min(pu, pl)
  }
}

#' Chisquare test for uniform distribution of rank histogram
#' 
#' @param r ranks
#' @param K ensemble size
#' 
#' @return 
#' P-value of the test
chisq.test.unif <- function(r, K) {
  tabr <- table(factor(r, levels = seq_len(K + 1)))
  chisq.test(x = tabr)$p.value
}

#' Apply hypothesis testing method in simulation example
#' 
#' Intended for usage with \code{do.call}
#' 
#' @param type type of simulation (usual PIT, quantile PIT, rank histogram)
#' @param method method used for testing. Either p-value based test, or the 
#'     given strategy for computing e-values
#' @param pit_sim usual PIT, list of upper and lower quantile PIT, or ranks
#' @param K ensemble size for rank histogram simulations
#' 
#' @return 
#' The p-value of the given test. In the case of e-values, stopping is applied.
apply_method <- function(
    type,
    method,
    pit_sim,
    K = NULL
  ) {
  if (type == "pit") {
    if (method == "ptest") {
      suppressWarnings(ks.test.unif(z = pit_sim$z))
    } else {
      1 / max(cumprod(e_pit(z = pit_sim$z, h = 1, strategy = method)$e))
    }
  } else if (type == "quantile_pit") {
    zu <- pit_sim$zu
    zl <- pit_sim$zl
    if (method == "ptest") {
      suppressWarnings(ks.test.unif(zu = zu, zl = zl))
    } else {
     1 /
      max(cumprod(e_quantile_pit(zu = zu, zl = zl, h = 1, strategy = method)$e))
    }
  } else if (type == "rank_histogram") {
    if (method == "ptest") {
      chisq.test.unif(r = pit_sim$r, K = K)
    } else {
      1 / max(cumprod(e_rank_histogram(r = pit_sim$r, m = K, h = 1, strategy = method)$e))
    }
  } else if (type == "stochtest") {
    if (method == "ptest") {
      ks.test(x = pit_sim$z, y = "punif", alternative = "greater")$p.value
    } else if (method == "grenander") {
      1 / max(cumprod(evalue_correct(grenander_e(z = pit_sim$z)$e)))
    } else if (method == "bernstein") {
      1 / max(cumprod(evalue_correct(bernstein_e(z = pit_sim$z)$e)))
    }
  }
}
