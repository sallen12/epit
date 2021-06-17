#include <Rcpp.h>
using namespace Rcpp;

//' Maximum likelihood estimation for beta distribution
//'
//' Estimates a beta-distribution by maximum likelihood, with pre-formatted
//' inputs for efficiency.
//'
//' @param tol tolerance for stopping (small positive number).
//' @param max_it maximum number of iterations in Newton method.
//' @param sly1 sum of log-observations.
//' @param sly2 sum of logs of 1 minus the observations.
//' @param sy sum of observations.
//' @param sy2 of squared observations.
//' @param n sample size.
//'
//' @details
//' Estimation is performed with a Newton method and the moment matching
//' estimators as starting values. Estimated parameters below 0.001 or above 100
//' are truncated for stability.
//'
//' This function is intended for internal use only.
//'
//' @return
//' NumericVector containing the parameters of the beta distribution.
//'
//' @author
//' Code is ported to Rcpp from the Rfast package (beta.mle) by Alexander Henzi.
//'
//' @keywords internal
NumericVector beta_mle_sequential(
    double tol,
    int max_it,
    double sly1,
    double sly2,
    double sy,
    double sy2,
    double n
) {
  //

  sly1 = sly1 / (1.0 * n);
  sly2 = sly2 / (1.0 * n);

  double iniphi = (sy - sy2) / (sy2 - pow(sy2, 2) / n);
  double a = sy * iniphi / n;
  double b = iniphi - a;
  double phi = a + b;
  double lik1 = -n * R::lbeta(a, b) + (a - 1.0) * n * sly1 + (b - 1.0) * sly2 * n;
  double dera = sly1 - R::digamma(a) + R::digamma(phi);
  double derb = sly2 - R::digamma(b) + R::digamma(phi);
  double derab = R::trigamma(phi);
  double dera2 = -R::trigamma(a) + derab;
  double derb2 = -R::trigamma(b) + derab;
  double denom = dera2 * derb2 - pow(derab, 2);
  a = a - (derb2 * dera - derab * derb) / denom;
  b = b - (dera2 * derb - derab * dera) / denom;
  phi = a + b;
  double lik2 = -n * R::lbeta(a, b) + (a - 1.0) * n * sly1 + (b - 1.0) * sly2 * n;
  int i = 2;
  while ((lik2 - lik1 > tol) && (i < max_it)) {
    i++;
    lik1 = lik2;
    dera = sly1 - R::digamma(a) + R::digamma(phi);
    derb = sly2 - R::digamma(b) + R::digamma(phi);
    derab = R::trigamma(phi);
    dera2 = -R::trigamma(a) + derab;
    derb2 = -R::trigamma(b) + derab;
    denom = dera2 * derb2 - pow(derab, 2);
    a = a - (derb2 * dera - derab * derb) / denom;
    b = b - (dera2 * derb - derab * dera) / denom;
    phi = a + b;
    lik2 = -n * R::lbeta(a, b) + (a - 1.0) * n * sly1 + (b - 1.0) * sly2 * n;
  }

  if (a < 0.001) a = 0.001;
  if (b < 0.001) b = 0.001;
  if (a > 100) a = 100;
  if (b > 100) b = 100;

  // export results
  NumericVector pars(2);
  pars[0] = a;
  pars[1] = b;
  return pars;
}

//' Compute e-values based on beta distrbution
//'
//' Sequentially estimate a beta distribution and evaluate the density at the
//' next point in the sample.
//'
//' @param z numbers in (0,1), should not be equal to 0 or 1.
//' @param tol tolerance for stopping (small positive number).
//' @param max_it maximum number of iterations in Newton method.
//' @param n0 starting sample size. The first \code{n0} observations are used
//'     only for estimation and the e-values are not computed.
//'
//' @details
//' This function is intended for internal use only.
//'
//' @return
//' List containing the e-values (density evaluated at the given value of
//' \code{z[i]}) and a matrix of the parameters for each index \code{i}.
//'
//' @keywords internal
//[[Rcpp::export]]
List beta_e_cpp(NumericVector z, double tol, int max_it, int n0) {
  int n = z.size();

  NumericVector z2 = pow(z, 2);
  NumericVector lz = log(z);
  NumericVector l1_z = log(1 - z);

  NumericVector e(n, 1.0);
  NumericMatrix par(n, 2);

  double sly = 0.0;
  double sly2 = 0.0;
  double sy = 0.0;
  double sy2 = 0.0;

  double shape1 = 1.0;
  double shape2 = 1.0;

  for (int i = 0; i < n0; i++) {
    e[i] = 1;
    par(i, 0) = 1.0;
    par(i, 1) = 1.0;
    sly += lz[i];
    sly2 += l1_z[i];
    sy += z[i];
    sy2 += z2[i];
  }

  for (int i = n0; i < n; i++) {
    par(i, _) = beta_mle_sequential(tol, max_it, sly, sly2, sy, sy2, i);
    shape1 = par(i, 0);
    shape2 = par(i, 1);
    e[i] = R::dbeta(z[i], shape1, shape2, 0);
    sly += lz[i];
    sly2 += l1_z[i];
    sy += z[i];
    sy2 += z2[i];
  }

  return List::create(_["e"] = e, _["par"] = par);
}
