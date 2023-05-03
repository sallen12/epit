#include <Rcpp.h>
using namespace Rcpp;

//' Probability function of betabinomial distribution
//'
//' Computes the probability function of the betabinomial distribution.
//'
//' @param x quantile (single number, not vectorized).
//' @param a first parameter of beta distribution (single number).
//' @param b first parameter of beta distribution (single number).
//' @param N number of trials (single number).
//'
//' @return
//' Probability of \code{x} under the distribution specified by the parameters.
//'
//' @name dbetabinom
//' @keywords internal
//[[Rcpp::export]]
double dbetabinom(double x, double a, double b, double N) {
  return exp(R::lchoose(N, x) + R::lbeta(x + a, N - x + b) - R::lbeta(a, b));
}

//' Wrapper for sums of digamma function
//'
//' Apply digamma function to vector and compute scalar product with another
//' vector
//'
//' @param u numeric vector (NumericVector).
//' @param v numeric vector (NumericVector, same length as \code{u})
//'
//' @return
//' Single number (double), the result.
//'
//' @name sdigammav
//' @keywords internal
double sdigammav(NumericVector u, NumericVector v) {
  // wrapper for sum(digamma(u) * v)
  int nu = u.length();
  double dgu = 0.0;
  for (int i = 0; i < nu; i++) dgu += R::digamma(u[i]) * v[i];
  return dgu;
}

//' Wrapper for sums of trigamma function
//'
//' Apply trigamma function to vector and compute scalar product with another
//' vector
//'
//' @param u numeric vector (NumericVector).
//' @param v numeric vector (NumericVector, same length as \code{u})
//'
//' @return
//' Single number (double), the result.
//'
//' @name strigammav
//' @keywords internal
double strigammav(NumericVector u, NumericVector v) {
  // wrapper for sum(trigamma(u) * v)
  int nu = u.length();
  double tgu = 0.0;
  for (int i = 0; i < nu; i++) tgu += R::trigamma(u[i]) * v[i];
  return tgu;
}

//' Maximum likelihood estimation for betabinomial distribution
//'
//' Estimates a betabinomial distribution by maximum likelihood, with
//' pre-formatted inputs for efficiency.
//'
//' @param N number of trials (int).
//' @param n total sample size (int).
//' @param z vector from 0 to N (NumericVector).
//' @param n1 vector of counts how often \code{z[i]} has been observed
//'     (NumericVector).
//' @param n2 vector of counts how often \code{N-z[i]} has been observed
//'     (NumericVector).
//' @param m1 sums of observations (double).
//' @param sums of squared observations (double).
//' @param tol tolerance for stopping (small positive number).
//' @param max_it maximum number of iterations in Newton method.
//'
//' @details
//' Estimation is performed with a Newton method and the moment matching
//' estimators as starting values. Estimated parameters below 0.001 or above 100
//' are truncated for stability.
//'
//' This function is intended for internal use only.
//'
//' @return
//' NumericVector containing the parameters of the betabinomial distribution.
//'
//' @author
//' Code is ported to Rcpp from the Rfast package (beta.mle) by Alexander Henzi.
//'
//' @name betabinom_mle_sequential
//' @keywords internal
NumericVector betabinom_mle_sequential(
    int N,
    int n,
    NumericVector z,
    NumericVector n1,
    NumericVector n2,
    double m1,
    double m2,
    double tol,
    int max_it
) {
  m1 = m1 / (1.0 * n);
  m2 = m2 / (1.0 * n);
  double down = N * m2 / m1 - N * m1 - N + m1;
  double a = (N * m1 - m2) / down;
  double b = (N - m1) * (N - m2 / m1) / down;
  if (a < 0) a = -a;
  if (b < 0) b = -b;

  double co = -n * R::digamma(N + a + b) + n * R::digamma(a + b);
  double dera = sdigammav(z + a, n1) + co - n * R::digamma(a);
  double derb = sdigammav(z + b, n2) + co - n * R::digamma(b);
  double derab = - n * R::trigamma(N + a + b) + n * R::trigamma(a + b);
  double dera2 = strigammav(z + a, n1) + derab - n * R::trigamma(a);
  double derb2 = strigammav(z + b, n2) + derab - n * R::trigamma(b);
  double det = dera2 * derb2 - pow(derab, 2);
  double anew = a - (derb2 * dera - derab * derb)/ det;
  double bnew = b - (dera2 * derb - derab * dera) / det;
  int s = 2;
  double diff = 0.0;
  if (anew < a) {
    diff += a - anew;
  } else {
    diff += anew - a;
  }
  if (bnew < b) {
    diff += b - bnew;
  } else {
    diff += bnew - b;
  }

  while ((diff > tol) && (s < max_it)) {
    s++;
    a = anew;
    b = bnew;
    co = -n * R::digamma(N + a + b) + n * R::digamma(a + b);
    dera = sdigammav(z + a, n1) + co - n * R::digamma(a);
    derb = sdigammav(z + b, n2) + co - n * R::digamma(b);
    derab = - n * R::trigamma(N + a + b) + n * R::trigamma(a + b);
    dera2 = strigammav(z + a, n1) + derab - n * R::trigamma(a);
    derb2 = strigammav(z + b, n2) + derab - n * R::trigamma(b);
    det = dera2 * derb2 - pow(derab, 2);
    anew = a - (derb2 * dera- derab * derb)/ det;
    bnew = b - (dera2 * derb - derab * dera) / det;
    diff = 0.0;
    if (anew < a) {
      diff += a - anew;
    } else {
      diff += anew - a;
    }
    if (bnew < b) {
      diff += b - bnew;
    } else {
      diff += bnew - b;
    }
  }

  if (anew < 0.001) anew = 0.001;
  if (bnew < 0.001) bnew = 0.001;
  if (anew > 100) anew = 100;
  if (bnew > 100) bnew = 100;

  NumericVector out (2);
  out[0] = anew;
  out[1] = bnew;
  return out;
}

//' Compute e-values based on betabinomial distrbution
//'
//' Sequentially estimate a betabinomial distribution and compute e-values.
//'
//' @param r numbers in \code{1,...,m}.
//' @param N number of trials (equal to \code{m} above).
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
//' \code{r[i]}) and a matrix of the parameters for each index \code{i}.
//'
//' @name betabinom_e_cpp
//' @keywords internal
//[[Rcpp::export]]
List betabinom_e_cpp(NumericVector r, int N, double tol, int max_it, int n0) {
  int n = r.size();
  NumericVector r1(n);
  for (int i = 0; i < n; i++) r1[i] = r[i] - 1;

  double m1 = 0.0;
  double m2 = 0.0;
  double a;
  double b;
  double N1 = 1.0 + N;

  NumericVector e(n, 1.0);
  NumericMatrix par(n, 2);
  NumericVector z(N + 1);
  NumericVector n1(N + 1);
  NumericVector n2(N + 1);

  for (int j = 0; j < N1; j++) {
    z[j] = j;
    n1[j] = 0.0;
    n2[j] = 0.0;
  }

  for (int i = 0; i < n0; i++) {
    e[i] = 1;
    par(i, 0) = 1.0;
    par(i, 1) = 1.0;
    n1[r1[i]] += 1.0;
    n2[N - r1[i]] += 1.0;
    m1 += r1[i];
    m2 += pow(r1[i], 2);
  }

  for (int i = n0; i < n; i++) {
    par(i, _) = betabinom_mle_sequential(N, i, z, n1, n2, m1, m2, tol, max_it);

    a = par(i, 0);
    b = par(i, 1);
    e[i] = dbetabinom(r1[i], a, b, N) * N1;
    m1 += r1[i];
    m2 += pow(r1[i], 2);
    n1[r1[i]] += 1.0;
    n2[N - r1[i]] += 1.0;
  }

  return List::create(_["e"] = e, _["pars"] = par);
}
