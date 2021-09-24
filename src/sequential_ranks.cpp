#include <Rcpp.h>
using namespace Rcpp;

//' Sequential computation of e-values for discrete uniform distribution
//' (empirical)
//'
//' For a sample \code{r} in \code{1,...,m}, compute e-values for testing the
//' null hypothesis that \code{r} follows the uniform distribution on
//' \code{1,...,m}.
//'
//' @param r observations.
//' @param m size of discrete uniform distribution PLUS ONE (positive integer,
//'     greater or equal to \code{max(r)+1}).
//' @param n0 minimum number of observations for starting. All e-values until
//'     \code{n0} (included) are equal to 1.
//' @param ... currently not used, only for calling the function with
//'     \code{do.call}.
//'
//' @return
//' List containing the of e-values.
//'
//' @name sequential_ranks
//' @keywords internal
//[[Rcpp::export]]
List sequential_ranks(IntegerVector r, int m, int n0) {
  int n = r.size();
  NumericVector r1(n);
  for (int i = 0; i < n; i++) r1[i] = r[i] - 1;
  NumericVector e (n);
  NumericVector w (m + 1);
  double denom = 0;
  double m1 = m + 1.0;

  for (int i = 0; i < n0; i++) {
    e[i] = 1.0;
    w[r1[i]] += 1.0;
    denom += 1.0;
  }

  if (min(w) < 0.5) {
    for (int j = 0; j <= m; j++) {
      w[j] += 1.0;
      denom++;
    }
  }

  for (int i = n0; i < n; i++) {
    e[i] = m1 * w[r1[i]] / denom;
    w[r1[i]] += 1.0;
    denom += 1.0;
  }

  return List::create(_["e"] = e);
}
