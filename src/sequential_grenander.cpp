#include <Rcpp.h>
using namespace Rcpp;

//' Sequential computation and evaluation of the Grenander estimator
//'
//' For observations \code{Z} and all indices \code{i=2,...,length(Z)}, the
//' Grenander estimator is computed with data \code{Z[1],...Z[i-1]} and
//' evaluated at \code{Z[i]}.
//'
//' @param z sorted distinct values of the data \code{Z}.
//' @param pos_Z integer vector giving for each index \code{i} the position of
//'    of the i-th data point \code{Z[i]} in \code{z}. Indexing starts at 1.
//'
//' @return
//' Vector of the density evaluated (same length as \code{pos_Z}, the first
//' entry is set to 1).
//'
//' @name sequential_grenander
//' @keywords internal
//'
// [[Rcpp::export]]
NumericVector sequential_grenander(NumericVector z, IntegerVector pos_Z) {
  // z: sorted, distinct values of Z, with 0 at beginning
  // pos_Z: position of Z[i] in z (indexing starts at 1)

  int m = z.size() - 1;
  int n = pos_Z.size();
  NumericVector N (m + 1);
  IntegerVector PP = rep(-1, m + 1);
  NumericVector WW (m + 1);
  NumericVector MM (m + 1);
  NumericVector e (n);
  e[0] = 1.0;
  N[0] = 0.0;
  MM[0] = R_PosInf;
  WW[0] = 0.0;

  // First iteration
  int d = 1;
  int j0 = pos_Z[0];
  N[j0] = 1.0;
  PP[1] = j0;
  WW[1] = z[j0];
  MM[1] = 1.0 / z[j0];
  if (j0 < m) {
    d++;
    PP[2] = m;
    WW[2] = z[m] - z[j0];
    MM[2] = 0.0;
  }

  // Prepare objects for loop
  int b0;
  int a0;
  int s0;
  int dz;
  double z0;
  IntegerVector remPP;
  NumericVector remWW;
  NumericVector remMM;

  for (int k = 1; k < n; k++) {
    // Update vector
    j0 = pos_Z[k];
    N[j0] = N[j0] + 1.0;

    // Find partition to which j0 belongs
    s0 = 1;
    while (j0 > PP[s0]) s0++;
    e[k] = MM[s0] / k; // store e-value
    a0 = PP[s0 - 1] + 1;
    b0 = PP[s0];
    dz = d;

    // Copy tail of vector, if necessary
    if (b0 < m) {
      remPP = PP[Range(s0 + 1, dz)];
      remWW = WW[Range(s0 + 1, dz)];
      remMM = MM[Range(s0 + 1, dz)];
    }

    // Update value on new partition
    d = s0;
    PP[s0] = j0;
    WW[s0] = z[j0] - z[a0 - 1];
    MM[s0] = sum(N[Range(a0, j0)]) / WW[s0];

    // Pooling
    while (MM[d-1] <= MM[d]) {
      d--;
      MM[d] = WW[d] * MM[d] + WW[d + 1] * MM[d + 1];
      WW[d] = WW[d] + WW[d + 1];
      MM[d] = MM[d] / WW[d];
      PP[d] = PP[d + 1];
    }

    // Add new partitions, pool
    if (j0 < b0) {
      z0 = z[j0];
      for (int i = j0 + 1; i <= b0; i++) {
        if (N[i] == 0 && i < b0) continue;
        d++;
        PP[d] = i;
        WW[d] = z[i] - z0;
        MM[d] = N[i] / WW[d];
        while (MM[d - 1] <= MM[d]) {
          d--;
          MM[d] = WW[d] * MM[d] + WW[d + 1] * MM[d + 1];
          WW[d] = WW[d] + WW[d + 1];
          MM[d] = MM[d]/WW[d];
          PP[d] = PP[d + 1];
        }
        z0 = z[i];
      }
    }

    // Copy (if necessary)
    if (b0 < m) {
      int l0 = dz - s0;
      IntegerVector ss = Range(d + 1, d + l0);
      PP[ss] = remPP;
      MM[ss] = remMM;
      WW[ss] = remWW;
      d = d + l0;
    }

    // Check for user interruption in R
    Rcpp::checkUserInterrupt();
  }

  return e;
}
