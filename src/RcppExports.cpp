// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// beta_e_cpp
List beta_e_cpp(NumericVector z, double tol, int max_it, int n0);
RcppExport SEXP _epit_beta_e_cpp(SEXP zSEXP, SEXP tolSEXP, SEXP max_itSEXP, SEXP n0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< int >::type n0(n0SEXP);
    rcpp_result_gen = Rcpp::wrap(beta_e_cpp(z, tol, max_it, n0));
    return rcpp_result_gen;
END_RCPP
}
// dbetabinom
double dbetabinom(double x, double a, double b, double N);
RcppExport SEXP _epit_dbetabinom(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(dbetabinom(x, a, b, N));
    return rcpp_result_gen;
END_RCPP
}
// betabinom_e_cpp
List betabinom_e_cpp(NumericVector r, int N, double tol, int max_it, int n0);
RcppExport SEXP _epit_betabinom_e_cpp(SEXP rSEXP, SEXP NSEXP, SEXP tolSEXP, SEXP max_itSEXP, SEXP n0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< int >::type n0(n0SEXP);
    rcpp_result_gen = Rcpp::wrap(betabinom_e_cpp(r, N, tol, max_it, n0));
    return rcpp_result_gen;
END_RCPP
}
// betabinom_e_agg_cpp
List betabinom_e_agg_cpp(NumericMatrix r, int N, double tol, int max_it, int n0);
RcppExport SEXP _epit_betabinom_e_agg_cpp(SEXP rSEXP, SEXP NSEXP, SEXP tolSEXP, SEXP max_itSEXP, SEXP n0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< int >::type n0(n0SEXP);
    rcpp_result_gen = Rcpp::wrap(betabinom_e_agg_cpp(r, N, tol, max_it, n0));
    return rcpp_result_gen;
END_RCPP
}
// sequential_grenander
NumericVector sequential_grenander(NumericVector z, IntegerVector pos_Z);
RcppExport SEXP _epit_sequential_grenander(SEXP zSEXP, SEXP pos_ZSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type pos_Z(pos_ZSEXP);
    rcpp_result_gen = Rcpp::wrap(sequential_grenander(z, pos_Z));
    return rcpp_result_gen;
END_RCPP
}
// sequential_ranks
List sequential_ranks(IntegerVector r, int m, int n0);
RcppExport SEXP _epit_sequential_ranks(SEXP rSEXP, SEXP mSEXP, SEXP n0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type n0(n0SEXP);
    rcpp_result_gen = Rcpp::wrap(sequential_ranks(r, m, n0));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_epit_beta_e_cpp", (DL_FUNC) &_epit_beta_e_cpp, 4},
    {"_epit_dbetabinom", (DL_FUNC) &_epit_dbetabinom, 4},
    {"_epit_betabinom_e_cpp", (DL_FUNC) &_epit_betabinom_e_cpp, 5},
    {"_epit_betabinom_e_agg_cpp", (DL_FUNC) &_epit_betabinom_e_agg_cpp, 5},
    {"_epit_sequential_grenander", (DL_FUNC) &_epit_sequential_grenander, 2},
    {"_epit_sequential_ranks", (DL_FUNC) &_epit_sequential_ranks, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_epit(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
