// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// ripser_cpp_dist
NumericVector ripser_cpp_dist(const NumericVector& dist_r, int dim, float thresh, int p);
RcppExport SEXP _TDAstats_ripser_cpp_dist(SEXP dist_rSEXP, SEXP dimSEXP, SEXP threshSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type dist_r(dist_rSEXP);
    Rcpp::traits::input_parameter< int >::type dim(dimSEXP);
    Rcpp::traits::input_parameter< float >::type thresh(threshSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(ripser_cpp_dist(dist_r, dim, thresh, p));
    return rcpp_result_gen;
END_RCPP
}
// ripser_cpp
NumericVector ripser_cpp(const NumericMatrix& input_points, int dim, float thresh, int p, int format);
RcppExport SEXP _TDAstats_ripser_cpp(SEXP input_pointsSEXP, SEXP dimSEXP, SEXP threshSEXP, SEXP pSEXP, SEXP formatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type input_points(input_pointsSEXP);
    Rcpp::traits::input_parameter< int >::type dim(dimSEXP);
    Rcpp::traits::input_parameter< float >::type thresh(threshSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type format(formatSEXP);
    rcpp_result_gen = Rcpp::wrap(ripser_cpp(input_points, dim, thresh, p, format));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TDAstats_ripser_cpp_dist", (DL_FUNC) &_TDAstats_ripser_cpp_dist, 4},
    {"_TDAstats_ripser_cpp", (DL_FUNC) &_TDAstats_ripser_cpp, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_TDAstats(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}