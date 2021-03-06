// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// baseflow_LyneHollick
NumericVector baseflow_LyneHollick(NumericVector Q, double k);
RcppExport SEXP _HSA_baseflow_LyneHollick(SEXP QSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Q(QSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(baseflow_LyneHollick(Q, k));
    return rcpp_result_gen;
END_RCPP
}
// baseflow_ChapmanMaxwell
NumericVector baseflow_ChapmanMaxwell(NumericVector Q, double k);
RcppExport SEXP _HSA_baseflow_ChapmanMaxwell(SEXP QSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Q(QSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(baseflow_ChapmanMaxwell(Q, k));
    return rcpp_result_gen;
END_RCPP
}
// baseflow_Boughton
NumericVector baseflow_Boughton(NumericVector Q, double k, double C);
RcppExport SEXP _HSA_baseflow_Boughton(SEXP QSEXP, SEXP kSEXP, SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Q(QSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type C(CSEXP);
    rcpp_result_gen = Rcpp::wrap(baseflow_Boughton(Q, k, C));
    return rcpp_result_gen;
END_RCPP
}
// baseflow_Eckhardt
NumericVector baseflow_Eckhardt(NumericVector Q, double a, double BFImax);
RcppExport SEXP _HSA_baseflow_Eckhardt(SEXP QSEXP, SEXP aSEXP, SEXP BFImaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Q(QSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type BFImax(BFImaxSEXP);
    rcpp_result_gen = Rcpp::wrap(baseflow_Eckhardt(Q, a, BFImax));
    return rcpp_result_gen;
END_RCPP
}
// ioh_min_pivots
IntegerVector ioh_min_pivots(NumericVector Q, int d, double k);
RcppExport SEXP _HSA_ioh_min_pivots(SEXP QSEXP, SEXP dSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Q(QSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(ioh_min_pivots(Q, d, k));
    return rcpp_result_gen;
END_RCPP
}
// rec_events
IntegerVector rec_events(IntegerVector imax, IntegerVector imin);
RcppExport SEXP _HSA_rec_events(SEXP imaxSEXP, SEXP iminSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type imax(imaxSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type imin(iminSEXP);
    rcpp_result_gen = Rcpp::wrap(rec_events(imax, imin));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_HSA_baseflow_LyneHollick", (DL_FUNC) &_HSA_baseflow_LyneHollick, 2},
    {"_HSA_baseflow_ChapmanMaxwell", (DL_FUNC) &_HSA_baseflow_ChapmanMaxwell, 2},
    {"_HSA_baseflow_Boughton", (DL_FUNC) &_HSA_baseflow_Boughton, 3},
    {"_HSA_baseflow_Eckhardt", (DL_FUNC) &_HSA_baseflow_Eckhardt, 3},
    {"_HSA_ioh_min_pivots", (DL_FUNC) &_HSA_ioh_min_pivots, 3},
    {"_HSA_rec_events", (DL_FUNC) &_HSA_rec_events, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_HSA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
