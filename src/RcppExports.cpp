// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// solver_wrapper
void solver_wrapper();
RcppExport SEXP _RFATODE_solver_wrapper() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    solver_wrapper();
    return R_NilValue;
END_RCPP
}
// explicit_RK
void explicit_RK();
RcppExport SEXP _RFATODE_explicit_RK() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    explicit_RK();
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RFATODE_solver_wrapper", (DL_FUNC) &_RFATODE_solver_wrapper, 0},
    {"_RFATODE_explicit_RK", (DL_FUNC) &_RFATODE_explicit_RK, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_RFATODE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
