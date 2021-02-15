// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// bg_ML
double bg_ML(int y0, int y1, NumericVector t, int n_samples, double alpha_mean, double alpha_sd, double alpha_prop_sd);
RcppExport SEXP _SimReg_bg_ML(SEXP y0SEXP, SEXP y1SEXP, SEXP tSEXP, SEXP n_samplesSEXP, SEXP alpha_meanSEXP, SEXP alpha_sdSEXP, SEXP alpha_prop_sdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type y0(y0SEXP);
    Rcpp::traits::input_parameter< int >::type y1(y1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< int >::type n_samples(n_samplesSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_mean(alpha_meanSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_sd(alpha_sdSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_prop_sd(alpha_prop_sdSEXP);
    rcpp_result_gen = Rcpp::wrap(bg_ML(y0, y1, t, n_samples, alpha_mean, alpha_sd, alpha_prop_sd));
    return rcpp_result_gen;
END_RCPP
}
// f_ML
double f_ML(NumericVector x_values, IntegerVector num_y0_phi, IntegerVector num_y1_phi, NumericVector t, double log_scale_tolerance, int min_samples, int max_samples, double min_log_ML, double alpha_mean, double alpha_sd, double log_beta_mean, double log_beta_sd, double logit_f_mean, double logit_f_sd, double log_f_a_plus_b_mean, double log_f_a_plus_b_sd, double alpha_prop_sd, double log_beta_prop_sd, double logit_f_mean_prop_sd, double log_f_a_plus_b_prop_sd);
RcppExport SEXP _SimReg_f_ML(SEXP x_valuesSEXP, SEXP num_y0_phiSEXP, SEXP num_y1_phiSEXP, SEXP tSEXP, SEXP log_scale_toleranceSEXP, SEXP min_samplesSEXP, SEXP max_samplesSEXP, SEXP min_log_MLSEXP, SEXP alpha_meanSEXP, SEXP alpha_sdSEXP, SEXP log_beta_meanSEXP, SEXP log_beta_sdSEXP, SEXP logit_f_meanSEXP, SEXP logit_f_sdSEXP, SEXP log_f_a_plus_b_meanSEXP, SEXP log_f_a_plus_b_sdSEXP, SEXP alpha_prop_sdSEXP, SEXP log_beta_prop_sdSEXP, SEXP logit_f_mean_prop_sdSEXP, SEXP log_f_a_plus_b_prop_sdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x_values(x_valuesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type num_y0_phi(num_y0_phiSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type num_y1_phi(num_y1_phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type log_scale_tolerance(log_scale_toleranceSEXP);
    Rcpp::traits::input_parameter< int >::type min_samples(min_samplesSEXP);
    Rcpp::traits::input_parameter< int >::type max_samples(max_samplesSEXP);
    Rcpp::traits::input_parameter< double >::type min_log_ML(min_log_MLSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_mean(alpha_meanSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_sd(alpha_sdSEXP);
    Rcpp::traits::input_parameter< double >::type log_beta_mean(log_beta_meanSEXP);
    Rcpp::traits::input_parameter< double >::type log_beta_sd(log_beta_sdSEXP);
    Rcpp::traits::input_parameter< double >::type logit_f_mean(logit_f_meanSEXP);
    Rcpp::traits::input_parameter< double >::type logit_f_sd(logit_f_sdSEXP);
    Rcpp::traits::input_parameter< double >::type log_f_a_plus_b_mean(log_f_a_plus_b_meanSEXP);
    Rcpp::traits::input_parameter< double >::type log_f_a_plus_b_sd(log_f_a_plus_b_sdSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_prop_sd(alpha_prop_sdSEXP);
    Rcpp::traits::input_parameter< double >::type log_beta_prop_sd(log_beta_prop_sdSEXP);
    Rcpp::traits::input_parameter< double >::type logit_f_mean_prop_sd(logit_f_mean_prop_sdSEXP);
    Rcpp::traits::input_parameter< double >::type log_f_a_plus_b_prop_sd(log_f_a_plus_b_prop_sdSEXP);
    rcpp_result_gen = Rcpp::wrap(f_ML(x_values, num_y0_phi, num_y1_phi, t, log_scale_tolerance, min_samples, max_samples, min_log_ML, alpha_mean, alpha_sd, log_beta_mean, log_beta_sd, logit_f_mean, logit_f_sd, log_f_a_plus_b_mean, log_f_a_plus_b_sd, alpha_prop_sd, log_beta_prop_sd, logit_f_mean_prop_sd, log_f_a_plus_b_prop_sd));
    return rcpp_result_gen;
END_RCPP
}
// ML
double ML(NumericVector s_phi_values, NumericVector s_x_values, IntegerVector num_y0_phi, IntegerVector num_y1_phi, NumericVector t, double log_scale_tolerance, int min_samples, int max_samples, double min_log_ML, double alpha_mean, double alpha_sd, double log_beta_mean, double log_beta_sd, double logit_f_mean, double logit_f_sd, double log_f_a_plus_b_mean, double log_f_a_plus_b_sd, double logit_g_mean, double logit_g_sd, double log_g_a_plus_b_mean, double log_g_a_plus_b_sd, double alpha_prop_sd, double log_beta_prop_sd, double logit_f_mean_prop_sd, double log_f_a_plus_b_prop_sd, double logit_g_mean_prop_sd, double log_g_a_plus_b_prop_sd);
RcppExport SEXP _SimReg_ML(SEXP s_phi_valuesSEXP, SEXP s_x_valuesSEXP, SEXP num_y0_phiSEXP, SEXP num_y1_phiSEXP, SEXP tSEXP, SEXP log_scale_toleranceSEXP, SEXP min_samplesSEXP, SEXP max_samplesSEXP, SEXP min_log_MLSEXP, SEXP alpha_meanSEXP, SEXP alpha_sdSEXP, SEXP log_beta_meanSEXP, SEXP log_beta_sdSEXP, SEXP logit_f_meanSEXP, SEXP logit_f_sdSEXP, SEXP log_f_a_plus_b_meanSEXP, SEXP log_f_a_plus_b_sdSEXP, SEXP logit_g_meanSEXP, SEXP logit_g_sdSEXP, SEXP log_g_a_plus_b_meanSEXP, SEXP log_g_a_plus_b_sdSEXP, SEXP alpha_prop_sdSEXP, SEXP log_beta_prop_sdSEXP, SEXP logit_f_mean_prop_sdSEXP, SEXP log_f_a_plus_b_prop_sdSEXP, SEXP logit_g_mean_prop_sdSEXP, SEXP log_g_a_plus_b_prop_sdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type s_phi_values(s_phi_valuesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s_x_values(s_x_valuesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type num_y0_phi(num_y0_phiSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type num_y1_phi(num_y1_phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type log_scale_tolerance(log_scale_toleranceSEXP);
    Rcpp::traits::input_parameter< int >::type min_samples(min_samplesSEXP);
    Rcpp::traits::input_parameter< int >::type max_samples(max_samplesSEXP);
    Rcpp::traits::input_parameter< double >::type min_log_ML(min_log_MLSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_mean(alpha_meanSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_sd(alpha_sdSEXP);
    Rcpp::traits::input_parameter< double >::type log_beta_mean(log_beta_meanSEXP);
    Rcpp::traits::input_parameter< double >::type log_beta_sd(log_beta_sdSEXP);
    Rcpp::traits::input_parameter< double >::type logit_f_mean(logit_f_meanSEXP);
    Rcpp::traits::input_parameter< double >::type logit_f_sd(logit_f_sdSEXP);
    Rcpp::traits::input_parameter< double >::type log_f_a_plus_b_mean(log_f_a_plus_b_meanSEXP);
    Rcpp::traits::input_parameter< double >::type log_f_a_plus_b_sd(log_f_a_plus_b_sdSEXP);
    Rcpp::traits::input_parameter< double >::type logit_g_mean(logit_g_meanSEXP);
    Rcpp::traits::input_parameter< double >::type logit_g_sd(logit_g_sdSEXP);
    Rcpp::traits::input_parameter< double >::type log_g_a_plus_b_mean(log_g_a_plus_b_meanSEXP);
    Rcpp::traits::input_parameter< double >::type log_g_a_plus_b_sd(log_g_a_plus_b_sdSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_prop_sd(alpha_prop_sdSEXP);
    Rcpp::traits::input_parameter< double >::type log_beta_prop_sd(log_beta_prop_sdSEXP);
    Rcpp::traits::input_parameter< double >::type logit_f_mean_prop_sd(logit_f_mean_prop_sdSEXP);
    Rcpp::traits::input_parameter< double >::type log_f_a_plus_b_prop_sd(log_f_a_plus_b_prop_sdSEXP);
    Rcpp::traits::input_parameter< double >::type logit_g_mean_prop_sd(logit_g_mean_prop_sdSEXP);
    Rcpp::traits::input_parameter< double >::type log_g_a_plus_b_prop_sd(log_g_a_plus_b_prop_sdSEXP);
    rcpp_result_gen = Rcpp::wrap(ML(s_phi_values, s_x_values, num_y0_phi, num_y1_phi, t, log_scale_tolerance, min_samples, max_samples, min_log_ML, alpha_mean, alpha_sd, log_beta_mean, log_beta_sd, logit_f_mean, logit_f_sd, log_f_a_plus_b_mean, log_f_a_plus_b_sd, logit_g_mean, logit_g_sd, log_g_a_plus_b_mean, log_g_a_plus_b_sd, alpha_prop_sd, log_beta_prop_sd, logit_f_mean_prop_sd, log_f_a_plus_b_prop_sd, logit_g_mean_prop_sd, log_g_a_plus_b_prop_sd));
    return rcpp_result_gen;
END_RCPP
}
// transform
NumericVector transform(NumericVector q, double shape1, double shape2, int windows);
RcppExport SEXP _SimReg_transform(SEXP qSEXP, SEXP shape1SEXP, SEXP shape2SEXP, SEXP windowsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type shape1(shape1SEXP);
    Rcpp::traits::input_parameter< double >::type shape2(shape2SEXP);
    Rcpp::traits::input_parameter< int >::type windows(windowsSEXP);
    rcpp_result_gen = Rcpp::wrap(transform(q, shape1, shape2, windows));
    return rcpp_result_gen;
END_RCPP
}
// sumgrid
IntegerMatrix sumgrid(NumericVector s_phi, NumericVector s_x, int breaks);
RcppExport SEXP _SimReg_sumgrid(SEXP s_phiSEXP, SEXP s_xSEXP, SEXP breaksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type s_phi(s_phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s_x(s_xSEXP);
    Rcpp::traits::input_parameter< int >::type breaks(breaksSEXP);
    rcpp_result_gen = Rcpp::wrap(sumgrid(s_phi, s_x, breaks));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SimReg_bg_ML", (DL_FUNC) &_SimReg_bg_ML, 7},
    {"_SimReg_f_ML", (DL_FUNC) &_SimReg_f_ML, 20},
    {"_SimReg_ML", (DL_FUNC) &_SimReg_ML, 27},
    {"_SimReg_transform", (DL_FUNC) &_SimReg_transform, 4},
    {"_SimReg_sumgrid", (DL_FUNC) &_SimReg_sumgrid, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_SimReg(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
