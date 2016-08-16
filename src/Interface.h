#include <Rcpp.h>
#include <cmath>
#include "Similarity.h"
#include "DeltaLikelihood.h"
#include "Chain.h"

#ifndef INCLUDED_INTERFACE_H
#define INCLUDED_INTERFACE_H

using namespace Rcpp;
using namespace std;

typedef NumericVector (*similarity_function)(NumericMatrix&,IntegerVector&,term_list&);

RcppExport SEXP leaf_matrix(
	SEXP R_row_is_column_anc,
	SEXP R_terms_matrix
);

RcppExport SEXP R_log_odds_trace(
	SEXP R_ttsm,
	SEXP R_row_is_column_anc,
	SEXP R_term_ids,
	SEXP R_case_ids,
	SEXP R_g,
	SEXP R_phi,
	SEXP R_gamma,
	SEXP R_alpha_star,
	SEXP R_alpha,
	SEXP R_log_beta,
	SEXP R_logit_mean_f,
	SEXP R_log_alpha_plus_beta_f,
	SEXP R_logit_mean_g,
	SEXP R_log_alpha_plus_beta_g
);

RcppExport SEXP R_asym_sim_func(
	SEXP R_ttsm,
	SEXP R_row_is_column_anc,
	SEXP R_num_cases,
	SEXP R_term_ids,
	SEXP R_case_ids,
	SEXP R_phi,
	SEXP R_average_across_phi
);

RcppExport SEXP R_sim_reg(
	SEXP R_its,
	SEXP R_thin,
	SEXP R_record_x,
	SEXP R_record_model_likelihoods,
	SEXP R_ttsm,
	SEXP R_term_ids,
	SEXP R_case_ids,
	SEXP R_y,
	SEXP R_g,

	SEXP R_gamma,
	SEXP R_alpha_star,
	SEXP R_alpha,
	SEXP R_log_beta,
	SEXP R_phi,
	SEXP R_logit_mean_f,
	SEXP R_log_alpha_plus_beta_f,
	SEXP R_logit_mean_g,
	SEXP R_log_alpha_plus_beta_g,

	SEXP R_gamma_prior_prob,
	SEXP R_alpha_star_mean,
	SEXP R_alpha_mean,
	SEXP R_alpha_star_sd,
	SEXP R_alpha_sd,
	SEXP R_log_beta_mean,
	SEXP R_log_beta_sd,
	SEXP R_logit_mean_f_mean,
	SEXP R_logit_mean_f_sd,
	SEXP R_log_alpha_plus_beta_f_mean,
	SEXP R_log_alpha_plus_beta_f_sd,
	SEXP R_logit_mean_g_mean,
	SEXP R_logit_mean_g_sd,
	SEXP R_log_alpha_plus_beta_g_mean,
	SEXP R_log_alpha_plus_beta_g_sd,
	SEXP R_pseudo_alpha_star_mean,
	SEXP R_pseudo_alpha_mean,
	SEXP R_pseudo_alpha_star_sd,
	SEXP R_pseudo_alpha_sd,
	SEXP R_pseudo_log_beta_mean,
	SEXP R_pseudo_log_beta_sd,
	SEXP R_pseudo_logit_mean_f_mean,
	SEXP R_pseudo_logit_mean_f_sd,
	SEXP R_pseudo_log_alpha_plus_beta_f_mean,
	SEXP R_pseudo_log_alpha_plus_beta_f_sd,
	SEXP R_pseudo_logit_mean_g_mean,
	SEXP R_pseudo_logit_mean_g_sd,
	SEXP R_pseudo_log_alpha_plus_beta_g_mean,
	SEXP R_pseudo_log_alpha_plus_beta_g_sd,
	SEXP R_pseudo_phi_marginal_prior,
	SEXP R_alpha_star_proposal_sd,
	SEXP R_alpha_proposal_sd,
	SEXP R_log_beta_proposal_sd,
	SEXP R_logit_mean_f_proposal_sd,
	SEXP R_log_alpha_plus_beta_f_proposal_sd,
	SEXP R_logit_mean_g_proposal_sd,
	SEXP R_log_alpha_plus_beta_g_proposal_sd,
	SEXP R_phi_jumps,

	SEXP R_lit_term_sims,
	SEXP R_row_is_column_anc,
	SEXP R_phi_num_leaves_geometric_rate,
	SEXP R_fix_phi,
	SEXP R_joint_proposal,
	SEXP R_annealing
);

#endif
