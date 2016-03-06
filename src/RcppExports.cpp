#include "RcppExports.h"

using namespace Rcpp;
using namespace std;

RcppExport SEXP leaf_matrix(
	SEXP R_row_is_column_anc,
	SEXP R_terms_matrix
) {
BEGIN_RCPP
	IntegerMatrix terms_matrix(R_terms_matrix);
	LogicalMatrix row_is_column_anc(R_row_is_column_anc);

	return hpo_leaves_matrix(
		row_is_column_anc,
		terms_matrix
	);
		
END_RCPP
}

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
) {
BEGIN_RCPP
	NumericMatrix ttsm(R_ttsm);

	IntegerVector term_ids(R_term_ids);
	IntegerVector case_ids(R_case_ids);
	NumericVector g(R_g);

	term_list h(term_ids, case_ids, g.length());

	LogicalMatrix row_is_column_anc(R_row_is_column_anc);

	IntegerMatrix phi_trace(R_phi);

	int its = phi_trace.nrow(); 

	LogicalVector gamma_trace(R_gamma);
	NumericVector alpha_star_trace(R_alpha_star);
	NumericVector alpha_trace(R_alpha);
	NumericVector log_beta_trace(R_log_beta);
	NumericVector logit_mean_f_trace(R_logit_mean_f);
	NumericVector log_alpha_plus_beta_f_trace(R_log_alpha_plus_beta_f);
	NumericVector logit_mean_g_trace(R_logit_mean_g);
	NumericVector log_alpha_plus_beta_g_trace(R_log_alpha_plus_beta_g);

	int phi_size = phi_trace.ncol();

	NumericMatrix log_odds_trace(its, h.num_cases);
	
	for (int i = 0; i < its; i++) {
		IntegerVector phi(phi_size);
		for (int t = 0; t < phi_size; t++)
			phi[t] = phi_trace(i, t);

		pair<NumericVector, NumericVector> s = get_each_way_sim(
			row_is_column_anc,
			ttsm,
			phi,
			h
		);

		NumericVector x = transform_each_way_sim(
			s,
			logit_mean_f_trace[i],
			log_alpha_plus_beta_f_trace[i],
			logit_mean_g_trace[i],
			log_alpha_plus_beta_g_trace[i]
		);

		for (int case_num = 0; case_num < h.num_cases; case_num++)
			log_odds_trace(i, case_num) = gamma_trace[i] ? (alpha_trace[i] + exp(log_beta_trace[i]) * x[case_num] + g[case_num]) : (alpha_star_trace[i] + g[case_num]);

	}

	return log_odds_trace;
END_RCPP
}

RcppExport SEXP R_asym_sim_func(
	SEXP R_ttsm,
	SEXP R_row_is_column_anc,
	SEXP R_num_cases,
	SEXP R_term_ids,
	SEXP R_case_ids,
	SEXP R_phi,
	SEXP R_average_across_phi
) {
BEGIN_RCPP
	similarity_function s_fun = as<bool>(R_average_across_phi) ? average_across_phi : average_across_h;

	NumericMatrix ttsm(R_ttsm);

	IntegerVector term_ids(R_term_ids);
	IntegerVector case_ids(R_case_ids);

	term_list h(term_ids, case_ids, as<int>(R_num_cases));

	LogicalMatrix row_is_column_anc(R_row_is_column_anc);

	IntegerMatrix phi_trace(R_phi);

	int its = phi_trace.nrow(); 

	int phi_size = phi_trace.ncol();

	NumericMatrix av_phi_trace(its, h.num_cases);
	
	for (int i = 0; i < its; i++) {
		IntegerVector phi(phi_size);
		for (int t = 0; t < phi_size; t++)
			phi[t] = phi_trace(i, t);

		IntegerVector phi_minimal = minimal(row_is_column_anc, phi);

		NumericVector s = s_fun(ttsm, phi_minimal, h);

		for (int case_ind = 0; case_ind < h.num_cases; case_ind++) {
			av_phi_trace(i, case_ind) = s[case_ind];
		}
	}

	return av_phi_trace;
END_RCPP
}

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
	SEXP R_H,
	SEXP R_adapt_block_size,
	SEXP R_max_tuning_batches,
	SEXP R_repeats
) {
BEGIN_RCPP
	NumericVector gamma_prior_prob = as<NumericVector>(R_gamma_prior_prob);

	RNGScope scope;

	bool fix_phi = as<bool>(R_fix_phi);
	int its = as<int>(R_its);
	int thin = as<int>(R_thin);

	NumericMatrix ttsm(R_ttsm);

	bool record_x = as<bool>(R_record_x);
	bool record_model_likelihoods = as<bool>(R_record_model_likelihoods);

	IntegerVector term_ids(R_term_ids);
	IntegerVector case_ids(R_case_ids);
	LogicalVector y(R_y);
	NumericVector g(R_g);

	term_list h(term_ids, case_ids, y.length());

	LogicalMatrix row_is_column_anc(R_row_is_column_anc);
	NumericVector lit_term_sims(R_lit_term_sims);

	geom_known_prior phi_lik_func(
		row_is_column_anc,
		lit_term_sims,
		as<double>(R_phi_num_leaves_geometric_rate)
	);

	bool gamma = as<bool>(R_gamma);
	double alpha_star = as<double>(R_alpha_star);
	double alpha = as<double>(R_alpha);
	double log_beta = as<double>(R_log_beta);
	IntegerVector phi = as<IntegerVector>(R_phi);
	double logit_mean_f = as<double>(R_logit_mean_f);
	double log_alpha_plus_beta_f = as<double>(R_log_alpha_plus_beta_f);
	double logit_mean_g = as<double>(R_logit_mean_g);
	double log_alpha_plus_beta_g = as<double>(R_log_alpha_plus_beta_g);

	double alpha_star_mean = as<double>(R_alpha_star_mean);
	double alpha_mean = as<double>(R_alpha_mean);
	double alpha_star_sd = as<double>(R_alpha_star_sd);
	double alpha_sd = as<double>(R_alpha_sd);
	double log_beta_mean = as<double>(R_log_beta_mean);
	double log_beta_sd = as<double>(R_log_beta_sd);
	double logit_mean_f_mean = as<double>(R_logit_mean_f_mean);
	double logit_mean_f_sd = as<double>(R_logit_mean_f_sd);
	double log_alpha_plus_beta_f_mean = as<double>(R_log_alpha_plus_beta_f_mean);
	double log_alpha_plus_beta_f_sd = as<double>(R_log_alpha_plus_beta_f_sd);
	double logit_mean_g_mean = as<double>(R_logit_mean_g_mean);
	double logit_mean_g_sd = as<double>(R_logit_mean_g_sd);
	double log_alpha_plus_beta_g_mean = as<double>(R_log_alpha_plus_beta_g_mean);
	double log_alpha_plus_beta_g_sd = as<double>(R_log_alpha_plus_beta_g_sd);
	double pseudo_alpha_star_mean = as<double>(R_pseudo_alpha_star_mean);
	double pseudo_alpha_mean = as<double>(R_pseudo_alpha_mean);
	double pseudo_alpha_star_sd = as<double>(R_pseudo_alpha_star_sd);
	double pseudo_alpha_sd = as<double>(R_pseudo_alpha_sd);
	double pseudo_log_beta_mean = as<double>(R_pseudo_log_beta_mean);
	double pseudo_log_beta_sd = as<double>(R_pseudo_log_beta_sd);
	double pseudo_logit_mean_f_mean = as<double>(R_pseudo_logit_mean_f_mean);
	double pseudo_logit_mean_f_sd = as<double>(R_pseudo_logit_mean_f_sd);
	double pseudo_log_alpha_plus_beta_f_mean = as<double>(R_pseudo_log_alpha_plus_beta_f_mean);
	double pseudo_log_alpha_plus_beta_f_sd = as<double>(R_pseudo_log_alpha_plus_beta_f_sd);
	double pseudo_logit_mean_g_mean = as<double>(R_pseudo_logit_mean_g_mean);
	double pseudo_logit_mean_g_sd = as<double>(R_pseudo_logit_mean_g_sd);
	double pseudo_log_alpha_plus_beta_g_mean = as<double>(R_pseudo_log_alpha_plus_beta_g_mean);
	double pseudo_log_alpha_plus_beta_g_sd = as<double>(R_pseudo_log_alpha_plus_beta_g_sd);
	IntegerVector pseudo_phi_marginal_prior = as<IntegerVector>(R_pseudo_phi_marginal_prior);
	double alpha_star_proposal_sd = as<double>(R_alpha_star_proposal_sd);
	double alpha_proposal_sd = as<double>(R_alpha_proposal_sd);
	double log_beta_proposal_sd = as<double>(R_log_beta_proposal_sd);
	double logit_mean_f_proposal_sd = as<double>(R_logit_mean_f_proposal_sd);
	double log_alpha_plus_beta_f_proposal_sd = as<double>(R_log_alpha_plus_beta_f_proposal_sd);
	double logit_mean_g_proposal_sd = as<double>(R_logit_mean_g_proposal_sd);
	double log_alpha_plus_beta_g_proposal_sd = as<double>(R_log_alpha_plus_beta_g_proposal_sd);
	IntegerVector phi_jumps = as<IntegerVector>(R_phi_jumps);

	int phi_size = phi.length();

	Likelihood likelihood(
		gamma_prior_prob[0],
		alpha_star_mean,
		alpha_mean,
		alpha_star_sd,
		alpha_sd,
		log_beta_mean,
		log_beta_sd,
		logit_mean_f_mean,
		logit_mean_f_sd,
		log_alpha_plus_beta_f_mean,
		log_alpha_plus_beta_f_sd,
		logit_mean_g_mean,
		logit_mean_g_sd,
		log_alpha_plus_beta_g_mean,
		log_alpha_plus_beta_g_sd,
		pseudo_alpha_star_mean,
		pseudo_alpha_mean,
		pseudo_alpha_star_sd,
		pseudo_alpha_sd,
		pseudo_log_beta_mean,
		pseudo_log_beta_sd,
		pseudo_logit_mean_f_mean,
		pseudo_logit_mean_f_sd,
		pseudo_log_alpha_plus_beta_f_mean,
		pseudo_log_alpha_plus_beta_f_sd,
		pseudo_logit_mean_g_mean,
		pseudo_logit_mean_g_sd,
		pseudo_log_alpha_plus_beta_g_mean,
		pseudo_log_alpha_plus_beta_g_sd,
		phi_lik_func,
		pseudo_phi_marginal_prior,
		row_is_column_anc.ncol(),
		as<double>(R_H)
	);

	Update update;

	update.alpha_star_proposal_sd = alpha_star_proposal_sd;
	update.alpha_proposal_sd = alpha_proposal_sd;
	update.log_beta_proposal_sd = log_beta_proposal_sd;
	update.logit_mean_f_proposal_sd = logit_mean_f_proposal_sd;
	update.log_alpha_plus_beta_f_proposal_sd = log_alpha_plus_beta_f_proposal_sd;
	update.logit_mean_g_proposal_sd = logit_mean_g_proposal_sd;
	update.log_alpha_plus_beta_g_proposal_sd = log_alpha_plus_beta_g_proposal_sd;
	update.phi_jumps = phi_jumps;
	update.joint_proposal = as<bool>(R_joint_proposal);
	IntegerVector phi_jumps_each(row_is_column_anc.ncol(), 0);

	for (int i = 0; i < update.phi_jumps.length(); i++)
		phi_jumps_each[update.phi_jumps[i]]++;
	update.phi_jumps_each = phi_jumps_each;

	Data d(y, h, g);

	State cur_state;
	cur_state.gamma = gamma;
	cur_state.alpha_star = alpha_star;
	cur_state.alpha = alpha;
	cur_state.log_beta = log_beta;
	cur_state.phi = phi;
	cur_state.logit_mean_f = logit_mean_f;
	cur_state.log_alpha_plus_beta_f = log_alpha_plus_beta_f;
	cur_state.logit_mean_g = logit_mean_g;
	cur_state.log_alpha_plus_beta_g = log_alpha_plus_beta_g;

	cur_state.initialise(
		likelihood,
		d,
		row_is_column_anc,
		ttsm
	);

	NumericVector target_range = NumericVector::create(0.3, 0.7);
	int adapt_block_size = as<int>(R_adapt_block_size);
	int max_tuning_batches = as<int>(R_max_tuning_batches);

	int batch_number = 0;
	double adaption_rate = 4;

	if (max_tuning_batches > 0) {
		double batch_mpg = 0.0;
		do {
			batch_mpg = 0.0;
			for (int adapting_it = 0; adapting_it < adapt_block_size; adapting_it++) {
				update.update(
					cur_state,
					likelihood,
					1.0,
					d,
					row_is_column_anc,
					ttsm,
					fix_phi
				);
				if (cur_state.gamma) batch_mpg += (double)1/(double)adapt_block_size;
			}

			if (batch_mpg < target_range[0]) {
				likelihood.H += adaption_rate/sqrt(1.0 + (double)batch_number);
			} 
			if (batch_mpg > target_range[1]) {
				likelihood.H -= adaption_rate/sqrt(1.0 + (double)batch_number);
			}

			batch_number++;
		}
		while ((batch_mpg < target_range[0] || batch_mpg > target_range[1]) && batch_number <= max_tuning_batches);
	}

	List result = Chain(
		likelihood,
		update,
		d, 
		row_is_column_anc, 
		ttsm,
		cur_state,
		thin,
		its,
		phi_size,
		1.0,
		record_x,
		fix_phi,
		record_model_likelihoods
	);

	return result;

END_RCPP
}


