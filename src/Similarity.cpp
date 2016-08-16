#include "Similarity.h"

using namespace Rcpp;
using namespace std;

Rcpp::NumericVector average_across_h(
	NumericMatrix &term_term_sim_mat,
	IntegerVector &phi,
	term_list &terms
) {
	NumericVector result(terms.num_cases);
	IntegerVector terms_per_case(terms.num_cases);

	int phi_terms = phi.length();
	for (int i = 0; i < terms.num_cases; i++) {
		terms_per_case[i] = 0;
		result[i] = 0.0;
	}

	for (int i = 0; i < terms.case_ids.length(); i++) {
		double best_phi_match = 0.0;
		for (int j = 0; j < phi_terms; j++)
			best_phi_match = (best_phi_match < term_term_sim_mat(terms.term_ids[i], phi[j])) ? term_term_sim_mat(terms.term_ids[i], phi[j]) : best_phi_match;

		result[terms.case_ids[i]] += best_phi_match;
		terms_per_case[terms.case_ids[i]]++;
	}

	for (int i = 0; i < terms.num_cases; i++)
		if (result[i] > 0) result[i] = result[i] / (double)terms_per_case[i];
	
	return result;
}

Rcpp::NumericVector average_across_phi(
	NumericMatrix &term_term_sim_mat,
	IntegerVector &phi,
	term_list &terms
) {
	NumericVector result(terms.num_cases);
	for (int case_index = 0; case_index < terms.num_cases; case_index++)
		result[case_index] = 0.0;

	int phi_terms = phi.length();
	for (int phi_term_index = 0; phi_term_index < phi_terms; phi_term_index++) {
		NumericVector best_match_in_cases(terms.num_cases);
		for (int case_index = 0; case_index < terms.num_cases; case_index++) 
			best_match_in_cases[case_index] = 0.0;
		for (int term = 0; term < terms.case_ids.length(); term++)
			best_match_in_cases[terms.case_ids[term]] = (best_match_in_cases[terms.case_ids[term]] < term_term_sim_mat(terms.term_ids[term], phi[phi_term_index])) ? term_term_sim_mat(terms.term_ids[term], phi[phi_term_index]) : best_match_in_cases[terms.case_ids[term]];
		for (int case_index = 0; case_index < terms.num_cases; case_index++)
			result[case_index] += best_match_in_cases[case_index] / (double)phi_terms;
	}
	return result;
}

NumericVector transform_one_way_sim(
	NumericVector s,
	double logit_mean,
	double log_alpha_plus_beta,
	bool interpolate
) {
	int n = s.length();

	NumericVector result(n);

	double mean = 1.0-1.0/(exp(logit_mean) + 1.0);
	double inv_var = exp(log_alpha_plus_beta);

	double alpha = inv_var * mean;
	double beta = inv_var * (1.0 - mean);
		
	if (interpolate) {
		int windows = 20;
		NumericVector f(windows+1);
		for (int i = 0; i <= windows; i++) f[i] = R::pbeta((double)i/windows, alpha, beta, 1, 0);
		result = linterpolate(f, s);
	} else {
		for (int i = 0; i < n; i++)
			result[i] = R::pbeta(s[i], alpha, beta, 1, 0);
	}

	return result;

}

NumericVector transform_each_way_sim(
	pair<NumericVector, NumericVector> &each_way_sims,
	double logit_mean_f,
	double log_alpha_plus_beta_f,
	double logit_mean_g,
	double log_alpha_plus_beta_g,
	bool interpolate
) {
	int n = each_way_sims.first.length();

	NumericVector result(n);

	double mean_f = 1.0-1.0/(exp(logit_mean_f) + 1.0);
	double mean_g = 1.0-1.0/(exp(logit_mean_g) + 1.0);
	double inv_var_f = exp(log_alpha_plus_beta_f);
	double inv_var_g = exp(log_alpha_plus_beta_g);

	double alpha_f = inv_var_f * mean_f;
	double beta_f = inv_var_f * (1.0 - mean_f);
	double alpha_g = inv_var_g * mean_g;
	double beta_g = inv_var_g * (1.0 - mean_g);
		
	if (interpolate) {
		int windows = 20;
		NumericVector f(windows+1);
		NumericVector g(windows+1);
		for (int i = 0; i <= windows; i++) {
			f[i] = R::pbeta((double)i/windows, alpha_f, beta_f, 1, 0);
			g[i] = R::pbeta((double)i/windows, alpha_g, beta_g, 1, 0);
		}
		result = linterpolate(f, each_way_sims.first) * linterpolate(g, each_way_sims.second);
	} else {
		for (int i = 0; i < n; i++) {
			result[i] = R::pbeta(each_way_sims.first[i], alpha_f, beta_f, 1, 0) * R::pbeta(each_way_sims.second[i], alpha_g, beta_g, 1, 0);
		}
	}

	return result;
}

pair<NumericVector, NumericVector> get_each_way_sim(
	LogicalMatrix &row_is_column_anc,
	NumericMatrix &term_term_sim_mat,
	IntegerVector &phi,
	term_list &terms
) {
	IntegerVector phi_minimal = minimal(row_is_column_anc, phi);
	NumericVector phi_av = average_across_phi(term_term_sim_mat, phi_minimal, terms);
	NumericVector h_av = average_across_h(term_term_sim_mat, phi_minimal, terms);
	
	pair<NumericVector, NumericVector> result;

	result.first = phi_av;
	result.second = h_av;

	return result;
}
