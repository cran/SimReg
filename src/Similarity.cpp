#include "Similarity.h"

using namespace Rcpp;
using namespace std;

NumericVector average_across_h(
	NumericMatrix term_term_sim_mat,
	IntegerVector phi,
	term_list terms
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
		result[i] = result[i] / (double)terms_per_case[i];
	
	return result;
}

NumericVector average_across_phi(
	NumericMatrix term_term_sim_mat,
	IntegerVector phi,
	term_list terms
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

NumericVector transform_each_way_sim(
	pair<NumericVector, NumericVector> each_way_sims,
	double logit_mean_f,
	double log_alpha_plus_beta_f,
	double logit_mean_g,
	double log_alpha_plus_beta_g,
	bool reparam
) {
	int n = each_way_sims.first.length();

	NumericVector result(n);

	double mean_f = 1.0-1.0/(exp(logit_mean_f) + 1.0);
	double mean_g = 1.0-1.0/(exp(logit_mean_g) + 1.0);
	double inv_var_f = exp(log_alpha_plus_beta_f);
	double inv_var_g = exp(log_alpha_plus_beta_g);

	double alpha_f;
	double beta_f;
	double alpha_g;
	double beta_g;

	if (reparam) {
		
		alpha_f = inv_var_f * mean_f;
		beta_f = inv_var_f * (1.0 - mean_f);
		alpha_g = inv_var_g * mean_g;
		beta_g = inv_var_g * (1.0 - mean_g);
		
	} else {
		alpha_f = exp(logit_mean_f);
		beta_f = exp(log_alpha_plus_beta_f);
		alpha_g = exp(logit_mean_g);
		beta_g = exp(log_alpha_plus_beta_g);
	}
		
	for (int i = 0; i < n; i++) {
		result[i] = R::pbeta(each_way_sims.first[i], alpha_f, beta_f, 1, 0) * R::pbeta(each_way_sims.second[i], alpha_g, beta_g, 1, 0);
	}

	return result;
}

pair<NumericVector, NumericVector> get_each_way_sim(
	LogicalMatrix row_is_column_anc,
	NumericMatrix term_term_sim_mat,
	IntegerVector phi,
	term_list terms,
	bool quantile_normalise
) {
	IntegerVector phi_minimal = minimal(row_is_column_anc, phi);
	NumericVector phi_av = average_across_phi(term_term_sim_mat, phi_minimal, terms);
	NumericVector h_av = average_across_h(term_term_sim_mat, phi_minimal, terms);
	
	pair<NumericVector, NumericVector> result;

	if (quantile_normalise) {
		NumericVector neg_phi = phi_av + 0.002 * runif(phi_av.length());
		NumericVector neg_h = h_av + 0.002 * runif(phi_av.length());

		result.first = as<NumericVector>(match(neg_phi, clone(neg_phi).sort())) / neg_phi.length();
		result.second = as<NumericVector>(match(neg_h, clone(neg_h).sort())) / neg_h.length();
	}
	else {
		result.first = phi_av;
		result.second = h_av;
	}

	return result;
}
