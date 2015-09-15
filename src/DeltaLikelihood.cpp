#include "DeltaLikelihood.h"

using namespace Rcpp;
using namespace std;

double get_log_leaf_set_representations(
	Rcpp::LogicalMatrix row_is_column_anc,
	Rcpp::IntegerVector leaves,
	int free_terms
) {
	int ancestors = 1;
	int num_leaves = leaves.length();
	int pot_ancestors = row_is_column_anc.nrow();
	for (int i = 0; i < pot_ancestors; i++) {
		bool any_ancs = false;
		for (int leaf = 0; leaf < num_leaves; leaf++)
			any_ancs |= row_is_column_anc(i, leaves[leaf]);
		ancestors += (int)any_ancs;
	}

	double result = 0.0;

	for (int i = 1; i <= num_leaves; i++)
		result += log((double)i);
	result += (double)(free_terms) * log((double)ancestors);
	return result;
}

geom_known_prior::geom_known_prior(
	Rcpp::LogicalMatrix in_row_is_column_anc,
	Rcpp::NumericVector in_lit_to_term_similarities,
	double in_rate
) {
	row_is_column_anc = in_row_is_column_anc;
	lit_to_term_similarities = in_lit_to_term_similarities;
	rate = in_rate;
	mean_lit_sim = 0.0;
	int n = lit_to_term_similarities.length();
	for (int i = 0; i < n; i++)
		mean_lit_sim += lit_to_term_similarities[i];
	mean_lit_sim /= (double)n;
}

double log_nCr(int n, int r) {
	double result = 0.0;
	for (int i = r + 1; i <= n; i++)
		result += log((double)i);
	for (int i = 1; i <= n-r; i++)
		result -= log((double)i);
	return result;
}

double geom_known_prior::operator()(Rcpp::IntegerVector phi) { 
	int lit_sim = 0.0;
	int N = phi.length();
	for (int i = 0; i < N; i++)
		lit_sim += lit_to_term_similarities(phi[i]);
	lit_sim = lit_sim / (double)N;

	Rcpp::LogicalVector is_leaf = hpo_leaves(row_is_column_anc, phi);

	int num_leaves = 0;
	for (int i = 0; i < N; i++) {
		num_leaves += (int)is_leaf[i];
	}
	Rcpp::IntegerVector leaves(num_leaves);
	int leaf_no = 0;
	for (int i = 0; i < N; i++) {
		if (is_leaf[i]) {
			leaves[leaf_no] = phi[i];
			leaf_no++;
		}
	}

	int num_exclusive_ancs = 0;
	for (int pot_anc = 0; pot_anc < row_is_column_anc.ncol(); pot_anc++)
		for (int leaf = 0; leaf < leaves.length(); leaf++) {
			if (row_is_column_anc(pot_anc, leaves[leaf])) {
				num_exclusive_ancs++;
				break;
			}
		}

	double total_mass = 0.0;
	for (int i = 1; i <= N; i++) {
		total_mass += exp(i * log(rate)) * exp(log_nCr(row_is_column_anc.nrow(), i));
	}

	double normalise_by = exp(num_leaves * log(rate)) / total_mass ;

	return (double)num_leaves * log(rate) + log(normalise_by) + log(lit_sim) - log(mean_lit_sim) - log((double)num_leaf_set_reps(
		num_leaves,
		num_exclusive_ancs,
		phi.length()
	));
}
