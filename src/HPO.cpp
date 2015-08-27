#include "HPO.h"

using namespace Rcpp;
using namespace std;

Rcpp::LogicalVector hpo_leaves(
	Rcpp::LogicalMatrix row_is_column_anc, 
	Rcpp::IntegerVector terms
) {
	int N = terms.length();
	Rcpp::LogicalVector is_leaf(N);
	for (int candidate_leaf = 0; candidate_leaf < N; candidate_leaf++) {
		is_leaf[candidate_leaf] = true;
		for (int possible_match = 0; possible_match < candidate_leaf; possible_match++)
			if (terms[possible_match] == terms[candidate_leaf])
				is_leaf[candidate_leaf] = false;
		for (int possible_descendant = 0; possible_descendant < N; possible_descendant++)
			if (
				(possible_descendant != candidate_leaf) 
					&& 
				row_is_column_anc(terms[candidate_leaf], terms[possible_descendant]) 
			)
				is_leaf[candidate_leaf] = false;
	}

	return is_leaf;
}

Rcpp::IntegerVector minimal(
	Rcpp::LogicalMatrix row_is_column_anc,
	Rcpp::IntegerVector terms
) {
	Rcpp::LogicalVector is_leaf = hpo_leaves(row_is_column_anc, terms);
	int num_leaves = 0;
	int num_terms = is_leaf.length();
	for (int i = 0; i < num_terms; i++)
		if (is_leaf[i]) 
			num_leaves++;

	Rcpp::IntegerVector result(num_leaves);
	int cursor = 0;
	for (int i = 0; i < num_terms; i++) {
		if (is_leaf[i]) {
			result[cursor] = terms[i];
			cursor++;
		}
	}
	return result;
}

Rcpp::LogicalMatrix hpo_leaves_matrix(
	Rcpp::LogicalMatrix row_is_column_anc,
	Rcpp::IntegerMatrix terms_matrix
) {
	int each_row = terms_matrix.ncol();
	int rows = terms_matrix.nrow();
	Rcpp::LogicalMatrix result(rows, each_row);
	for (int row_num = 0; row_num < rows; row_num++) {
		Rcpp::IntegerVector terms = terms_matrix(row_num, _);
		Rcpp::LogicalVector is_leaf = hpo_leaves(row_is_column_anc, terms);
		for (int col = 0; col < each_row; col++)
			result(row_num, col) = is_leaf[col];
	}
	return result;
}

bool next_capped_combo(IntegerVector item, int n, int cap) {
	int used = 0;
	for (size_t i = 0; i < n; i++)
		used += item[i];
		
	int focus = 0;
    while (focus < n) {
        if (used < cap) {
            item[focus]++;
			return true;
		} else {
			used -= (item[focus]-1);
			item[focus] = 1;
		}
		focus++;
    }
    return false;
}

int int_pow(int base, int pow) {
	int result = 1;
	for (int i = 0; i < pow; i++)
		result *= base;
	return result;
}

int ncr(int n, int r) {
	int result = 1;
	for (int i = r + 1; i <= n; i++)
		result *= i;
	for (int i = 1; i <= n-r; i++)
		result /= i;
	return result;
}

int factorial(int n) {
	int result = 1;
	for (int i = 1; i <= n; i++)
		result *= i;
	return result;
}

int num_orderings(IntegerVector items_of_each_kind) {
	int total_items = 0;
	for (int i = 0; i < items_of_each_kind.length(); i++)
		total_items += items_of_each_kind[i];

	int result = factorial(total_items);
	for (int i = 0; i < items_of_each_kind.length(); i++)
		result /= factorial(items_of_each_kind[i]);
	
	return result;
}

int num_leaf_set_reps(int leaves, int num_exc_ancs, int vector_length) {
	int result = 0;

	IntegerVector anc_combs(vector_length - leaves + 1);
	IntegerVector leaf_placement_combos(vector_length - leaves + 1);
	for (int i = 0; i <= vector_length - leaves; i++) {
		anc_combs[i] = int_pow(num_exc_ancs, i);
		leaf_placement_combos[i] = ncr(vector_length, vector_length-i);
	}

	IntegerVector leaf_combo(leaves, 1);		
	do {
		int leaf_places = 0;
		for (int i = 0; i < leaf_combo.length(); i++)
			leaf_places += leaf_combo[i];
		int free_terms = vector_length - leaf_places;

		result += num_orderings(leaf_combo) * anc_combs[free_terms] * leaf_placement_combos[free_terms];
	} 
	while (next_capped_combo(leaf_combo, leaves, vector_length));

	return result;
}


