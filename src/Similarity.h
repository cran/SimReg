#include <Rcpp.h>
#include <cmath>
#include "TermList.h"
#include "Utils.h"

#ifndef SIMILARITY_H
#define SIMILARITY_H

using namespace Rcpp;
using namespace std;

Rcpp::NumericVector average_across_h(
	NumericMatrix term_term_sim_mat,
	IntegerVector phi,
	term_list terms
);

Rcpp::NumericVector average_across_phi(
	NumericMatrix term_term_sim_mat,
	IntegerVector phi,
	term_list terms
);

NumericVector transform_each_way_sim(
	pair<NumericVector, NumericVector> each_way_sims,
	double logit_mean_f,
	double log_alpha_plus_beta_f,
	double logit_mean_g,
	double log_alpha_plus_beta_g
);

pair<NumericVector, NumericVector> get_each_way_sim(
	LogicalMatrix row_is_column_anc,
	NumericMatrix term_term_sim_mat,
	IntegerVector phi,
	term_list terms
);

#endif
