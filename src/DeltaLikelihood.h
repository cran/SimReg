#include <Rcpp.h>
#include <cmath>
#include "HPO.h"

#ifndef DELTA_LIKELIHOOD_H
#define DELTA_LIKELIHOOD_H

using namespace Rcpp;
using namespace std;

class geom_known_prior
{
private:
public:
	typedef Rcpp::IntegerVector type;
	geom_known_prior() {}
	geom_known_prior(Rcpp::LogicalMatrix, Rcpp::NumericVector, double);
	double operator()(Rcpp::IntegerVector);
	Rcpp::LogicalMatrix row_is_column_anc;
	Rcpp::NumericVector lit_to_term_similarities;
	double rate;
	double mean_lit_sim;
};

double get_log_leaf_set_representations(
	Rcpp::LogicalMatrix row_is_column_anc,
	Rcpp::IntegerVector leaves,
	int free_terms
);

#endif
