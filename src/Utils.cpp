#include "Utils.h"

using namespace Rcpp;
using namespace std;

double log_likelihood_beta(
	double value,
	double alpha,
	double beta 
) {
	Rcpp::NumericVector val_vec(1);
	val_vec[0] = value;
	return Rcpp::dbeta(
		val_vec,
		alpha,
		beta,
		true
	)[0];
}

double sum_vec(
	Rcpp::NumericVector v
) {
	int n = v.length();
	double result = 0;
	for (int i = 0; i < n; i++)
		result += v[i];
	return result;
}

