#include "Utils.h"

using namespace Rcpp;
using namespace std;

NumericVector linterpolate(
	NumericVector &scaffold,
	NumericVector &x
) {
	int n = x.length();
	int scaf_len = scaffold.length();
	int windows = scaf_len - 1;
	double window_length = 1.0/windows;
	NumericVector result(n);
	for (int i = 0; i < n; i++) {
		int low = (int)(x[i]/window_length);
		if (low == windows)
			result[i] = scaffold[windows];
		else
			result[i] = scaffold[low] + (scaffold[low+1]-scaffold[low])*(x[i]-low*window_length)/window_length;
	}
	return result;
}

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


