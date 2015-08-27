#include <Rcpp.h>
#include <cmath>

#ifndef INCLUDED_UTILS_H
#define INCLUDED_UTILS_H

using namespace Rcpp;
using namespace std;

const double LOG_ODDS_TAYLOR_SERIES_THRESHOLD = 18.0;
const double MINUS_LOG_SQRT_2_PI = -0.9189385;

inline double log_ncr(int n, int r) {
	double result = 0.0;
	for (int i = r + 1; i <= n; i++)
		result += log((double)i);
	for (int i = 1; i <= n - r; i++)
		result -= log((double)i);
	return result;
}

inline int random_integer(int exc_max)
{
	return (int)(unif_rand() * (double)exc_max) % exc_max;
}

inline double sq(double x) { return x * x; }

inline double log_likelihood_normal(
	double mean,
	double sd,
	double value
){ return MINUS_LOG_SQRT_2_PI - log(sd) - 0.5 * sq((value - mean)/sd); }

double log_likelihood_beta(
	double value,
	double alpha,
	double beta
);

double sum_vec(
	Rcpp::NumericVector v
);

template<typename T>
void displayRcppVector(T vector) {
	int n = vector.length();
	for (int i = 0; i < n; i++)
		cout << vector[i] << " ";
	cout << endl;
}

inline double log_prob(double log_odds)
{
	return abs(log_odds) < LOG_ODDS_TAYLOR_SERIES_THRESHOLD ? log(1.0 - 1.0/(1.0 + exp(log_odds))) : (log_odds > 0.0 ? -exp(-log_odds) : log_odds - exp(log_odds));
}

template<typename value_type, typename vector_type>
vector_type filledRcppVector(value_type value, int n) {
	vector_type result(n);
	for (int i = 0; i < n; i++)
		result[i] = value;
	return result;
}

#endif

