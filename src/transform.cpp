#include <Rcpp.h>
#include <cmath>
#include "utils.h"

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

// [[Rcpp::export]]
NumericVector transform(
	NumericVector q,
	double shape1,
	double shape2,
	int windows=5
) {
	NumericVector f(windows+1);
	NumericVector result(q.length());
	for (int i = 0; i <= windows; i++) {
		f[i] = R::pbeta((double)i/windows, shape1, shape2, 1, 0);
	}
	result = linterpolate(f, q);
	return result;
}

IntegerMatrix sq_backwards(IntegerMatrix sq) {
	int N = sq.nrow();
	IntegerMatrix result(N,N);
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			result(i,j) = sq(N-i-1,N-j-1);
	return result;
}

// [[Rcpp::export]]
IntegerMatrix sumgrid(
	NumericVector s_phi,
	NumericVector s_x,
	int breaks
) {
	IntegerMatrix count(breaks, breaks);
	int N = s_phi.length();
	double scaler = (double)(breaks-1);
	for (int i = 0; i < N; i++) {
		count((int)(s_phi[i] * scaler), (int)(s_x[i] * scaler)) += 1;
	}

	IntegerMatrix backwards = sq_backwards(count);
	
	for (int i = 1; i < breaks; i++) {
		for (int j = 0; j < breaks; j++) {
			backwards(i, j) += backwards(i-1,j);
		}
	}
	for (int i = 1; i < breaks; i++) {
		for (int j = 0; j < breaks; j++) {
			backwards(j, i) += backwards(j,i-1);
		}
	}
	
//	cumsumgrid <- function(x) apply(apply(x, 1, cumsum), 1, cumsum)
//	cumsumgrid(tab[sq,sq])[sq,sq]

	return sq_backwards(backwards);
}


