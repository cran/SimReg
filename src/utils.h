#include <Rcpp.h>
#include <cmath>

#ifndef INCLUDED_UTILS_H
#define INCLUDED_UTILS_H

using namespace std;

const double LOG_ODDS_TAYLOR_SERIES_THRESHOLD = 18.0;

inline double log_prob(double log_odds)
{
	return abs(log_odds) < LOG_ODDS_TAYLOR_SERIES_THRESHOLD ? log(1.0 - 1.0/(1.0 + exp(log_odds))) : (log_odds > 0.0 ? -exp(-log_odds) : log_odds - exp(log_odds));
}

#endif


