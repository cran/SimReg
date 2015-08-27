#include <Rcpp.h>
#include <cmath>
#include "HPO.h"
#include "Utils.h"

#ifndef TERM_LIST_H
#define TERM_LIST_H

using namespace Rcpp;
using namespace std;

struct term_list
{
	term_list(Rcpp::IntegerVector, Rcpp::IntegerVector, int);
	IntegerVector case_ids;
	IntegerVector term_ids;
	int num_cases;
};

#endif

