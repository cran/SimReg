#include <Rcpp.h>
#include <cmath>
#include "State.h"
#include "Utils.h"
#include "Update.h"

#ifndef CHAIN_H
#define CHAIN_H

using namespace Rcpp;
using namespace std;

List Chain(
	Likelihood& likelihood,
	Update& updater,
	Data& d, 
	LogicalMatrix row_is_column_anc, 
	NumericMatrix ttsm,
	State& cur_state,
	int thin,
	int iterations,
	int phi_size,
	double temperature,
	bool record_x,
	bool fix_phi,
	bool record_model_likelhoods
);

List SA_Chain(
	Likelihood& likelihood,
	Update& updater,
	Data& d, 
	LogicalMatrix row_is_column_anc, 
	NumericMatrix ttsm,
	State& cur_state,
	int thin,
	int iterations,
	int pre_anneal_its,
	bool fix_phi,
	double initial_temperature,
	int phi_size
);

#endif
