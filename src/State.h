#include <Rcpp.h>
#include <cmath>
#include "Likelihood.h"
#include "TermList.h"
#include "Similarity.h"

#ifndef STATE_H
#define STATE_H

using namespace Rcpp;
using namespace std;

class Data
{
private:
public:
	Data(LogicalVector, term_list, NumericVector);
	LogicalVector y;
	term_list h;
	NumericVector g;
};

class State
{
private:
public:
	void set_random(
		Likelihood& likelihood,
		Data& d,
		LogicalMatrix row_is_column_anc,
		NumericMatrix ttsm,
		int phi_terms
	);

	static State random(
		Likelihood& likelihood,
		Data& d,
		LogicalMatrix row_is_column_anc,
		NumericMatrix ttsm,
		int phi_terms
	);

	double get_gamma_lik(Likelihood lik, double temperature);
	double get_alpha_star_lik(Likelihood lik, double temperature);
	double get_alpha_lik(Likelihood lik, double temperature);
	double get_log_beta_lik(Likelihood lik, double temperature);
	double get_phi_lik(Likelihood lik, double temperature);
	double get_logit_mean_f_lik(Likelihood lik, double temperature);
	double get_log_alpha_plus_beta_f_lik(Likelihood lik, double temperature);
	double get_logit_mean_g_lik(Likelihood lik, double temperature);
	double get_log_alpha_plus_beta_g_lik(Likelihood lik, double temperature);
	double get_y_lik(Likelihood lik, NumericMatrix ttsm, LogicalMatrix row_is_column_anc, Data& d, double temperature, bool reparameterise, bool quantile_normalise);
	NumericVector get_x(NumericMatrix ttsm, LogicalMatrix row_is_column_anc, Data& d, bool reparameterise, bool quantile_normalise);

	double get_total_lik(LogicalMatrix row_is_column_anc, NumericMatrix ttsm, Likelihood& lik, Data& d, double temperature, bool reparameterise, bool quantile_normalise);
	
	bool gamma;
	double alpha_star;
	double alpha;
	double log_beta;
	IntegerVector phi;
	double logit_mean_f;
	double log_alpha_plus_beta_f;
	double logit_mean_g;
	double log_alpha_plus_beta_g;
	double inflexion;
};

#endif
