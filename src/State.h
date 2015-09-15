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
	void initialise(
		Likelihood& likelihood,
		Data& d,
		LogicalMatrix row_is_column_anc,
		NumericMatrix ttsm
	);

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

	double get_gamma_lik(Likelihood lik, double temperature = 1.0);
	double get_alpha_star_lik(Likelihood lik, double temperature = 1.0);
	double get_alpha_lik(Likelihood lik, double temperature = 1.0);
	double get_log_beta_lik(Likelihood lik, double temperature = 1.0);
	double get_phi_lik(Likelihood lik, double temperature = 1.0);
	double get_logit_mean_f_lik(Likelihood lik, double temperature = 1.0);
	double get_log_alpha_plus_beta_f_lik(Likelihood lik, double temperature = 1.0);
	double get_logit_mean_g_lik(Likelihood lik, double temperature = 1.0);
	double get_log_alpha_plus_beta_g_lik(Likelihood lik, double temperature = 1.0);
	double get_y_lik(Likelihood lik, NumericMatrix ttsm, LogicalMatrix row_is_column_anc, Data& d, double temperature = 1.0, bool reparameterise = true, bool quantile_normalise = false);
	NumericVector get_x(NumericMatrix ttsm, LogicalMatrix row_is_column_anc, Data& d, bool reparameterise = true, bool quantile_normalise = false);
	pair<NumericVector, NumericVector> get_s(NumericMatrix ttsm, LogicalMatrix row_is_column_anc, Data& d, bool reparameterise = true, bool quantile_normalise = false);

	double get_total_lik(LogicalMatrix row_is_column_anc, NumericMatrix ttsm, Likelihood& lik, Data& d, double temperature = 1.0, bool reparameterise = true, bool quantile_normalise = false);
	
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

	double cur_gamma_lik;
	double cur_alpha_star_lik;
	double cur_alpha_lik;
	double cur_log_beta_lik;
	double cur_logit_mean_f_lik;
	double cur_log_alpha_plus_beta_f_lik;
	double cur_logit_mean_g_lik;
	double cur_log_alpha_plus_beta_g_lik;
	double cur_phi_lik;
	double cur_y_lik;

	pair<NumericVector, NumericVector> _s;
	NumericVector _x;

};

#endif
