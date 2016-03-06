#include <Rcpp.h>
#include <cmath>
#include "State.h"
#include "Utils.h"
#include "Similarity.h"

#ifndef UPDATE_H
#define UPDATE_H

using namespace Rcpp;
using namespace std;

class Update
{
private:
public:
	double alpha_star_proposal_sd;
	double alpha_proposal_sd;
	double log_beta_proposal_sd;
	double logit_mean_f_proposal_sd;
	double log_alpha_plus_beta_f_proposal_sd;
	double logit_mean_g_proposal_sd;
	double log_alpha_plus_beta_g_proposal_sd;
	IntegerVector phi_jumps;
	IntegerVector phi_jumps_each;
	bool joint_proposal;

	void update_gamma(State& in_state, Likelihood& likelihood, double temperature, Data& d, LogicalMatrix row_is_column_anc, NumericMatrix ttsm);
	void update_alpha_star(State& in_state, Likelihood& likelihood, double temperature, Data& d, LogicalMatrix row_is_column_anc, NumericMatrix ttsm);
	void update_alpha(State& in_state, Likelihood& likelihood, double temperature, Data& d, LogicalMatrix row_is_column_anc, NumericMatrix ttsm);
	void update_log_beta(State& in_state, Likelihood& likelihood, double temperature, Data& d, LogicalMatrix row_is_column_anc, NumericMatrix ttsm);
	void update_phi(State& in_state, Likelihood& likelihood, double temperature, Data& d, LogicalMatrix row_is_column_anc, NumericMatrix ttsm);
	void update_logit_mean_f(State& in_state, Likelihood& likelihood, double temperature, Data& d, LogicalMatrix row_is_column_anc, NumericMatrix ttsm);
	void update_log_alpha_plus_beta_f(State& in_state, Likelihood& likelihood, double temperature, Data& d, LogicalMatrix row_is_column_anc, NumericMatrix ttsm);
	void update_logit_mean_g(State& in_state, Likelihood& likelihood, double temperature, Data& d, LogicalMatrix row_is_column_anc, NumericMatrix ttsm);
	void update_log_alpha_plus_beta_g(State& in_state, Likelihood& likelihood, double temperature, Data& d, LogicalMatrix row_is_column_anc, NumericMatrix ttsm);
	void update(State& in_state, Likelihood& likelihood, double temperature, Data& d, LogicalMatrix row_is_column_anc, NumericMatrix ttsm, bool fix_phi);
	void update_phi_and_inflexion(State& in_state, Likelihood& likelihood, double temperature, Data& d, LogicalMatrix row_is_column_anc, NumericMatrix ttsm);
};


#endif
