#include "State.h"

using namespace Rcpp;
using namespace std;

Data::Data(LogicalVector in_y, term_list in_h, NumericVector in_g) : h(in_h) {
	y = in_y;
	h = in_h;
	g = in_g;
}

double State::get_total_lik(LogicalMatrix row_is_column_anc, NumericMatrix ttsm, Likelihood& lik, Data& d, double temperature, bool reparameterise, bool quantile_normalise) {
	return
		get_gamma_lik(lik, temperature) + 
		get_alpha_star_lik(lik, temperature) + 
		get_alpha_lik(lik, temperature) + 
		get_log_beta_lik(lik, temperature) + 
		get_phi_lik(lik, temperature) + 
		get_logit_mean_f_lik(lik, temperature) + 
		get_log_alpha_plus_beta_f_lik(lik, temperature) + 
		get_logit_mean_g_lik(lik, temperature) + 
		get_log_alpha_plus_beta_g_lik(lik, temperature) + 
		get_y_lik(lik, ttsm, row_is_column_anc, d, temperature, reparameterise, quantile_normalise);
}

void State::set_random(
	Likelihood& likelihood,
	Data& d,
	LogicalMatrix row_is_column_anc,
	NumericMatrix ttsm,
	int phi_terms
) {
	gamma = unif_rand() < likelihood.gamma_prior_prob;
	IntegerVector _phi(phi_terms);
	for (int i = 0; i < phi_terms; i++)
		_phi[i] = random_integer(ttsm.ncol());
	phi = _phi;
	alpha_star = norm_rand() * likelihood.alpha_star_sd + likelihood.alpha_star_mean;
	alpha = norm_rand() * likelihood.alpha_sd + likelihood.alpha_mean;
	log_beta = norm_rand() * likelihood.log_beta_sd + likelihood.log_beta_mean;
	logit_mean_f = norm_rand() * likelihood.logit_mean_f_sd + likelihood.logit_mean_f_mean;
	log_alpha_plus_beta_f = norm_rand() * likelihood.log_alpha_plus_beta_f_sd + likelihood.log_alpha_plus_beta_f_mean;
	logit_mean_g = norm_rand() * likelihood.logit_mean_g_sd + likelihood.logit_mean_g_mean;
	log_alpha_plus_beta_g = norm_rand() * likelihood.log_alpha_plus_beta_g_sd + likelihood.log_alpha_plus_beta_g_mean;
}

State State::random(
	Likelihood& likelihood,
	Data& d,
	LogicalMatrix row_is_column_anc,
	NumericMatrix ttsm,
	int phi_terms
) {
	State s;
	s.gamma = unif_rand() < likelihood.gamma_prior_prob;
	IntegerVector _phi(phi_terms);
	for (int i = 0; i < phi_terms; i++)
		_phi[i] = random_integer(ttsm.ncol());
	s.phi = _phi;
	s.alpha_star = norm_rand() * likelihood.alpha_star_sd + likelihood.alpha_star_mean;
	s.alpha = norm_rand() * likelihood.alpha_sd + likelihood.alpha_mean;
	s.log_beta = norm_rand() * likelihood.log_beta_sd + likelihood.log_beta_mean;
	s.logit_mean_f = norm_rand() * likelihood.logit_mean_f_sd + likelihood.logit_mean_f_mean;
	s.log_alpha_plus_beta_f = norm_rand() * likelihood.log_alpha_plus_beta_f_sd + likelihood.log_alpha_plus_beta_f_mean;
	s.logit_mean_g = norm_rand() * likelihood.logit_mean_g_sd + likelihood.logit_mean_g_mean;
	s.log_alpha_plus_beta_g = norm_rand() * likelihood.log_alpha_plus_beta_g_sd + likelihood.log_alpha_plus_beta_g_mean;
	return s;
}

NumericVector State::get_x(NumericMatrix ttsm, LogicalMatrix row_is_column_anc, Data& d, bool reparameterise, bool quantile_normalise) {
	pair<NumericVector, NumericVector> s = get_each_way_sim(
		row_is_column_anc,
		ttsm,
		phi,
		d.h,
		quantile_normalise
	);

	return transform_each_way_sim(
		s,
		logit_mean_f,
		log_alpha_plus_beta_f,
		logit_mean_g,
		log_alpha_plus_beta_g,
		reparameterise
	);
}

double State::get_gamma_lik(Likelihood lik, double temperature) { return lik.get_gamma_lik(gamma) / temperature; }

double State::get_alpha_star_lik(Likelihood lik, double temperature) { return lik.get_alpha_star_lik(alpha_star, gamma) / temperature; }

double State::get_alpha_lik(Likelihood lik, double temperature) { return lik.get_alpha_lik(alpha, gamma) / temperature; }

double State::get_log_beta_lik(Likelihood lik, double temperature) { return lik.get_log_beta_lik(log_beta, gamma) / temperature; }

double State::get_phi_lik(Likelihood lik, double temperature) { return lik.get_phi_lik(phi, gamma) / temperature; }

double State::get_logit_mean_f_lik(Likelihood lik, double temperature) { return lik.get_logit_mean_f_lik(logit_mean_f, gamma) / temperature; }

double State::get_log_alpha_plus_beta_f_lik(Likelihood lik, double temperature) { return lik.get_log_alpha_plus_beta_f_lik(log_alpha_plus_beta_f, gamma) / temperature; }

double State::get_logit_mean_g_lik(Likelihood lik, double temperature) { return lik.get_logit_mean_g_lik(logit_mean_g, gamma) / temperature; }

double State::get_log_alpha_plus_beta_g_lik(Likelihood lik, double temperature) { return lik.get_log_alpha_plus_beta_g_lik(log_alpha_plus_beta_g, gamma) / temperature; }

double State::get_y_lik(Likelihood lik, NumericMatrix ttsm, LogicalMatrix row_is_column_anc, Data& d, double temperature, bool reparameterise, bool quantile_normalise) { return lik.get_y_lik(d.y, transform_each_way_sim(get_each_way_sim(
			row_is_column_anc,
			ttsm,
			phi,
			d.h,
			quantile_normalise
		),
		logit_mean_f,
		log_alpha_plus_beta_f,
		logit_mean_g,
		log_alpha_plus_beta_g, reparameterise), d.g, alpha_star, alpha, log_beta, gamma) / temperature; }


