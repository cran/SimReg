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
		cur_gamma_lik +
		cur_alpha_star_lik +
		cur_alpha_lik +
		cur_log_beta_lik +
		cur_logit_mean_f_lik +
		cur_log_alpha_plus_beta_f_lik +
		cur_logit_mean_g_lik +
		cur_log_alpha_plus_beta_g_lik +
		cur_phi_lik +
		cur_y_lik
	;
}

void State::initialise(
	Likelihood& likelihood,
	Data& d,
	LogicalMatrix row_is_column_anc,
	NumericMatrix ttsm
) {
	cur_gamma_lik = get_gamma_lik(likelihood, 1.0);
	cur_alpha_star_lik = get_alpha_star_lik(likelihood, 1.0);
	cur_alpha_lik = get_alpha_lik(likelihood, 1.0);
	cur_log_beta_lik = get_log_beta_lik(likelihood, 1.0);
	cur_logit_mean_f_lik = get_logit_mean_f_lik(likelihood, 1.0);
	cur_log_alpha_plus_beta_f_lik = get_log_alpha_plus_beta_f_lik(likelihood, 1.0);
	cur_logit_mean_g_lik = get_logit_mean_g_lik(likelihood, 1.0);
	cur_log_alpha_plus_beta_g_lik = get_log_alpha_plus_beta_g_lik(likelihood, 1.0);
	cur_phi_lik = get_phi_lik(likelihood, 1.0);
	cur_y_lik = get_y_lik(likelihood, ttsm, row_is_column_anc, d, 1.0, true, false);
	_s = get_s(ttsm, row_is_column_anc, d);
	_x = get_x(ttsm, row_is_column_anc, d);
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

	initialise(likelihood, d, row_is_column_anc, ttsm);
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

	s.initialise(likelihood, d, row_is_column_anc, ttsm);
	return s;
}

pair<NumericVector, NumericVector> State::get_s(NumericMatrix ttsm, LogicalMatrix row_is_column_anc, Data& d, bool reparameterise, bool quantile_normalise) {
	pair<NumericVector, NumericVector> s = get_each_way_sim(
		row_is_column_anc,
		ttsm,
		phi,
		d.h,
		quantile_normalise
	);
	return s;
}

NumericVector State::get_x(NumericMatrix ttsm, LogicalMatrix row_is_column_anc, Data& d, bool reparameterise, bool quantile_normalise) {
	NumericVector x = transform_each_way_sim(
		_s,
		logit_mean_f,
		log_alpha_plus_beta_f,
		logit_mean_g,
		log_alpha_plus_beta_g,
		reparameterise
	);
	return x;
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


