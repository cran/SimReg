#include "Likelihood.h"

using namespace Rcpp;
using namespace std;

Likelihood::Likelihood(
	double in_gamma_prior_prob,
	double in_alpha_star_mean,
	double in_alpha_mean,
	double in_alpha_star_sd,
	double in_alpha_sd,
	double in_log_beta_mean,
	double in_log_beta_sd,
	double in_logit_mean_f_mean,
	double in_logit_mean_f_sd,
	double in_log_alpha_plus_beta_f_mean,
	double in_log_alpha_plus_beta_f_sd,
	double in_logit_mean_g_mean,
	double in_logit_mean_g_sd,
	double in_log_alpha_plus_beta_g_mean,
	double in_log_alpha_plus_beta_g_sd,
	double in_pseudo_alpha_star_mean,
	double in_pseudo_alpha_mean,
	double in_pseudo_alpha_star_sd,
	double in_pseudo_alpha_sd,
	double in_pseudo_log_beta_mean,
	double in_pseudo_log_beta_sd,
	double in_pseudo_logit_mean_f_mean,
	double in_pseudo_logit_mean_f_sd,
	double in_pseudo_log_alpha_plus_beta_f_mean,
	double in_pseudo_log_alpha_plus_beta_f_sd,
	double in_pseudo_logit_mean_g_mean,
	double in_pseudo_logit_mean_g_sd,
	double in_pseudo_log_alpha_plus_beta_g_mean,
	double in_pseudo_log_alpha_plus_beta_g_sd,
	geom_known_prior& in_phi_lik_func,
	IntegerVector in_pseudo_phi_marginal_prior,
	int total_terms,
	double in_H
) : phi_lik_func(in_phi_lik_func) {
	gamma_prior_prob = in_gamma_prior_prob;
	alpha_star_mean = in_alpha_star_mean;
	alpha_mean = in_alpha_mean;
	alpha_star_sd = in_alpha_star_sd;
	alpha_sd = in_alpha_sd;
	log_beta_mean = in_log_beta_mean;
	log_beta_sd = in_log_beta_sd;
	logit_mean_f_mean = in_logit_mean_f_mean;
	logit_mean_f_sd = in_logit_mean_f_sd;
	log_alpha_plus_beta_f_mean = in_log_alpha_plus_beta_f_mean;
	log_alpha_plus_beta_f_sd = in_log_alpha_plus_beta_f_sd;
	logit_mean_g_mean = in_logit_mean_g_mean;
	logit_mean_g_sd = in_logit_mean_g_sd;
	log_alpha_plus_beta_g_mean = in_log_alpha_plus_beta_g_mean;
	log_alpha_plus_beta_g_sd = in_log_alpha_plus_beta_g_sd;
	pseudo_alpha_star_mean = in_pseudo_alpha_star_mean;
	pseudo_alpha_mean = in_pseudo_alpha_mean;
	pseudo_alpha_star_sd = in_pseudo_alpha_star_sd;
	pseudo_alpha_sd = in_pseudo_alpha_sd;
	pseudo_log_beta_mean = in_pseudo_log_beta_mean;
	pseudo_log_beta_sd = in_pseudo_log_beta_sd;
	pseudo_logit_mean_f_mean = in_pseudo_logit_mean_f_mean;
	pseudo_logit_mean_f_sd = in_pseudo_logit_mean_f_sd;
	pseudo_log_alpha_plus_beta_f_mean = in_pseudo_log_alpha_plus_beta_f_mean;
	pseudo_log_alpha_plus_beta_f_sd = in_pseudo_log_alpha_plus_beta_f_sd;
	pseudo_logit_mean_g_mean = in_pseudo_logit_mean_g_mean;
	pseudo_logit_mean_g_sd = in_pseudo_logit_mean_g_sd;
	pseudo_log_alpha_plus_beta_g_mean = in_pseudo_log_alpha_plus_beta_g_mean;
	pseudo_log_alpha_plus_beta_g_sd = in_pseudo_log_alpha_plus_beta_g_sd;

	pseudo_phi_marginal_prior = in_pseudo_phi_marginal_prior;

	IntegerVector _temp(total_terms);
	for (int i = 0; i < pseudo_phi_marginal_prior.length(); i++)
		_temp[pseudo_phi_marginal_prior[i]]++;
	pseudo_phi_marginal_prior_each = _temp;
	
	H = in_H;
}

double Likelihood::get_gamma_lik(bool gamma) {
	return gamma ? (H + log(gamma_prior_prob)) : (-H + log(1.0-gamma_prior_prob));		
}

double Likelihood::get_alpha_star_lik(double val, bool gamma) {
	return (!gamma) ? log_likelihood_normal(alpha_star_mean, alpha_star_sd, val) : log_likelihood_normal(pseudo_alpha_star_mean, pseudo_alpha_star_sd, val);
}

double Likelihood::get_alpha_lik(double val, bool gamma) {
	return gamma ? log_likelihood_normal(alpha_mean, alpha_sd, val) : log_likelihood_normal(pseudo_alpha_mean, pseudo_alpha_sd, val);
}

double Likelihood::get_log_beta_lik(double val, bool gamma) {
	return gamma ? log_likelihood_normal(log_beta_mean, log_beta_sd, val) : log_likelihood_normal(pseudo_log_beta_mean, pseudo_log_beta_sd, val);
}

double Likelihood::get_phi_lik(IntegerVector phi, bool gamma) {
	if (gamma) {
		return phi_lik_func(phi);
	}
 else {
		double phi_lik_false = - (double)phi.length() * log((double)pseudo_phi_marginal_prior.length());
		for (int node = 0; node < phi.length(); node++) {
			phi_lik_false += log((double)pseudo_phi_marginal_prior_each[phi[node]]);
		}

		return phi_lik_false;
	}

}

double Likelihood::get_logit_mean_f_lik(double val, bool gamma) {
	return gamma ? log_likelihood_normal(logit_mean_f_mean, logit_mean_f_sd, val) : log_likelihood_normal(pseudo_logit_mean_f_mean, pseudo_logit_mean_f_sd, val);
}

double Likelihood::get_log_alpha_plus_beta_f_lik(double val, bool gamma) {
	return gamma ? log_likelihood_normal(log_alpha_plus_beta_f_mean, log_alpha_plus_beta_f_sd, val) : log_likelihood_normal(pseudo_log_alpha_plus_beta_f_mean, pseudo_log_alpha_plus_beta_f_sd, val);
}

double Likelihood::get_logit_mean_g_lik(double val, bool gamma) {
	return gamma ? log_likelihood_normal(logit_mean_g_mean, logit_mean_g_sd, val) : log_likelihood_normal(pseudo_logit_mean_g_mean, pseudo_logit_mean_g_sd, val);
}

double Likelihood::get_log_alpha_plus_beta_g_lik(double val, bool gamma) {
	return gamma ? log_likelihood_normal(log_alpha_plus_beta_g_mean, log_alpha_plus_beta_g_sd, val) : log_likelihood_normal(pseudo_log_alpha_plus_beta_g_mean, pseudo_log_alpha_plus_beta_g_sd, val);
}

double Likelihood::get_y_lik(LogicalVector y, NumericVector x, NumericVector g, double alpha_star, double alpha, double log_beta, bool gamma) {
	double loglik = 0.0;
	int n = y.length();
	if (gamma) {
		for (int i = 0; i < n; i++) {
			double log_odds = alpha + exp(log_beta) * x[i] + g[i];
			loglik += (double)y[i] * log_prob(log_odds) + (1.0 - (double)y[i]) * log_prob(-log_odds);
		}
	}
	else {
		for (int i = 0; i < n; i++) {
			double log_odds = alpha_star + g[i];
			loglik += (double)y[i] * log_prob(log_odds) + (1.0 - (double)y[i]) * log_prob(-log_odds);
		}
	}

	return loglik;
}

double Likelihood::get_total_lik(LogicalVector y, NumericVector x, NumericVector g, double alpha_star, double alpha, double log_beta, double logit_mean_f, double log_alpha_plus_beta_f, double logit_mean_g, double log_alpha_plus_beta_g, IntegerVector phi, bool gamma) {
	return
	get_gamma_lik(gamma) +
	get_alpha_star_lik(alpha_star, gamma) +
	get_alpha_lik(alpha, gamma) +
	get_log_beta_lik(log_beta, gamma) +
	get_phi_lik(phi, gamma) +
	get_logit_mean_f_lik(logit_mean_f, gamma) +
	get_log_alpha_plus_beta_f_lik(log_alpha_plus_beta_f, gamma) +
	get_logit_mean_g_lik(logit_mean_g, gamma) +
	get_log_alpha_plus_beta_g_lik(log_alpha_plus_beta_g, gamma) +
	get_y_lik(y, x, g, alpha_star, alpha, log_beta, gamma);
}


