#include <Rcpp.h>
#include <cmath>
#include "Utils.h"
#include "DeltaLikelihood.h"

#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

using namespace Rcpp;
using namespace std;

class Likelihood
{
private:
public:
	Likelihood(
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
		double in_H = 0.0
	);

	double gamma_prior_prob;
	double alpha_star_mean;
	double alpha_mean;
	double alpha_star_sd;
	double alpha_sd;
	double log_beta_mean;
	double log_beta_sd;
	double logit_mean_f_mean;
	double logit_mean_f_sd;
	double log_alpha_plus_beta_f_mean;
	double log_alpha_plus_beta_f_sd;
	double logit_mean_g_mean;
	double logit_mean_g_sd;
	double log_alpha_plus_beta_g_mean;
	double log_alpha_plus_beta_g_sd;
	double pseudo_alpha_star_mean;
	double pseudo_alpha_mean;
	double pseudo_alpha_star_sd;
	double pseudo_alpha_sd;
	double pseudo_log_beta_mean;
	double pseudo_log_beta_sd;
	double pseudo_logit_mean_f_mean;
	double pseudo_logit_mean_f_sd;
	double pseudo_log_alpha_plus_beta_f_mean;
	double pseudo_log_alpha_plus_beta_f_sd;
	double pseudo_logit_mean_g_mean;
	double pseudo_logit_mean_g_sd;
	double pseudo_log_alpha_plus_beta_g_mean;
	double pseudo_log_alpha_plus_beta_g_sd;
	geom_known_prior phi_lik_func;
	IntegerVector pseudo_phi_marginal_prior_each;
	IntegerVector pseudo_phi_marginal_prior;
	double H;

	double get_gamma_lik(bool gamma);
	double get_alpha_star_lik(double alpha_star, bool gamma);
	double get_alpha_lik(double alpha, bool gamma);
	double get_log_beta_lik(double log_beta, bool gamma);
	double get_phi_lik(IntegerVector phi, bool gamma);
	double get_logit_mean_f_lik(double logit_mean_f, bool gamma);
	double get_log_alpha_plus_beta_f_lik(double log_alpha_plus_beta_f, bool gamma);
	double get_logit_mean_g_lik(double logit_mean_g, bool gamma);
	double get_log_alpha_plus_beta_g_lik(double log_alpha_plus_beta_g, bool gamma);
	double get_y_lik(LogicalVector y, NumericVector x, NumericVector g, double alpha_star, double alpha, double log_beta, bool gamma);
	double get_total_lik(LogicalVector y, NumericVector x, NumericVector g, double alpha_star, double alpha, double log_beta, double logit_mean_f, double log_alpha_plus_beta_f, double logit_mean_g, double log_alpha_plus_beta_g, IntegerVector phi, bool gamma);

};

#endif
