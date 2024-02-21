#include <Rcpp.h>
#include <cmath>
#include "utils.h"

using namespace Rcpp;
using namespace std;

struct FPrior
{
	double alpha_mean;
	double alpha_sd;
	double log_beta_mean;
	double log_beta_sd;
	double logit_f_mean_mean;
	double logit_f_mean_sd;
	double log_f_a_plus_b_mean;
	double log_f_a_plus_b_sd;

	FPrior(
		double in_alpha_mean,
		double in_alpha_sd,
		double in_log_beta_mean,
		double in_log_beta_sd,
		double in_logit_f_mean_mean,
		double in_logit_f_mean_sd,
		double in_log_f_a_plus_b_mean,
		double in_log_f_a_plus_b_sd
	) {
		alpha_mean = in_alpha_mean;
		alpha_sd = in_alpha_sd;
		log_beta_mean = in_log_beta_mean;
		log_beta_sd = in_log_beta_sd;
		logit_f_mean_mean = in_logit_f_mean_mean;
		logit_f_mean_sd = in_logit_f_mean_sd;
		log_f_a_plus_b_mean = in_log_f_a_plus_b_mean;
		log_f_a_plus_b_sd = in_log_f_a_plus_b_sd;
	}

	FPrior(FPrior& p) {
		alpha_mean = p.alpha_mean;
		alpha_sd = p.alpha_sd;
		log_beta_mean = p.log_beta_mean;
		log_beta_sd = p.log_beta_sd;
		logit_f_mean_mean = p.logit_f_mean_mean;
		logit_f_mean_sd = p.logit_f_mean_sd;
		log_f_a_plus_b_mean = p.log_f_a_plus_b_mean;
		log_f_a_plus_b_sd = p.log_f_a_plus_b_sd;
	}
};

struct FState
{
	NumericVector x_values;
	IntegerVector num_y0_phi;
	IntegerVector num_y1_phi;

	NumericVector s;
	NumericVector lo;

	double likelihood;
	double prior_density;

	double alpha;
	double log_beta;
	double logit_f_mean;
	double log_f_a_plus_b;

	FPrior prior;
	int gran;

	FState(
		NumericVector x,
		IntegerVector y0phi,
		IntegerVector y1phi,
		FPrior in_prior
	) : x_values(x), num_y0_phi(y0phi), num_y1_phi(y1phi), prior(in_prior) {
		gran = x_values.length();

		alpha = norm_rand() * prior.alpha_sd + prior.alpha_mean;
		log_beta = norm_rand() * prior.log_beta_sd + prior.log_beta_mean;
		logit_f_mean = norm_rand() * prior.logit_f_mean_sd + prior.logit_f_mean_mean;
		log_f_a_plus_b = norm_rand() * prior.log_f_a_plus_b_sd + prior.log_f_a_plus_b_mean;

		prior_density = 0.0;
		prior_density += R::dnorm(alpha, prior.alpha_mean, prior.alpha_sd, 1);
		prior_density += R::dnorm(log_beta, prior.log_beta_mean, prior.log_beta_sd, 1);
		prior_density += R::dnorm(logit_f_mean, prior.logit_f_mean_mean, prior.logit_f_mean_sd, 1);
		prior_density += R::dnorm(log_f_a_plus_b, prior.log_f_a_plus_b_mean, prior.log_f_a_plus_b_sd, 1);

		s = NumericVector(gran);
		lo = NumericVector(gran);
		set_s();
	}

	void set_s() {
		double mean = 1.0-1.0/(exp(logit_f_mean) + 1.0);
		double inv_var = exp(log_f_a_plus_b);
		for (int i = 0; i < gran; i++)
			s[i] = max(0.0, min(1.0, inv_var * x_values[i] + 0.5 - inv_var * mean)); 
		set_lo();
	}

	void set_lo() {
		double beta = exp(log_beta);
		for (int i = 0; i < gran; i++)
			lo[i] = alpha + beta * s[i];
		set_likelihood();
	}

	void set_likelihood() {
		double result = 0.0;
		for (int i = 0; i < gran; i++) {
			result += (double)num_y1_phi[i] * log_prob(lo[i]);
			result += (double)num_y0_phi[i] * log_prob(-lo[i]);
		}
		likelihood = result;
	}

	void set_alpha(double value) {
		prior_density -= R::dnorm(alpha, prior.alpha_mean, prior.alpha_sd, 1);
		alpha = value;
		prior_density += R::dnorm(alpha, prior.alpha_mean, prior.alpha_sd, 1);
		set_lo();
	}

	void set_log_beta(double value) {
		prior_density -= R::dnorm(log_beta, prior.log_beta_mean, prior.log_beta_sd, 1);
		log_beta = value;
		prior_density += R::dnorm(log_beta, prior.log_beta_mean, prior.log_beta_sd, 1);
		set_lo();
	}

	void set_logit_f_mean(double value) {
		prior_density -= R::dnorm(logit_f_mean, prior.logit_f_mean_mean, prior.logit_f_mean_sd, 1);
		logit_f_mean = value;
		prior_density += R::dnorm(logit_f_mean, prior.logit_f_mean_mean, prior.logit_f_mean_sd, 1);
		set_s();
	}

	void set_log_f_a_plus_b(double value) {
		prior_density -= R::dnorm(log_f_a_plus_b, prior.log_f_a_plus_b_mean, prior.log_f_a_plus_b_sd, 1);
		log_f_a_plus_b = value;
		prior_density += R::dnorm(log_f_a_plus_b, prior.log_f_a_plus_b_mean, prior.log_f_a_plus_b_sd, 1);
		set_s();
	}
};

// [[Rcpp::export]]
double f_ML(
	NumericVector x_values,
	IntegerVector num_y0_phi,
	IntegerVector num_y1_phi,
	NumericVector t,
	double log_scale_tolerance,
	int min_samples,
	int max_samples,
	double min_log_ML,
	double alpha_mean,
	double alpha_sd,
	double log_beta_mean,
	double log_beta_sd,
	double logit_f_mean,
	double logit_f_sd,
	double log_f_a_plus_b_mean,
	double log_f_a_plus_b_sd,
	double alpha_prop_sd,
	double log_beta_prop_sd,
	double logit_f_mean_prop_sd,
	double log_f_a_plus_b_prop_sd
) {
	FPrior prior(
		alpha_mean,
		alpha_sd,
		log_beta_mean,
		log_beta_sd,
		logit_f_mean,
		logit_f_sd,
		log_f_a_plus_b_mean,
		log_f_a_plus_b_sd
	);

	int n = t.length();

	double current_factor = 1.0;
	double highest_multiplier = -INFINITY;
	double current_factor2 = 1.0;
	double highest_multiplier2 = -INFINITY;

	int i = 0;
	double log_sum_so_far = -INFINITY;

	double E_x;
	double sd;

	do {

		double sum_ratios = 0;
		FState state(
			x_values,
			num_y0_phi,
			num_y1_phi,
			prior
		);

		for (int j = 1; j < n; j++) {
			double cur;
			double prp;

			sum_ratios += (t[j]-t[j-1]) * state.likelihood;

			double old_alpha = state.alpha;
			cur = state.prior_density + t[j] * state.likelihood;
			state.set_alpha(state.alpha + norm_rand() * alpha_prop_sd);
			prp = state.prior_density + t[j] * state.likelihood;
			
			if (log(unif_rand()) > prp-cur) {
				state.set_alpha(old_alpha);
			}
			
			double old_log_beta = state.log_beta;
			cur = state.prior_density + t[j] * state.likelihood;
			state.set_log_beta(state.log_beta + norm_rand() * log_beta_prop_sd);
			prp = state.prior_density + t[j] * state.likelihood;
			
			if (log(unif_rand()) > prp-cur) {
				state.set_log_beta(old_log_beta);
			}
			
			double old_logit_f_mean = state.logit_f_mean;
			cur = state.prior_density + t[j] * state.likelihood;
			state.set_logit_f_mean(state.logit_f_mean + norm_rand() * logit_f_mean_prop_sd);
			prp = state.prior_density + t[j] * state.likelihood;
			
			if (log(unif_rand()) > prp-cur) {
				state.set_logit_f_mean(old_logit_f_mean);
			}

			double old_log_f_a_plus_b = state.log_f_a_plus_b;
			cur = state.prior_density + t[j] * state.likelihood;
			state.set_log_f_a_plus_b(state.log_f_a_plus_b + norm_rand() * log_f_a_plus_b_prop_sd);
			prp = state.prior_density + t[j] * state.likelihood;
			
			if (log(unif_rand()) > prp-cur) {
				state.set_log_f_a_plus_b(old_log_f_a_plus_b);
			}
		}

		if (sum_ratios > highest_multiplier) {
			current_factor = 1.0 + current_factor * exp(highest_multiplier-sum_ratios);
			highest_multiplier = sum_ratios;
			current_factor2 = 1.0 + current_factor2 * exp(highest_multiplier2-sum_ratios*2);
			highest_multiplier2 = sum_ratios * 2;
		}
		else {
			current_factor = exp(sum_ratios - highest_multiplier) + current_factor;
			current_factor2 = exp(sum_ratios * 2 - highest_multiplier2) + current_factor2;
		}

		i++;
		log_sum_so_far = log(current_factor) + highest_multiplier;
		
		E_x = (log_sum_so_far - log((double)i));
		double E_x_2 = 2 * E_x;
		double E_x2 = log(current_factor2) + highest_multiplier2 - log((double)i);
		sd = (E_x_2 + log(exp(E_x2-E_x_2)-1.0) - log((double)i)) / 2;
	}
	while ((i < min_samples) || ((i < max_samples) && ((log_sum_so_far-log((double)i)) > min_log_ML) && (sd-E_x > log_scale_tolerance)));

	return log_sum_so_far - log((double)i);
}
