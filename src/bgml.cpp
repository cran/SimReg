#include <Rcpp.h>
#include <cmath>
#include "utils.h"

using namespace Rcpp;
using namespace std;

struct BGPrior
{
	double alpha_mean;
	double alpha_sd;

	BGPrior(
		double in_alpha_mean,
		double in_alpha_sd
	) {
		alpha_mean = in_alpha_mean;
		alpha_sd = in_alpha_sd;
	}

	BGPrior(BGPrior& p) {
		alpha_mean = p.alpha_mean;
		alpha_sd = p.alpha_sd;
	}
};

struct BGState
{
	double likelihood;
	double prior_density;

	int y0;
	int y1;
	double alpha;

	BGPrior prior;

	BGState(
		int in_y0,
		int in_y1,
		BGPrior in_prior
	) : prior(in_prior) {
		y0 = in_y0;
		y1 = in_y1;
		alpha = norm_rand() * prior.alpha_sd + prior.alpha_mean;

		prior_density = 0.0;
		prior_density += R::dnorm(alpha, prior.alpha_mean, prior.alpha_sd, 1);
	}

	void set_likelihood() {
		likelihood = (double)y1 * log_prob(alpha) + y0 * log_prob(-alpha);
	}

	void set_alpha(double value) {
		prior_density -= R::dnorm(alpha, prior.alpha_mean, prior.alpha_sd, 1);
		alpha = value;
		prior_density += R::dnorm(alpha, prior.alpha_mean, prior.alpha_sd, 1);
		set_likelihood();
	}
};

// [[Rcpp::export]]
double bg_ML(
	int y0,
	int y1,
	NumericVector t,
	int n_samples,
	double alpha_mean,
	double alpha_sd,
	double alpha_prop_sd
) {
	BGPrior prior(
		alpha_mean,
		alpha_sd
	);

	int n = t.length();

	double current_factor = 1.0;
	double highest_multiplier = -INFINITY;

	for (int i = 0; i < n_samples; i++) {
		double sum_ratios = 0;
		BGState state(
			y0,
			y1,
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
		}

		if (sum_ratios > highest_multiplier) {
			current_factor = 1.0 + current_factor * exp(highest_multiplier-sum_ratios);
			highest_multiplier = sum_ratios;
		}
		else {
			current_factor = exp(sum_ratios - highest_multiplier) + current_factor;
		}
	}

	return log(current_factor) + highest_multiplier - log((double)n_samples);
}

