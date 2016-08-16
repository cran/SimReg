#include "Chain.h"

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
	double temperature = 1.0,
	bool record_x = false,
	bool fix_phi = false,
	bool record_model_likelihoods = false
) {

	int recorded_iterations = (iterations - 1)/thin + 1;
	LogicalVector gamma_trace(recorded_iterations);
	NumericVector likelihood_trace(recorded_iterations);
	NumericVector alpha_trace(recorded_iterations);
	NumericVector log_beta_trace(recorded_iterations);
	NumericVector alpha_star_trace(recorded_iterations);
	NumericVector logit_mean_f_trace(recorded_iterations);
	NumericVector log_alpha_plus_beta_f_trace(recorded_iterations);
	NumericVector logit_mean_g_trace(recorded_iterations);
	NumericVector log_alpha_plus_beta_g_trace(recorded_iterations);
    IntegerMatrix phi_trace(recorded_iterations, phi_size);
    IntegerMatrix proposed_phi_trace(recorded_iterations, phi_size);

	NumericMatrix x_trace;
	NumericMatrix s_phi_trace;
	NumericMatrix s_x_trace;
	NumericVector g0_likelihood;
	NumericVector g1_likelihood;

	if (record_x) {
		NumericMatrix _x_trace(recorded_iterations, d.h.num_cases);
		NumericMatrix _s_phi_trace(recorded_iterations, d.h.num_cases);
		NumericMatrix _s_x_trace(recorded_iterations, d.h.num_cases);
		x_trace = _x_trace;
		s_phi_trace = _s_phi_trace;
		s_x_trace = _s_x_trace;
	}

	if (record_model_likelihoods) {
		NumericVector _g0_likelihood(recorded_iterations);
		NumericVector _g1_likelihood(recorded_iterations);
		g0_likelihood = _g0_likelihood;
		g1_likelihood = _g1_likelihood;
	}

	NumericVector gamma_lik_trace(recorded_iterations);
	NumericVector alpha_star_lik_trace(recorded_iterations);
	NumericVector alpha_lik_trace(recorded_iterations);
	NumericVector log_beta_lik_trace(recorded_iterations);
	NumericVector phi_lik_trace(recorded_iterations);
	NumericVector logit_mean_f_lik_trace(recorded_iterations);
	NumericVector log_alpha_plus_beta_f_lik_trace(recorded_iterations);
	NumericVector logit_mean_g_lik_trace(recorded_iterations);
	NumericVector log_alpha_plus_beta_g_lik_trace(recorded_iterations);
	NumericVector y_lik_trace(recorded_iterations);
	NumericVector iteration(recorded_iterations);

	for (int i = 0; i < iterations; i++) {
		updater.update(
			cur_state,
			likelihood,
			temperature,
			d,
			row_is_column_anc,
			ttsm,
			fix_phi
		);

		if (i % thin == 0) {
			int recording_it = i / thin;
			iteration[recording_it] = i;

			likelihood_trace[recording_it] = cur_state.get_total_lik(row_is_column_anc, ttsm, likelihood, d, temperature);
			gamma_trace[recording_it] = cur_state.gamma;
			alpha_trace[recording_it] = cur_state.alpha;
			log_beta_trace[recording_it] = cur_state.log_beta;
			alpha_star_trace[recording_it] = cur_state.alpha_star;
			logit_mean_f_trace[recording_it] = cur_state.logit_mean_f;
			log_alpha_plus_beta_f_trace[recording_it] = cur_state.log_alpha_plus_beta_f;
			logit_mean_g_trace[recording_it] = cur_state.logit_mean_g;
			log_alpha_plus_beta_g_trace[recording_it] = cur_state.log_alpha_plus_beta_g;
			for (int j = 0; j < phi_size; j++) {
				phi_trace(recording_it, j) = cur_state.phi[j];
				proposed_phi_trace(recording_it, j) = cur_state.phi_proposed[j];
			}

			if (record_model_likelihoods) {
				g0_likelihood[recording_it] = likelihood.get_total_lik(
					d.y,
					cur_state._x,
					d.g,
					cur_state.alpha_star,
					cur_state.alpha,
					cur_state.log_beta,
					cur_state.logit_mean_f,
					cur_state.log_alpha_plus_beta_f,
					cur_state.logit_mean_g,
					cur_state.log_alpha_plus_beta_g,
					cur_state.phi,
					false
				);
				g1_likelihood[recording_it] = likelihood.get_total_lik(
					d.y,
					cur_state._x,
					d.g,
					cur_state.alpha_star,
					cur_state.alpha,
					cur_state.log_beta,
					cur_state.logit_mean_f,
					cur_state.log_alpha_plus_beta_f,
					cur_state.logit_mean_g,
					cur_state.log_alpha_plus_beta_g,
					cur_state.phi,
					true	
				);
			}

			if (record_x) {
				for (int j = 0; j < d.h.num_cases; j++) {
					x_trace(recording_it, j) = cur_state._x[j];
					s_phi_trace(recording_it, j) = cur_state._s.first[j];
					s_x_trace(recording_it, j) = cur_state._s.second[j];
				}
			}

			gamma_lik_trace[recording_it] = cur_state.cur_gamma_lik;
			alpha_star_lik_trace[recording_it] = cur_state.cur_alpha_star_lik;
			alpha_lik_trace[recording_it] = cur_state.cur_alpha_lik;
			log_beta_lik_trace[recording_it] = cur_state.cur_log_beta_lik;
			phi_lik_trace[recording_it] = cur_state.cur_phi_lik;
			logit_mean_f_lik_trace[recording_it] = cur_state.cur_logit_mean_f_lik;
			log_alpha_plus_beta_f_lik_trace[recording_it] = cur_state.cur_log_alpha_plus_beta_f_lik;
			logit_mean_g_lik_trace[recording_it] = cur_state.cur_logit_mean_g_lik;
			log_alpha_plus_beta_g_lik_trace[recording_it] = cur_state.cur_log_alpha_plus_beta_g_lik;
			y_lik_trace[recording_it] = cur_state.cur_y_lik;
		}
	}

	return List::create(
		Named("iteration") = iteration,
		Named("likelihoods") = List::create(
			Named("total_gamma0_likelihood") = g0_likelihood,
			Named("total_gamma1_likelihood") = g1_likelihood,
			Named("total_likelihood") = likelihood_trace,
			Named("gamma_likelihood") = gamma_lik_trace,
			Named("alpha_star_likelihood") = alpha_star_lik_trace,
			Named("alpha_likelihood") = alpha_lik_trace,
			Named("log_beta_likelihood") = log_beta_lik_trace,
			Named("phi_likelihood") = phi_lik_trace,
			Named("logit_mean_f_likelihood") = logit_mean_f_lik_trace,
			Named("log_alpha_plus_beta_f_likelihood") = log_alpha_plus_beta_f_lik_trace,
			Named("logit_mean_g_likelihood") = logit_mean_g_lik_trace,
			Named("log_alpha_plus_beta_g_likelihood") = log_alpha_plus_beta_g_lik_trace,
			Named("y_likelihood") = y_lik_trace
		),
		Named("alpha_star") = alpha_star_trace,
		Named("alpha") = alpha_trace,
		Named("log_beta") = log_beta_trace,
		Named("logit_mean_f") = logit_mean_f_trace,
		Named("log_alpha_plus_beta_f") = log_alpha_plus_beta_f_trace,
		Named("logit_mean_g") = logit_mean_g_trace,
		Named("log_alpha_plus_beta_g") = log_alpha_plus_beta_g_trace,
		Named("phi_vector") = phi_trace,
		Named("proposed_phi_vector") = proposed_phi_trace,
		Named("s") = x_trace,
		Named("s_phi") = s_phi_trace,
		Named("s_x") = s_x_trace,
		Named("gamma") = gamma_trace
	);
}

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
) {
	int recorded_iterations = (iterations + pre_anneal_its - 1)/thin + 1;
	LogicalVector gamma_trace(recorded_iterations);

	NumericVector likelihood_trace(recorded_iterations);
	NumericVector alpha_trace(recorded_iterations);
	NumericVector log_beta_trace(recorded_iterations);
	NumericVector alpha_star_trace(recorded_iterations);
	NumericVector logit_mean_f_trace(recorded_iterations);
	NumericVector log_alpha_plus_beta_f_trace(recorded_iterations);
	NumericVector logit_mean_g_trace(recorded_iterations);
	NumericVector log_alpha_plus_beta_g_trace(recorded_iterations);
    IntegerMatrix phi_trace(recorded_iterations, phi_size);
    IntegerMatrix proposed_phi_trace(recorded_iterations, phi_size);

	NumericVector gamma_lik_trace(recorded_iterations);
	NumericVector alpha_star_lik_trace(recorded_iterations);
	NumericVector alpha_lik_trace(recorded_iterations);
	NumericVector log_beta_lik_trace(recorded_iterations);
	NumericVector phi_lik_trace(recorded_iterations);
	NumericVector logit_mean_f_lik_trace(recorded_iterations);
	NumericVector log_alpha_plus_beta_f_lik_trace(recorded_iterations);
	NumericVector logit_mean_g_lik_trace(recorded_iterations);
	NumericVector log_alpha_plus_beta_g_lik_trace(recorded_iterations);
	NumericVector y_lik_trace(recorded_iterations);
	NumericVector iteration(recorded_iterations);
	NumericVector temperature_trace(recorded_iterations);

	for (int i = 0; i < pre_anneal_its + iterations; i++) {
		double i_temp = i < pre_anneal_its ? initial_temperature : initial_temperature * ((double)(iterations + pre_anneal_its - i)/(double)iterations);

		updater.update(
			cur_state,
			likelihood,
			i_temp,
			d,
			row_is_column_anc,
			ttsm,
			fix_phi
		);

		if (i % thin == 0) {
			int recording_it = i / thin;
			iteration[recording_it] = i;

			temperature_trace[recording_it] = i_temp;
			likelihood_trace[recording_it] = cur_state.get_total_lik(row_is_column_anc, ttsm, likelihood, d, i_temp);
			gamma_trace[recording_it] = cur_state.gamma;
			alpha_trace[recording_it] = cur_state.alpha;
			log_beta_trace[recording_it] = cur_state.log_beta;
			alpha_star_trace[recording_it] = cur_state.alpha_star;
			logit_mean_f_trace[recording_it] = cur_state.logit_mean_f;
			log_alpha_plus_beta_f_trace[recording_it] = cur_state.log_alpha_plus_beta_f;
			logit_mean_g_trace[recording_it] = cur_state.logit_mean_g;
			log_alpha_plus_beta_g_trace[recording_it] = cur_state.log_alpha_plus_beta_g;
			for (int j = 0; j < phi_size; j++) {
				phi_trace(recording_it, j) = cur_state.phi[j];
				proposed_phi_trace(recording_it, j) = cur_state.phi_proposed[j];
			}

			gamma_lik_trace[recording_it] = cur_state.cur_gamma_lik;
			alpha_star_lik_trace[recording_it] = cur_state.cur_alpha_star_lik;
			alpha_lik_trace[recording_it] = cur_state.cur_alpha_lik;
			log_beta_lik_trace[recording_it] = cur_state.cur_log_beta_lik;
			phi_lik_trace[recording_it] = cur_state.cur_phi_lik;
			logit_mean_f_lik_trace[recording_it] = cur_state.cur_logit_mean_f_lik;
			log_alpha_plus_beta_f_lik_trace[recording_it] = cur_state.cur_log_alpha_plus_beta_f_lik;
			logit_mean_g_lik_trace[recording_it] = cur_state.cur_logit_mean_g_lik;
			log_alpha_plus_beta_g_lik_trace[recording_it] = cur_state.cur_log_alpha_plus_beta_g_lik;
			y_lik_trace[recording_it] = cur_state.cur_y_lik;
		}
	}

	return List::create(
		Named("iteration") = iteration,
		Named("likelihoods") = List::create(
			Named("total_likelihood") = likelihood_trace,
			Named("gamma_likelihood") = gamma_lik_trace,
			Named("alpha_star_likelihood") = alpha_star_lik_trace,
			Named("alpha_likelihood") = alpha_lik_trace,
			Named("log_beta_likelihood") = log_beta_lik_trace,
			Named("phi_likelihood") = phi_lik_trace,
			Named("logit_mean_f_likelihood") = logit_mean_f_lik_trace,
			Named("log_alpha_plus_beta_f_likelihood") = log_alpha_plus_beta_f_lik_trace,
			Named("logit_mean_g_likelihood") = logit_mean_g_lik_trace,
			Named("log_alpha_plus_beta_g_likelihood") = log_alpha_plus_beta_g_lik_trace,
			Named("y_likelihood") = y_lik_trace
		),
		Named("alpha_star") = alpha_star_trace,
		Named("alpha") = alpha_trace,
		Named("log_beta") = log_beta_trace,
		Named("logit_mean_f") = logit_mean_f_trace,
		Named("log_alpha_plus_beta_f") = log_alpha_plus_beta_f_trace,
		Named("logit_mean_g") = logit_mean_g_trace,
		Named("log_alpha_plus_beta_g") = log_alpha_plus_beta_g_trace,
		Named("phi_vector") = phi_trace,
		Named("proposed_phi_vector") = proposed_phi_trace,
		Named("temperature") = temperature_trace,
		Named("gamma") = gamma_trace
	);
}


