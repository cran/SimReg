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
	bool reparameterise = false,
	bool quantile_normalise = false
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

	NumericMatrix x_trace;
	NumericMatrix s_phi_trace;
	NumericMatrix s_h_trace;

	if (record_x) {
		NumericMatrix _x_trace(recorded_iterations, d.h.num_cases);
		NumericMatrix _s_phi_trace(recorded_iterations, d.h.num_cases);
		NumericMatrix _s_h_trace(recorded_iterations, d.h.num_cases);
		x_trace = _x_trace;
		s_phi_trace = _s_phi_trace;
		s_h_trace = _s_h_trace;
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

	for (int i = 0; i < iterations; i++) {
		updater.update(
			cur_state,
			likelihood,
			temperature,
			d,
			row_is_column_anc,
			ttsm,
			fix_phi,
			reparameterise,
			quantile_normalise
		);

		if (i % thin == 0) {
			int recording_it = i / thin;
			likelihood_trace[recording_it] = cur_state.get_total_lik(row_is_column_anc, ttsm, likelihood, d, temperature, reparameterise, quantile_normalise);
			gamma_trace[recording_it] = cur_state.gamma;
			alpha_trace[recording_it] = cur_state.alpha;
			log_beta_trace[recording_it] = cur_state.log_beta;
			alpha_star_trace[recording_it] = cur_state.alpha_star;
			logit_mean_f_trace[recording_it] = cur_state.logit_mean_f;
			log_alpha_plus_beta_f_trace[recording_it] = cur_state.log_alpha_plus_beta_f;
			logit_mean_g_trace[recording_it] = cur_state.logit_mean_g;
			log_alpha_plus_beta_g_trace[recording_it] = cur_state.log_alpha_plus_beta_g;
			for (int j = 0; j < phi_size; j++)
				phi_trace(recording_it, j) = cur_state.phi[j];

			if (record_x) {
				for (int j = 0; j < d.h.num_cases; j++) {
					x_trace(recording_it, j) = cur_state._x[j];
					s_phi_trace(recording_it, j) = cur_state._s.first[j];
					s_h_trace(recording_it, j) = cur_state._s.second[j];
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
		Named("liks") = List::create(
			Named("likelihood") = likelihood_trace,
			Named("gamma") = gamma_lik_trace,
			Named("alpha_star") = alpha_star_lik_trace,
			Named("alpha") = alpha_lik_trace,
			Named("log_beta") = log_beta_lik_trace,
			Named("phi") = phi_lik_trace,
			Named("logit_mean_f") = logit_mean_f_lik_trace,
			Named("log_alpha_plus_beta_f") = log_alpha_plus_beta_f_lik_trace,
			Named("logit_mean_g") = logit_mean_g_lik_trace,
			Named("log_alpha_plus_beta_g") = log_alpha_plus_beta_g_lik_trace,
			Named("y") = y_lik_trace
		),
		Named("alpha_star") = alpha_star_trace,
		Named("alpha") = alpha_trace,
		Named("log_beta") = log_beta_trace,
		Named("logit_mean_f") = logit_mean_f_trace,
		Named("log_alpha_plus_beta_f") = log_alpha_plus_beta_f_trace,
		Named("logit_mean_g") = logit_mean_g_trace,
		Named("log_alpha_plus_beta_g") = log_alpha_plus_beta_g_trace,
		Named("phi") = phi_trace,
		Named("x") = x_trace,
		Named("s_phi") = s_phi_trace,
		Named("s_h") = s_h_trace,
		Named("gamma") = gamma_trace,
		Named("alpha_star_lik") = alpha_star_lik_trace,
		Named("log_beta_lik") = log_beta_lik_trace,
		Named("phi_lik") = phi_lik_trace,
		Named("logit_mean_f_lik") = logit_mean_f_lik_trace,
		Named("y_lik") = y_lik_trace,
		Named("H") = likelihood.H
	);

}


