#include "Update.h"

using namespace Rcpp;
using namespace std;

void Update::update_gamma(State& in_state, Likelihood& likelihood, double temperature, Data& d, LogicalMatrix row_is_column_anc, NumericMatrix ttsm) {
	double alt_gamma_lik = likelihood.get_gamma_lik(!in_state.gamma) / temperature;
	double alt_alpha_star_lik = likelihood.get_alpha_star_lik(in_state.alpha_star, !in_state.gamma) / temperature;
	double alt_alpha_lik = likelihood.get_alpha_lik(in_state.alpha, !in_state.gamma) / temperature;
	double alt_log_beta_lik = likelihood.get_log_beta_lik(in_state.log_beta, !in_state.gamma) / temperature;
	double alt_phi_lik = likelihood.get_phi_lik(in_state.phi, !in_state.gamma) / temperature;
	double alt_logit_mean_f_lik = likelihood.get_logit_mean_f_lik(in_state.logit_mean_f, !in_state.gamma) / temperature;
	double alt_log_alpha_plus_beta_f_lik = likelihood.get_log_alpha_plus_beta_f_lik(in_state.log_alpha_plus_beta_f, !in_state.gamma) / temperature;
	double alt_logit_mean_g_lik = likelihood.get_logit_mean_g_lik(in_state.logit_mean_g, !in_state.gamma) / temperature;
	double alt_log_alpha_plus_beta_g_lik = likelihood.get_log_alpha_plus_beta_g_lik(in_state.log_alpha_plus_beta_g, !in_state.gamma) / temperature;
	double alt_y_lik = likelihood.get_y_lik(d.y, in_state.get_x(ttsm, row_is_column_anc, d), d.g, in_state.alpha_star, in_state.alpha, in_state.log_beta, !in_state.gamma) / temperature;

	double alt_total_lik = 
		alt_gamma_lik +
		alt_alpha_star_lik +
		alt_alpha_lik +
		alt_log_beta_lik +
		alt_logit_mean_f_lik +
		alt_log_alpha_plus_beta_f_lik +
		alt_logit_mean_g_lik +
		alt_log_alpha_plus_beta_g_lik +
		alt_phi_lik +
		alt_y_lik
	;
	
	double cur_total_lik = 
		in_state.cur_gamma_lik +
		in_state.cur_alpha_star_lik +
		in_state.cur_alpha_lik +
		in_state.cur_log_beta_lik +
		in_state.cur_logit_mean_f_lik +
		in_state.cur_log_alpha_plus_beta_f_lik +
		in_state.cur_logit_mean_g_lik +
		in_state.cur_log_alpha_plus_beta_g_lik +
		in_state.cur_phi_lik +
		in_state.cur_y_lik
	;

	if (log(unif_rand()) < log_prob(alt_total_lik - cur_total_lik)) {
		in_state.gamma = !in_state.gamma;
		in_state.cur_gamma_lik = alt_gamma_lik;
		in_state.cur_alpha_star_lik = alt_alpha_star_lik;
		in_state.cur_alpha_lik = alt_alpha_lik;
		in_state.cur_log_beta_lik = alt_log_beta_lik;
		in_state.cur_logit_mean_f_lik = alt_logit_mean_f_lik;
		in_state.cur_log_alpha_plus_beta_f_lik = alt_log_alpha_plus_beta_f_lik;
		in_state.cur_logit_mean_g_lik = alt_logit_mean_g_lik;
		in_state.cur_log_alpha_plus_beta_g_lik = alt_log_alpha_plus_beta_g_lik;
		in_state.cur_phi_lik = alt_phi_lik;
		in_state.cur_y_lik = alt_y_lik;
	}
}

void Update::update_alpha_star(State& in_state, Likelihood& likelihood, double temperature, Data& d, LogicalMatrix row_is_column_anc, NumericMatrix ttsm) {
	if (in_state.gamma) {
		in_state.alpha_star = norm_rand() * likelihood.pseudo_alpha_star_sd + likelihood.pseudo_alpha_star_mean;
		in_state.cur_alpha_star_lik = in_state.get_alpha_star_lik(likelihood);
	} else {

		double proposed_alpha_star = in_state.alpha_star + norm_rand() * alpha_star_proposal_sd;
		
		double proposed_alpha_star_lik = likelihood.get_alpha_star_lik(proposed_alpha_star, in_state.gamma) / temperature;
		
		double proposed_y_lik = likelihood.get_y_lik(d.y, in_state.get_x(ttsm, row_is_column_anc, d), d.g, proposed_alpha_star, in_state.alpha, in_state.log_beta, in_state.gamma) / temperature;

		if (log(unif_rand()) < (proposed_alpha_star_lik + proposed_y_lik - (in_state.get_alpha_star_lik(likelihood, temperature) + in_state.get_y_lik(likelihood, ttsm, row_is_column_anc, d, temperature)))) {
			in_state.alpha_star = proposed_alpha_star;
			in_state.cur_alpha_star_lik = proposed_alpha_star_lik;
			in_state.cur_y_lik = proposed_y_lik;
		}
	}
}

void Update::update_alpha(State& in_state, Likelihood& likelihood, double temperature, Data& d, LogicalMatrix row_is_column_anc, NumericMatrix ttsm) {
	if (in_state.gamma) {
		double proposed_alpha = in_state.alpha + norm_rand() * alpha_proposal_sd;
		
		double proposed_alpha_lik = likelihood.get_alpha_lik(proposed_alpha, in_state.gamma) / temperature;
		
		double proposed_y_lik = likelihood.get_y_lik(d.y, in_state.get_x(ttsm, row_is_column_anc, d), d.g, in_state.alpha_star, proposed_alpha, in_state.log_beta, in_state.gamma) / temperature;

		if (log(unif_rand()) < ((proposed_alpha_lik + proposed_y_lik) - (in_state.get_alpha_lik(likelihood, temperature) + in_state.get_y_lik(likelihood, ttsm, row_is_column_anc, d, temperature)))) {
			in_state.alpha = proposed_alpha;
			in_state.cur_alpha_lik = proposed_alpha_lik;
			in_state.cur_y_lik = proposed_y_lik;

		}
	} else {
		in_state.alpha = norm_rand() * likelihood.pseudo_alpha_sd + likelihood.pseudo_alpha_mean;
		in_state.cur_alpha_lik = in_state.get_alpha_lik(likelihood);
	}
}

void Update::update_log_beta(State& in_state, Likelihood& likelihood, double temperature, Data& d, LogicalMatrix row_is_column_anc, NumericMatrix ttsm) {
	if (in_state.gamma) {
		double current_log_beta_lik = in_state.get_log_beta_lik(likelihood, temperature);
		double current_y_lik = in_state.get_y_lik(likelihood, ttsm, row_is_column_anc, d, temperature);
	
		double proposed_log_beta = in_state.log_beta + norm_rand() * log_beta_proposal_sd;
		
		double proposed_log_beta_lik = likelihood.get_log_beta_lik(proposed_log_beta, in_state.gamma) / temperature;
		
		double proposed_y_lik = likelihood.get_y_lik(d.y, in_state.get_x(ttsm, row_is_column_anc, d), d.g, in_state.alpha_star, in_state.alpha, proposed_log_beta, in_state.gamma) / temperature;

		if (log(unif_rand()) < ((proposed_log_beta_lik + proposed_y_lik) - (current_log_beta_lik + current_y_lik))) {
			in_state.log_beta = proposed_log_beta;
			in_state.cur_log_beta_lik = proposed_log_beta_lik;
			in_state.cur_y_lik = proposed_y_lik;
		}
	} else {
		in_state.log_beta = norm_rand() * likelihood.pseudo_log_beta_sd + likelihood.pseudo_log_beta_mean;
		in_state.cur_log_beta_lik = in_state.get_log_beta_lik(likelihood);
	}
}

void Update::update_phi_and_inflexion(State& in_state, Likelihood& likelihood, double temperature, Data& d, LogicalMatrix row_is_column_anc, NumericMatrix ttsm) {
	if (in_state.gamma) {	
		IntegerVector phi_proposed(in_state.phi.length());
		double proposed_logit_mean_f = in_state.logit_mean_f + norm_rand() * logit_mean_f_proposal_sd;
		double proposed_logit_mean_g = in_state.logit_mean_g + norm_rand() * logit_mean_g_proposal_sd;

		for (int k = 0; k < in_state.phi.length(); k++)
			phi_proposed[k] = in_state.phi[k];
		
		int to_change = random_integer(in_state.phi.length());

		phi_proposed[to_change] = phi_jumps[random_integer(phi_jumps.length())];

		double proposed_phi_lik = likelihood.get_phi_lik(phi_proposed, in_state.gamma) / temperature;
		double proposed_logit_mean_f_lik = likelihood.get_logit_mean_f_lik(proposed_logit_mean_f, in_state.gamma) / temperature;
		double proposed_logit_mean_g_lik = likelihood.get_logit_mean_g_lik(proposed_logit_mean_g, in_state.gamma) / temperature;

		double jump_back_lik = log((double)phi_jumps_each[in_state.phi[to_change]]);
		double jump_here_lik = log((double)phi_jumps_each[phi_proposed[to_change]]);

		pair<NumericVector, NumericVector> prop_s = get_each_way_sim(
			row_is_column_anc,
			ttsm,
			phi_proposed,
			d.h
		);

		NumericVector prop_x = transform_each_way_sim(
			prop_s,
			proposed_logit_mean_f,
			in_state.log_alpha_plus_beta_f,
			proposed_logit_mean_g,
			in_state.log_alpha_plus_beta_g
		);

		double proposed_y_lik = likelihood.get_y_lik(
			d.y, 
			prop_x, 
			d.g, 
			in_state.alpha_star, 
			in_state.alpha, 
			in_state.log_beta, 
			in_state.gamma
		);

		double numerator = 
			proposed_phi_lik + 
			proposed_logit_mean_f_lik + 
			proposed_logit_mean_g_lik + 
			proposed_y_lik + 
			jump_back_lik
		;

		double denominator = 
			in_state.cur_logit_mean_f_lik + 
			in_state.cur_logit_mean_g_lik + 
			in_state.cur_phi_lik + 
			in_state.cur_y_lik + 
			jump_here_lik
		;

		if (log(unif_rand()) < (numerator - denominator)) {
			in_state.phi = phi_proposed;
			in_state.logit_mean_f = proposed_logit_mean_f;
			in_state.logit_mean_g = proposed_logit_mean_g;
			in_state.cur_phi_lik = proposed_phi_lik;
			in_state.cur_y_lik = proposed_y_lik;
			in_state.cur_logit_mean_f_lik = proposed_logit_mean_f_lik;
			in_state.cur_logit_mean_g_lik = proposed_logit_mean_g_lik;
			in_state._s = prop_s;
			in_state._x = prop_x;
		}
	} else { 
		for (int k = 0; k < in_state.phi.length(); k++)
			in_state.phi[k] = likelihood.pseudo_phi_marginal_prior[random_integer(likelihood.pseudo_phi_marginal_prior.length())];
			in_state.logit_mean_f = norm_rand() * likelihood.pseudo_logit_mean_f_sd + likelihood.pseudo_logit_mean_f_mean;
			in_state.logit_mean_g = norm_rand() * likelihood.pseudo_logit_mean_g_sd + likelihood.pseudo_logit_mean_g_mean;
			in_state.cur_phi_lik = in_state.get_phi_lik(likelihood);
			in_state.cur_logit_mean_f_lik = in_state.get_logit_mean_f_lik(likelihood);
			in_state.cur_logit_mean_g_lik = in_state.get_logit_mean_g_lik(likelihood);

			in_state._s = in_state.get_s(ttsm, row_is_column_anc, d);
			in_state._x = in_state.get_x(ttsm, row_is_column_anc, d);

	}
}

void Update::update_phi(State& in_state, Likelihood& likelihood, double temperature, Data& d, LogicalMatrix row_is_column_anc, NumericMatrix ttsm) {
	if (in_state.gamma) {	
		IntegerVector phi_proposed(in_state.phi.length());

		for (int k = 0; k < in_state.phi.length(); k++)
			phi_proposed[k] = in_state.phi[k];
		
		int to_change = random_integer(in_state.phi.length());

		phi_proposed[to_change] = phi_jumps[random_integer(phi_jumps.length())];

		in_state.phi_proposed = phi_proposed;

		double proposed_phi_lik = likelihood.get_phi_lik(phi_proposed, in_state.gamma) / temperature;
		
		double jump_back_lik = log((double)phi_jumps_each[in_state.phi[to_change]]);
		double jump_here_lik = log((double)phi_jumps_each[phi_proposed[to_change]]);

		pair<NumericVector, NumericVector> prop_s = get_each_way_sim(
			row_is_column_anc,
			ttsm,
			phi_proposed,
			d.h
		);

		NumericVector prop_x = transform_each_way_sim(
			prop_s,
			in_state.logit_mean_f,
			in_state.log_alpha_plus_beta_f,
			in_state.logit_mean_g,
			in_state.log_alpha_plus_beta_g
		);

		double proposed_y_lik = likelihood.get_y_lik(
			d.y, 
			prop_x, 
			d.g, 
			in_state.alpha_star, 
			in_state.alpha, 
			in_state.log_beta, 
			in_state.gamma
		);

		double numerator = 
			proposed_phi_lik + 
			proposed_y_lik + 
			jump_back_lik
		;

		double denominator = 
			in_state.cur_phi_lik + 
			in_state.cur_y_lik + 
			jump_here_lik
		;

		if (log(unif_rand()) < (numerator - denominator)) {
			in_state.phi = phi_proposed;
			in_state.cur_phi_lik = proposed_phi_lik;
			in_state.cur_y_lik = proposed_y_lik;
			in_state._s = prop_s;
			in_state._x = prop_x;
		}
	} else { 
		for (int k = 0; k < in_state.phi.length(); k++) {
			in_state.phi[k] = likelihood.pseudo_phi_marginal_prior[random_integer(likelihood.pseudo_phi_marginal_prior.length())];
			in_state.phi_proposed[k] = in_state.phi[k];
		}
		in_state.cur_phi_lik = in_state.get_phi_lik(likelihood);

		in_state._s = in_state.get_s(ttsm, row_is_column_anc, d);
		in_state._x = in_state.get_x(ttsm, row_is_column_anc, d);
	}
}

void Update::update_logit_mean_f(State& in_state, Likelihood& likelihood, double temperature, Data& d, LogicalMatrix row_is_column_anc, NumericMatrix ttsm) {
	if (in_state.gamma) {
		double proposed_logit_mean_f = in_state.logit_mean_f + norm_rand() * logit_mean_f_proposal_sd;
		
		double proposed_logit_mean_f_lik = likelihood.get_logit_mean_f_lik(proposed_logit_mean_f, in_state.gamma) / temperature;
		
		NumericVector x_proposed = transform_one_way_sim(
			in_state._s.first,
			proposed_logit_mean_f,
			in_state.log_alpha_plus_beta_f
		) * in_state._s.second;
		
		double proposed_y_lik = likelihood.get_y_lik(d.y, x_proposed, d.g, in_state.alpha_star, in_state.alpha, in_state.log_beta, in_state.gamma) / temperature;

		if (log(unif_rand()) < ((proposed_logit_mean_f_lik + proposed_y_lik) - (in_state.cur_logit_mean_f_lik + in_state.cur_y_lik))) {
			in_state.logit_mean_f = proposed_logit_mean_f;
			in_state.cur_y_lik = proposed_y_lik;
			in_state.cur_logit_mean_f_lik = proposed_logit_mean_f_lik;
			in_state._x = x_proposed;
		}
	} else {
		in_state.logit_mean_f = norm_rand() * likelihood.pseudo_logit_mean_f_sd + likelihood.pseudo_logit_mean_f_mean;
		in_state.cur_logit_mean_f_lik = in_state.get_logit_mean_f_lik(likelihood);
		in_state._x = in_state.get_x(ttsm, row_is_column_anc, d);
	}
}

void Update::update_log_alpha_plus_beta_f(State& in_state, Likelihood& likelihood, double temperature, Data& d, LogicalMatrix row_is_column_anc, NumericMatrix ttsm) {
	if (in_state.gamma) {
		double proposed_log_alpha_plus_beta_f = in_state.log_alpha_plus_beta_f + norm_rand() * log_alpha_plus_beta_f_proposal_sd;
		
		double proposed_log_alpha_plus_beta_f_lik = likelihood.get_log_alpha_plus_beta_f_lik(proposed_log_alpha_plus_beta_f, in_state.gamma) / temperature;
		
		NumericVector x_proposed = transform_one_way_sim(
			in_state._s.first,
			in_state.logit_mean_f,
			proposed_log_alpha_plus_beta_f 
		) * in_state._s.second;
		
		double proposed_y_lik = likelihood.get_y_lik(d.y, x_proposed, d.g, in_state.alpha_star, in_state.alpha, in_state.log_beta, in_state.gamma) / temperature;

		if (log(unif_rand()) < ((proposed_log_alpha_plus_beta_f_lik + proposed_y_lik) - (in_state.cur_log_alpha_plus_beta_f_lik + in_state.cur_y_lik))) {
			in_state.log_alpha_plus_beta_f = proposed_log_alpha_plus_beta_f;
			in_state.cur_y_lik = proposed_y_lik;
			in_state.cur_log_alpha_plus_beta_f_lik = proposed_log_alpha_plus_beta_f_lik;
			in_state._x = x_proposed;
		}
	} else {
		in_state.log_alpha_plus_beta_f = norm_rand() * likelihood.pseudo_log_alpha_plus_beta_f_sd + likelihood.pseudo_log_alpha_plus_beta_f_mean;
		in_state.cur_log_alpha_plus_beta_f_lik = in_state.get_log_alpha_plus_beta_f_lik(likelihood);
		in_state._x = in_state.get_x(ttsm, row_is_column_anc, d);
	}
}

void Update::update_logit_mean_g(State& in_state, Likelihood& likelihood, double temperature, Data& d, LogicalMatrix row_is_column_anc, NumericMatrix ttsm) {
	if (in_state.gamma) {
		double proposed_logit_mean_g = in_state.logit_mean_g + norm_rand() * logit_mean_g_proposal_sd;
		
		double proposed_logit_mean_g_lik = likelihood.get_logit_mean_g_lik(proposed_logit_mean_g, in_state.gamma) / temperature;
		
		NumericVector x_proposed = transform_one_way_sim(
			in_state._s.second,
			proposed_logit_mean_g,
			in_state.log_alpha_plus_beta_g
		) * in_state._s.first;
		
		double proposed_y_lik = likelihood.get_y_lik(d.y, x_proposed, d.g, in_state.alpha_star, in_state.alpha, in_state.log_beta, in_state.gamma) / temperature;

		if (log(unif_rand()) < ((proposed_logit_mean_g_lik + proposed_y_lik) - (in_state.cur_logit_mean_g_lik + in_state.cur_y_lik))) {
			in_state.logit_mean_g = proposed_logit_mean_g;
			in_state.cur_y_lik = proposed_y_lik;
			in_state.cur_logit_mean_g_lik = proposed_logit_mean_g_lik;
			in_state._x = x_proposed;
		}
	} else {
		in_state.logit_mean_g = norm_rand() * likelihood.pseudo_logit_mean_g_sd + likelihood.pseudo_logit_mean_g_mean;
		in_state.cur_logit_mean_g_lik = in_state.get_logit_mean_g_lik(likelihood);
		in_state._x = in_state.get_x(ttsm, row_is_column_anc, d);
	}
}

void Update::update_log_alpha_plus_beta_g(State& in_state, Likelihood& likelihood, double temperature, Data& d, LogicalMatrix row_is_column_anc, NumericMatrix ttsm) {
	if (in_state.gamma) {
		double proposed_log_alpha_plus_beta_g = in_state.log_alpha_plus_beta_g + norm_rand() * log_alpha_plus_beta_g_proposal_sd;
		
		double proposed_log_alpha_plus_beta_g_lik = likelihood.get_log_alpha_plus_beta_g_lik(proposed_log_alpha_plus_beta_g, in_state.gamma) / temperature;
		
		NumericVector x_proposed = transform_one_way_sim(
			in_state._s.second,
			in_state.logit_mean_g,
			proposed_log_alpha_plus_beta_g 
		) * in_state._s.first;
		
		double proposed_y_lik = likelihood.get_y_lik(d.y, x_proposed, d.g, in_state.alpha_star, in_state.alpha, in_state.log_beta, in_state.gamma) / temperature;

		if (log(unif_rand()) < ((proposed_log_alpha_plus_beta_g_lik + proposed_y_lik) - (in_state.cur_log_alpha_plus_beta_g_lik + in_state.cur_y_lik))) {
			in_state.log_alpha_plus_beta_g = proposed_log_alpha_plus_beta_g;
			in_state.cur_y_lik = proposed_y_lik;
			in_state.cur_log_alpha_plus_beta_g_lik = proposed_log_alpha_plus_beta_g_lik;
			in_state._x = x_proposed;
		}
	} else {
		in_state.log_alpha_plus_beta_g = norm_rand() * likelihood.pseudo_log_alpha_plus_beta_g_sd + likelihood.pseudo_log_alpha_plus_beta_g_mean;
		in_state.cur_log_alpha_plus_beta_g_lik = in_state.get_log_alpha_plus_beta_g_lik(likelihood);
		in_state._x = in_state.get_x(ttsm, row_is_column_anc, d);
	}
}

void Update::update(State& in_state, Likelihood& likelihood, double temperature, Data& d, LogicalMatrix row_is_column_anc, NumericMatrix ttsm, bool fix_phi) {
	update_gamma(in_state, likelihood, temperature, d, row_is_column_anc, ttsm);
	update_alpha_star(in_state, likelihood, temperature, d, row_is_column_anc, ttsm);
	update_alpha(in_state, likelihood, temperature, d, row_is_column_anc, ttsm);
	update_log_beta(in_state, likelihood, temperature, d, row_is_column_anc, ttsm);
	update_logit_mean_f(in_state, likelihood, temperature, d, row_is_column_anc, ttsm);
	update_log_alpha_plus_beta_f(in_state, likelihood, temperature, d, row_is_column_anc, ttsm);
	update_logit_mean_g(in_state, likelihood, temperature, d, row_is_column_anc, ttsm);
	update_log_alpha_plus_beta_g(in_state, likelihood, temperature, d, row_is_column_anc, ttsm);

	if (joint_proposal) {
		update_phi_and_inflexion(in_state, likelihood, temperature, d, row_is_column_anc, ttsm);
	}
	else {
		update_logit_mean_g(in_state, likelihood, temperature, d, row_is_column_anc, ttsm);
		update_logit_mean_f(in_state, likelihood, temperature, d, row_is_column_anc, ttsm);
		if (!fix_phi) update_phi(in_state, likelihood, temperature, d, row_is_column_anc, ttsm);
	}
}


