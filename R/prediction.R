fg_step_pp <- function(
	s_phi_new,
	s_x_new,
	s_phi,
	s_x,
	y,
	x_values=seq(from=0, to=1, length.out=30),
	tf_s1=fg_step_defaults$tf_s1,
	tf_s2=fg_step_defaults$tf_s2,
	tg_s1=fg_step_defaults$tg_s1,
	tg_s2=fg_step_defaults$tg_s2,
	q_s1=mean(y),
	q_s2=fg_step_defaults$q_s2,
	p_s1=fg_step_defaults$p_s1,
	p_s2=fg_step_defaults$p_s2
) {
	tf_probs <- discrete_beta(x_values, tf_s1, tf_s2)
	tg_probs <- discrete_beta(x_values, tg_s1, tg_s2)

	x1y1 <- sims_to_sumgrid(s_phi[y], s_x[y], x_values)
	x1y0 <- sims_to_sumgrid(s_phi[!y], s_x[!y], x_values)

	x0y1 <- sum(y)-x1y1
	x0y0 <- sum(!y)-x1y0
	liks <- outer(FUN="+", log(tf_probs), log(tg_probs)) + lbeta(x1y1 + p_s1, x1y0 + p_s2) - lbeta(p_s1, p_s2) + lbeta(x1y1[1]-x1y1 + q_s1, x1y0[1]-x1y0 + q_s2) - lbeta(q_s1, q_s2)
	t_probs <- exp(liks)/sum(exp(liks))
	eq_over_t <- (q_s1 + x0y1) / (q_s1 + x0y1 + q_s2 + x0y0)
	ep_over_t <- (p_s1 + x1y1) / (p_s1 + x1y1 + p_s2 + x1y0)

	p_grd <- cumsumgrid(ep_over_t * t_probs)
	q_grd <- cumsumgrid(eq_over_t * t_probs)
	neg_q_grd <- max(q_grd) - q_grd

	levels_coords <- cbind(as.integer(cut(s_phi_new, breaks=c(-Inf, x_values))), as.integer(cut(s_x_new, breaks=c(-Inf, x_values))))
	
	p_grd[levels_coords]+neg_q_grd[levels_coords]
}

fg_step_tab_pp <- function(
	N, 
	q_s1=fg_step_defaults$q_s1,
	q_s2=fg_step_defaults$q_s2,
	p_s1=fg_step_defaults$p_s1,
	p_s2=fg_step_defaults$p_s2
) {
	q_s1_gt <- c(0, cumsum(log(seq(from=0, to=N-1) + q_s1)))
	q_s2_gt <- c(0, cumsum(log(seq(from=0, to=N-1) + q_s2)))
	q_t_gt <- c(0, cumsum(log(seq(from=0, to=N-1) + q_s1 + q_s2)))
	p_s1_gt <- c(0, cumsum(log(seq(from=0, to=N-1) + p_s1)))
	p_s2_gt <- c(0, cumsum(log(seq(from=0, to=N-1) + p_s2)))
	p_t_gt <- c(0, cumsum(log(seq(from=0, to=N-1) + p_s1 + p_s2)))

	function(
		s_phi_new,
		s_x_new,
		s_phi,
		s_x,
		y,
		min_log_ML=-Inf,
		x_values=30,
		tf_s1=fg_step_defaults$tf_s1,
		tf_s2=fg_step_defaults$tf_s2,
		tg_s1=fg_step_defaults$tg_s1,
		tg_s2=fg_step_defaults$tg_s2
	) {
		xs <- seq(from=0, to=1, length.out=x_values)
		tf_probs <- discrete_beta(xs, tf_s1, tf_s2)
		tg_probs <- discrete_beta(xs, tg_s1, tg_s2)

		x1y1 <- sumgrid(s_phi[y], s_x[y], x_values)
		x1y0 <- sumgrid(s_phi[!y], s_x[!y], x_values)

		liks <- outer(FUN="+", log(tf_probs), log(tg_probs)) + p_s1_gt[x1y1+1L] + p_s2_gt[x1y0+1L] - p_t_gt[x1y0 + x1y1 + 1L] + q_s1_gt[x1y1[1]-x1y1 + 1L] + q_s2_gt[x1y0[1]-x1y0 + 1L] - q_t_gt[x1y1[1]-x1y1 + x1y0[1]-x1y0 + 1L]

		x0y1 <- sum(y)-x1y1
		x0y0 <- sum(!y)-x1y0

		t_probs <- exp(liks)/sum(exp(liks))
		eq_over_t <- (q_s1 + x0y1) / (q_s1 + x0y1 + q_s2 + x0y0)
		ep_over_t <- (p_s1 + x1y1) / (p_s1 + x1y1 + p_s2 + x1y0)

		p_grd <- cumsumgrid(ep_over_t * t_probs)
		q_grd <- cumsumgrid(eq_over_t * t_probs)
		neg_q_grd <- max(q_grd) - q_grd

		levels_coords <- cbind(as.integer(cut(s_phi_new, breaks=c(-Inf, xs))), as.integer(cut(s_x_new, breaks=c(-Inf, xs))))
		
		p_grd[levels_coords]+neg_q_grd[levels_coords]
	}
}



f_step_pp <- function(
	x_new,
	x,
	y,
	x_values=seq(from=0, to=1, length.out=30),
	t_s1=f_step_defaults$t_s1,
	t_s2=f_step_defaults$t_s2,
	q_s1=mean(y),
	q_s2=f_step_defaults$q_s2,
	p_s1=f_step_defaults$p_s1,
	p_s2=f_step_defaults$p_s2
) {
	t_prior_probs <- local({ h0 <- dbeta(x_values[-1], shape1=t_s1, shape2=t_s2); h <- c(h0[1], h0); h/sum(h) })
	x1y1 <- rev(cumsum(rev(as.integer(table(cut(x[y], breaks=c(-Inf, x_values)))))))
	x1y0 <- rev(cumsum(rev(as.integer(table(cut(x[!y], breaks=c(-Inf, x_values)))))))

	x0y1 <- sum(y)-x1y1
	x0y0 <- sum(!y)-x1y0

	liks <- log(t_probs) + lbeta(x1y1 + p_s1, x1y0 + p_s2) - lbeta(p_s1, p_s2) + lbeta(x1y1[1]-x1y1 + q_s1, x1y0[1]-x1y0 + q_s2) - lbeta(q_s1, q_s2)

	t_probs <- exp(liks)/sum(exp(liks))
	eq_over_t <- (q_s1 + x0y1) / (q_s1 + x0y1 + q_s2 + x0y0)
	ep_over_t <- (p_s1 + x1y1) / (p_s1 + x1y1 + p_s2 + x1y0)

	p_vec <- cumsum(ep_over_t * t_probs)
	q_vec <- cumsum(eq_over_t * t_probs)
	neg_q_vec <- max(q_vec) - q_vec

	levels_coords <- as.integer(cut(x_new, breaks=c(-Inf, x_values)))
	p_vec[levels_coords]+neg_q_vec[levels_coords]
}

#' Predicted probability of \code{y} given \code{x} conditional on association and given data.
#'
#' @template ontology
#' @template x
#' @template y
#' @template sim_reg_out
#' @param x_new New \code{list} of ontological term sets to perform prediction on. Defaults to \code{x}.
#' @template information_content
#' @template sim_params
#' @template two_way
#' @param prediction_fn Function for computing predicted probabilities for \code{y[i]=TRUE}.
#' @param min_ratio Threshold for fraction of posterior probability which sampled phi must hold in order to be included in sum. 
#' @param ... Additional arguments to pass to \code{prediction_fn}.
#' @return Vector of predicted probabilities corresponding to term sets in \code{x_new}.
#' @export
posterior_prediction <- function(
	ontology,
	x,
	y,
	sim_reg_out,
	x_new=x,
	information_content=get_term_info_content(ontology, x),
	sim_params=list(ontology=ontology, information_content=information_content),
	two_way=TRUE,
	prediction_fn=SimReg:::fg_step_tab_pp(N=length(y)),
	min_ratio=1e-3,
	...
) {
	use_phis <- which(sim_reg_out$ml >= max(sim_reg_out$ml) + log(min_ratio))
	
	mls <- exp(sim_reg_out$ml[use_phis])
	phis <- sim_reg_out$phis[use_phis]

	sims <- if (two_way) {
		s_phis <- split(do.call(what=get_asym_sim_grid, c(sim_params, list(B=x, A=phis))), seq(length(phis)))
		s_xs <- split(t(do.call(what=get_asym_sim_grid, c(sim_params, list(A=x, B=phis)))), seq(length(phis)))
		s_phis_new <- split(do.call(what=get_asym_sim_grid, c(sim_params, list(B=x_new, A=phis))), seq(length(phis)))
		s_xs_new <- split(t(do.call(what=get_asym_sim_grid, c(sim_params, list(A=x_new, B=phis)))), seq(length(phis)))
		mapply(SIMPLIFY=FALSE, FUN=list, s_phi=s_phis, s_x=s_xs, s_phi_new=s_phis_new, s_x_new=s_xs_new)
	} else {
		xs <- split(do.call(what=get_asym_sim_grid, c(sim_params, list(B=x, A=phis))), seq(length(phis)))
		xs_new <- split(do.call(what=get_asym_sim_grid, c(sim_params, list(B=x_new, A=phis))), seq(length(phis)))
		mapply(SIMPLIFY=FALSE, FUN=list, x=xs, x_new=xs_new)
	}

	apply(mapply(
		FUN=function(pr_cond_on_phi, pr_phi) pr_phi * pr_cond_on_phi,
		lapply(sims, function(sim) do.call(what=prediction_fn, c(sim, list(y=y, ...)))),
		mls
	), 1, sum)/sum(mls)
}
