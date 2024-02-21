#' Calculate sum of log probabilities on log scale without over/under-flow
#'
#' @param log_probs Numeric vector of probabilities on log scale.
#' @return Numeric value on log scale.
#' @export
sum_log_probs <- function(log_probs) log(sum(exp(log_probs-max(log_probs))))+max(log_probs)

cumsumgrid <- function(x) apply(apply(x, 1, cumsum), 1, cumsum)

#' @importFrom stats dbeta
discrete_beta <- function(x, s1, s2) {
	h0 <- dbeta(x[-1], shape1=s1, shape2=s2)
	h <- c(h0[1], h0)
	h/sum(h)
}

sims_to_sumgrid <- function(s_phi, s_x, x_values) {
	tab <- table(
		data.frame(
			phi=cut(s_phi, breaks=c(-Inf, x_values)),
			x=cut(s_x, breaks=c(-Inf, x_values))))

	sq <- seq(from=length(x_values), to=1, by=-1)
	cumsumgrid(tab[sq,sq])[sq,sq]
}

logit <- function(x) log(x)-log(1-x)
expit <- function(x) exp(x)/(1+exp(x))

fg_step_defaults <- list(
	tf_s1=3,
	tf_s2=1,
	tg_s1=1,
	tg_s2=1,
	q_s1=0.01,
	q_s2=1,
	p_s1=1,
	p_s2=1
)

fg_step_tab <- function(
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
		tf_probs <- discrete_beta(seq(from=0, to=1, length.out=x_values), tf_s1, tf_s2)
		tg_probs <- discrete_beta(seq(from=0, to=1, length.out=x_values), tg_s1, tg_s2)

		x1y1 <- sumgrid(s_phi[y], s_x[y], x_values)
		x1y0 <- sumgrid(s_phi[!y], s_x[!y], x_values)

		liks <- outer(FUN="+", log(tf_probs), log(tg_probs)) + p_s1_gt[x1y1+1L] + p_s2_gt[x1y0+1L] - p_t_gt[x1y0 + x1y1 + 1L] + q_s1_gt[x1y1[1]-x1y1 + 1L] + q_s2_gt[x1y0[1]-x1y0 + 1L] - q_t_gt[x1y1[1]-x1y1 + x1y0[1]-x1y0 + 1L]
		sum_log_probs(liks)
	}
}

fg_step <- function(
	s_phi,
	s_x,
	y,
	min_log_ML=-Inf,
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

	liks <- outer(FUN="+", log(tf_probs), log(tg_probs)) + lbeta(x1y1 + p_s1, x1y0 + p_s2) - lbeta(p_s1, p_s2) + lbeta(x1y1[1]-x1y1 + q_s1, x1y0[1]-x1y0 + q_s2) - lbeta(q_s1, q_s2)
	sum_log_probs(liks)
}

#' @importFrom stats dnorm pbeta
full_fg_lik <- function(
	params,
	y,
	s_phi,
	s_x,
	alpha_mean=0,
	alpha_sd=1,
	log_beta_mean=1, 
	log_beta_sd=1, 
	logit_f_mean=1.5, 
	logit_f_sd=0.5, 
	log_f_a_plus_b_mean=1, 
	log_f_a_plus_b_sd=0.5, 
	logit_g_mean=0, 
	logit_g_sd=1, 
	log_g_a_plus_b_mean=1, 
	log_g_a_plus_b_sd=1
) {
	alpha <- params[1]
	beta <- exp(params[2])
	f_mean <- expit(params[3])
	f_a_plus_b <- exp(params[4])
	g_mean <- expit(params[5])
	g_a_plus_b <- exp(params[6])

	f_shape1 <- f_a_plus_b * f_mean
	f_shape2 <- f_a_plus_b * (1.0 - f_mean)
	g_shape1 <- g_a_plus_b * g_mean
	g_shape2 <- g_a_plus_b * (1.0 - g_mean)

	x <- transform(s_phi, f_shape1, f_shape2)*transform(s_x, g_shape1, g_shape2)
	
	p <- expit(alpha + beta * x)
	sum(
		ifelse(y, log(p), log(1-p)),
		dnorm(log=TRUE, x=alpha, mean=alpha_mean, sd=alpha_sd),
		dnorm(log=TRUE, x=log(beta), mean=log_beta_mean, sd=log_beta_sd),
		dnorm(log=TRUE, x=logit(f_mean), mean=logit_f_mean, sd=logit_f_sd),
		dnorm(log=TRUE, x=logit(g_mean), mean=logit_g_mean, sd=logit_g_sd),
		dnorm(log=TRUE, x=log(f_a_plus_b), mean=log_f_a_plus_b_mean, sd=log_f_a_plus_b_sd),
		dnorm(log=TRUE, x=log(g_a_plus_b), mean=log_g_a_plus_b_mean, sd=log_g_a_plus_b_sd)
	)
}

#' @importFrom stats optim
fg_lap <- function(
	s_phi,
	s_x,
	y,
	min_log_ML=-Inf,
	alpha_mean=log(mean(y)), 
	alpha_sd=1, 
	log_beta_mean=1, 
	log_beta_sd=1, 
	logit_f_mean=1.5, 
	logit_f_sd=0.5, 
	log_f_a_plus_b_mean=1, 
	log_f_a_plus_b_sd=0.5, 
	logit_g_mean=0, 
	logit_g_sd=1, 
	log_g_a_plus_b_mean=1, 
	log_g_a_plus_b_sd=1
) {
	opt <- optim(
		fn=full_fg_lik, 
		par=c(
			alpha_mean,
			log_beta_mean,
			logit_f_mean,
			log_f_a_plus_b_mean,
			logit_g_mean,
			log_g_a_plus_b_mean
		),
		y=y, 
		s_x=s_x, 
		s_phi=s_phi,
		alpha_mean=alpha_mean, 
		alpha_sd=alpha_sd, 
		log_beta_mean=log_beta_mean,
		log_beta_sd=log_beta_sd,
		logit_f_mean=logit_f_mean,
		logit_f_sd=logit_f_sd,
		log_f_a_plus_b_mean=log_f_a_plus_b_mean,
		log_f_a_plus_b_sd=log_f_a_plus_b_sd,
		logit_g_mean=logit_g_mean,
		logit_g_sd=logit_g_sd,
		log_g_a_plus_b_mean=log_g_a_plus_b_mean,
		log_g_a_plus_b_sd=log_g_a_plus_b_sd,
		hessian=TRUE,
		control=list(fnscale=-1)
	)

	(log(2) + log(pi)) * 3 + 0.5 * log(det(-solve(opt$hessian))) + opt$value
}

f_step_defaults <- list(
	t_s1=3,
	t_s2=1,
	q_s1=0.01,
	q_s2=1,
	p_s1=3,
	p_s2=1
)

f_step <- function(
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
	t_probs <- local({ h0 <- dbeta(x_values[-1], shape1=t_s1, shape2=t_s2); h <- c(h0[1], h0); h/sum(h) })
	y1 <- rev(cumsum(rev(as.integer(table(cut(x[y], breaks=c(-Inf, x_values)))))))
	y0 <- rev(cumsum(rev(as.integer(table(cut(x[!y], breaks=c(-Inf, x_values)))))))
	liks <- log(t_probs) + lbeta(y1 + p_s1, y0 + p_s2) - lbeta(p_s1, p_s2) + lbeta(y1[1]-y1 + q_s1, y0[1]-y0 + q_s2) - lbeta(q_s1, q_s2)
	sum_log_probs(liks)
}

#' @importFrom Rcpp evalCpp
#' @useDynLib SimReg 
fg <- function(
	s_phi, 
	s_x, 
	y, 
	x_values=seq(from=0, to=1, length.out=20), 
	t=(0:100/100)^2, 
	max_samples=30, 
	min_samples=7, 
	min_log_ML=-Inf,
	log_scale_tolerance=-2.3,
	alpha_mean=log(mean(y)), 
	alpha_sd=1, 
	log_beta_mean=1, 
	log_beta_sd=1, 
	logit_f_mean=1.5, 
	logit_f_sd=0.5, 
	log_f_a_plus_b_mean=1, 
	log_f_a_plus_b_sd=0.5, 
	logit_g_mean=0, 
	logit_g_sd=1, 
	log_g_a_plus_b_mean=1, 
	log_g_a_plus_b_sd=1, 
	alpha_prop_sd=0.8, 
	log_beta_prop_sd=0.8, 
	logit_f_mean_prop_sd=0.8, 
	log_f_a_plus_b_prop_sd=0.8, 
	logit_g_mean_prop_sd=0.8, 
	log_g_a_plus_b_prop_sd=0.8
) { 

	y1_tab <- table(
		data.frame(
			phi=cut(s_phi[y], breaks=c(-Inf, x_values)),
			x=cut(s_x[y], breaks=c(-Inf, x_values))))
	y0_tab <- table(
		data.frame(
			phi=cut(s_phi[!y], breaks=c(-Inf, x_values)),
			x=cut(s_x[!y], breaks=c(-Inf, x_values))))

	tab <- y1_tab + y0_tab

	who <- which(arr.ind=TRUE, tab > 0)

	ML(
		s_phi_values=x_values[who[,1]],
		s_x_values=x_values[who[,2]],
		num_y0_phi=y0_tab[who],
		num_y1_phi=y1_tab[who],
		t=t,
		log_scale_tolerance=log_scale_tolerance,
		max_samples=max_samples,
		min_samples=min_samples,
		min_log_ML=min_log_ML,
		alpha_mean=alpha_mean,
		alpha_sd=alpha_sd,
		log_beta_mean=log_beta_mean,
		log_beta_sd=log_beta_sd,
		logit_f_mean=logit_f_mean,
		logit_f_sd=logit_f_sd,
		log_f_a_plus_b_mean=log_f_a_plus_b_mean,
		log_f_a_plus_b_sd=log_f_a_plus_b_sd,
		logit_g_mean=logit_g_mean,
		logit_g_sd=logit_g_sd,
		log_g_a_plus_b_mean=log_g_a_plus_b_mean,
		log_g_a_plus_b_sd=log_g_a_plus_b_sd,
		alpha_prop_sd=alpha_prop_sd,
		log_beta_prop_sd=log_beta_prop_sd,
		logit_f_mean_prop_sd=logit_f_mean_prop_sd,
		log_f_a_plus_b_prop_sd=log_f_a_plus_b_prop_sd,
		logit_g_mean_prop_sd=logit_g_mean_prop_sd,
		log_g_a_plus_b_prop_sd=log_g_a_plus_b_prop_sd
	)
}

#' @useDynLib SimReg 
f <- function(
	x, 
	y, 
	x_values=seq(from=0, to=1, length.out=20), 
	t=(0:100/100)^2, 
	max_samples=50, 
	min_samples=10, 
	min_log_ML=-Inf,
	log_scale_tolerance=-2.3,
	alpha_mean=log(mean(y)), 
	alpha_sd=1, 
	log_beta_mean=1, 
	log_beta_sd=1, 
	logit_f_mean=1.5, 
	logit_f_sd=0.5, 
	log_f_a_plus_b_mean=1, 
	log_f_a_plus_b_sd=0.5, 
	alpha_prop_sd=0.8, 
	log_beta_prop_sd=0.8, 
	logit_f_mean_prop_sd=0.8, 
	log_f_a_plus_b_prop_sd=0.8
) { 
	y1 <- as.integer(table(cut(x[y], breaks=c(-Inf, x_values))))
	y0 <- as.integer(table(cut(x[!y], breaks=c(-Inf, x_values))))
	f_ML(
		x_values=x_values,
		num_y0_phi=y0,
		num_y1_phi=y1,
		t=t,
		log_scale_tolerance=log_scale_tolerance,
		max_samples=max_samples,
		min_samples=min_samples,
		min_log_ML=min_log_ML,
		alpha_mean=alpha_mean,
		alpha_sd=alpha_sd,
		log_beta_mean=log_beta_mean,
		log_beta_sd=log_beta_sd,
		logit_f_mean=logit_f_mean,
		logit_f_sd=logit_f_sd,
		log_f_a_plus_b_mean=log_f_a_plus_b_mean,
		log_f_a_plus_b_sd=log_f_a_plus_b_sd,
		alpha_prop_sd=alpha_prop_sd,
		log_beta_prop_sd=log_beta_prop_sd,
		logit_f_mean_prop_sd=logit_f_mean_prop_sd,
		log_f_a_plus_b_prop_sd=log_f_a_plus_b_prop_sd
	)
}

bg_rate <- function(y, q_s1=1, q_s2=1) {
	liks <- lbeta(sum(y) + q_s1, sum(!y) + q_s2) - lbeta(q_s1, q_s2)
	log(sum(exp(liks-min(liks))))+min(liks)
}

#' @useDynLib SimReg 
bg <- function(y, alpha_mean=log(mean(y)), alpha_sd=1, t=(0:100/100)^2, n_samples=300, alpha_prop_sd=1) bg_ML(y0=sum(!y), y1=sum(y), t=t, n_samples=n_samples, alpha_mean=alpha_mean, alpha_sd=alpha_sd, alpha_prop_sd=alpha_prop_sd)
