.sim.reg <- function(
	hpo.terms,
	y,
	x=NULL, 
	g=rep(0, length(y)),
	its=10000,
	thin=1,
	burn=2000,
	record_sims=FALSE,

	information.content=get.term.info.content(hpo.terms, hpo.phenotypes=x),
	ttsm=get.term.term.matrix(hpo.terms, information.content, TRUE),

	row_is_column_anc=get.term.descendancy.matrix(hpo.terms, names(information.content)),
	case_ids=unlist(mapply(SIMPLIFY=FALSE, FUN=rep, 0:(length(x)-1), sapply(x, length))),
	term_ids=as.integer(match(unlist(x), colnames(row_is_column_anc)))-1,

	gamma=(runif(1) < gamma_prior_prob)[1],
	alpha_star=rnorm(n=1, mean=alpha_star_mean, sd=alpha_star_sd),
	alpha=rnorm(n=1, mean=alpha_mean, sd=alpha_sd),
	log_beta=rnorm(n=1, mean=log_beta_mean, sd=log_beta_sd),
	phi=sample.int(n=ncol(row_is_column_anc), size=3, replace=TRUE)-1,
	logit_mean_f=rnorm(n=1, mean=logit_mean_f_mean, sd=logit_mean_f_sd),
	log_alpha_plus_beta_f=rnorm(n=1, mean=log_alpha_plus_beta_f_mean, sd=log_alpha_plus_beta_f_sd),
	logit_mean_g=rnorm(n=1, mean=logit_mean_g_mean, sd=logit_mean_g_sd),
	log_alpha_plus_beta_g=rnorm(n=1, mean=log_alpha_plus_beta_g_mean, sd=log_alpha_plus_beta_g_sd),

	gamma_prior_prob=0.05,
	alpha_star_mean=0,
	alpha_mean=0,
	alpha_star_sd=5,
	alpha_sd=5,
	log_beta_mean=2,
	log_beta_sd=1,
	logit_mean_f_mean=1,
	logit_mean_f_sd=1,
	log_alpha_plus_beta_f_mean=2,
	log_alpha_plus_beta_f_sd=1,
	logit_mean_g_mean=0,
	logit_mean_g_sd=1.5,
	log_alpha_plus_beta_g_mean=2,
	log_alpha_plus_beta_g_sd=1,
	pseudo_alpha_star_mean=0,
	pseudo_alpha_mean=0,
	pseudo_alpha_star_sd=5,
	pseudo_alpha_sd=5,
	pseudo_log_beta_mean=2,
	pseudo_log_beta_sd=2,
	pseudo_logit_mean_f_mean=2,
	pseudo_logit_mean_f_sd=2,
	pseudo_log_alpha_plus_beta_f_mean=2,
	pseudo_log_alpha_plus_beta_f_sd=2,
	pseudo_logit_mean_g_mean=2,
	pseudo_logit_mean_g_sd=2,
	pseudo_log_alpha_plus_beta_g_mean=2,
	pseudo_log_alpha_plus_beta_g_sd=2,
	alpha_star_proposal_sd=2,
	alpha_proposal_sd=2,
	log_beta_proposal_sd=2,
	logit_mean_f_proposal_sd=2,
	log_alpha_plus_beta_f_proposal_sd=2,
	logit_mean_g_proposal_sd=2,
	log_alpha_plus_beta_g_proposal_sd=2,
	phi_jumps=c(0:(ncol(row_is_column_anc)-1), rep(match(unlist(lapply(x[y], get.ancestors, hpo.terms=hpo.terms)), colnames(row_is_column_anc))-1, 50)),
	pseudo_phi_marginal_prior=c(0:(ncol(row_is_column_anc)-1), rep(match(unlist(lapply(x[y], get.ancestors, hpo.terms=hpo.terms)), colnames(row_is_column_anc))-1, 50)),
	phi_num_leaves_geometric_rate=1
) {
	result <- .Call(
		"R_sim_reg",
		PACKAGE="SimReg",
		its,
		thin,
		record_sims,
		ttsm,
		term_ids,
		case_ids,
		y,
		g,

		gamma,
		alpha_star,
		alpha,
		log_beta,
		phi,
		logit_mean_f,
		log_alpha_plus_beta_f,
		logit_mean_g,
		log_alpha_plus_beta_g,

		gamma_prior_prob,
		alpha_star_mean,
		alpha_mean,
		alpha_star_sd,
		alpha_sd,
		log_beta_mean,
		log_beta_sd,
		logit_mean_f_mean,
		logit_mean_f_sd,
		log_alpha_plus_beta_f_mean,
		log_alpha_plus_beta_f_sd,
		logit_mean_g_mean,
		logit_mean_g_sd,
		log_alpha_plus_beta_g_mean,
		log_alpha_plus_beta_g_sd,
		pseudo_alpha_star_mean,
		pseudo_alpha_mean,
		pseudo_alpha_star_sd,
		pseudo_alpha_sd,
		pseudo_log_beta_mean,
		pseudo_log_beta_sd,
		pseudo_logit_mean_f_mean,
		pseudo_logit_mean_f_sd,
		pseudo_log_alpha_plus_beta_f_mean,
		pseudo_log_alpha_plus_beta_f_sd,
		pseudo_logit_mean_g_mean,
		pseudo_logit_mean_g_sd,
		pseudo_log_alpha_plus_beta_g_mean,
		pseudo_log_alpha_plus_beta_g_sd,
		pseudo_phi_marginal_prior,
		alpha_star_proposal_sd,
		alpha_proposal_sd,
		log_beta_proposal_sd,
		logit_mean_f_proposal_sd,
		log_alpha_plus_beta_f_proposal_sd,
		logit_mean_g_proposal_sd,
		log_alpha_plus_beta_g_proposal_sd,
		phi_jumps,

		rep(1, ncol(row_is_column_anc)),
		row_is_column_anc,
		phi_num_leaves_geometric_rate,
		FALSE,
		TRUE,
		FALSE,
		FALSE	
	)


	trunc <- function(y) if (class(y) == "matrix") { if (nrow(y) == its) y[(burn+1):(nrow(y)),,drop=FALSE] else y } else { y[(burn+1):(length(y))] }

	result <- c(
		list(liks=lapply(result$liks, trunc)),
		lapply(result[setdiff(names(result), "liks")], trunc)
	)

	result$phi <- apply(result$phi, 2, function(x) colnames(row_is_column_anc)[x+1])

	c(	
		list(
			priors=list(
				gamma_prior_prob=gamma_prior_prob,
				alpha_star_mean=alpha_star_mean,
				alpha_mean=alpha_mean,
				alpha_star_sd=alpha_star_sd,
				alpha_sd=alpha_sd,
				log_beta_mean=log_beta_mean,
				log_beta_sd=log_beta_sd,
				logit_mean_f_mean=logit_mean_f_mean,
				logit_mean_f_sd=logit_mean_f_sd,
				log_alpha_plus_beta_f_mean=log_alpha_plus_beta_f_mean,
				log_alpha_plus_beta_f_sd=log_alpha_plus_beta_f_sd,
				logit_mean_g_mean=logit_mean_g_mean,
				logit_mean_g_sd=logit_mean_g_sd,
				log_alpha_plus_beta_g_mean=log_alpha_plus_beta_g_mean,
				log_alpha_plus_beta_g_sd=log_alpha_plus_beta_g_sd,
				pseudo_alpha_star_mean=pseudo_alpha_star_mean,
				pseudo_alpha_star_sd=pseudo_alpha_star_sd,
				pseudo_alpha_mean=pseudo_alpha_mean,
				pseudo_alpha_sd=pseudo_alpha_sd,
				pseudo_log_beta_mean=pseudo_log_beta_mean,
				pseudo_log_beta_sd=pseudo_log_beta_sd,
				pseudo_logit_mean_f_mean=pseudo_logit_mean_f_mean,
				pseudo_logit_mean_f_sd=pseudo_logit_mean_f_sd,
				pseudo_log_alpha_plus_beta_f_mean=pseudo_log_alpha_plus_beta_f_mean,
				pseudo_log_alpha_plus_beta_f_sd=pseudo_log_alpha_plus_beta_f_sd,
				pseudo_logit_mean_g_mean=pseudo_logit_mean_g_mean,
				pseudo_logit_mean_g_sd=pseudo_logit_mean_g_sd,
				pseudo_log_alpha_plus_beta_g_mean=pseudo_log_alpha_plus_beta_g_mean,
				pseudo_log_alpha_plus_beta_g_sd=pseudo_log_alpha_plus_beta_g_sd,
				pseudo_phi_marginal_prior=pseudo_phi_marginal_prior,
				alpha_star_proposal_sd=alpha_star_proposal_sd,
				alpha_proposal_sd=alpha_proposal_sd,
				log_beta_proposal_sd=log_beta_proposal_sd,
				logit_mean_f_proposal_sd=logit_mean_f_proposal_sd,
				log_alpha_plus_beta_f_proposal_sd=log_alpha_plus_beta_f_proposal_sd,
				logit_mean_g_proposal_sd=logit_mean_g_proposal_sd,
				log_alpha_plus_beta_g_proposal_sd=log_alpha_plus_beta_g_proposal_sd
			)
		), 
		result, 
		c(
			result,
			c(
				list(phi_accept=sapply(1:(nrow(result$phi)-1), function(x) !identical(result$phi[x,], result$phi[x+1,]))),
				lapply(
					c("alpha_star", "alpha", "log_beta", "logit_mean_f", "log_alpha_plus_beta_f", "logit_mean_g", "log_alpha_plus_beta_g") %>%
					(function(x) setNames(x, paste(x, "_accept", sep=""))),
					function(x) result[[x]][1:(length(result[[x]])-1)] != result[[x]][2:length(result[[x]])]
				)
			)
		)
	)
}

#' Similarity regression
#'
#' Performins Bayesian `similarity regression' on given binary genotype \code{y} (logical vector) against HPO encoded phenotype \code{x} (list of character vectors of HPO term IDs). It returns a list of traces for the various estimated parameters. Of particular interest are the estimated mean posterior of \code{gamma} (the model selection indicator, thus giving an estimate of the probability of an association under the model assumptions - obtained with \code{mean(result$gamma)}) and the posterior distribution of the characteristic phenotype (which can be visualised by the functions \code{\link{hpo.plot.marginal.freqs}}, \code{\link{two.term.marginals.plot}}, and \code{\link{single.term.marginals.plot}}).
#'
#' @template hpo.terms
#' @param y Logical vector of genotypes (typically 1 for rare genotype, 0 for common genotype)
#' @param x List of character vectors of HPO phenotypes of cases
#' @param g Genotype log odds offsets per individual
#' @param its Number of update cycles to perform 
#' @param thin Factor by which to thin resultant chains of parameter samples
#' @param record_sims Logical indicating whether to record trace of similarities
#' @param information.content Numeric vector, named by HPO IDs, containing the information content of corresponding terms
#' @param ttsm The `term-term' similarity matrix, a numeric matrix whose dimensions are named by the terms so that cell i,j contains the similarity of term i to term j
#' @param row_is_column_anc Logical matrix, whose dimensions are named by HPO term IDs so that cell i,j is TRUE if i is an ancestor term of j
#' @param case_ids IDs for the N cases from 0 to N-1, indicating which case terms in \code{term_ids} belong to (automatically determined given x)
#' @param term_ids Vector of HPO term IDs belonging to cases
#' @param return_pseudo_runs Logical indicating whether to return the MCMC output of the tuning phase of the inference procedure
#' @param tuning_its Number of update cycles to perform in the tuning phase of the inference procedure
#' @param tuning_burn Number of update cycles to discard in tuning phase
#' @param burn Number of update cycles to discard 
#' @param gamma Initial value of model selection indicator gamma.
#' @param alpha_star  Initial value of alpha_star, the rate of observing the rare genotype y = 1 under gamma = 0, i.e. the no association model 
#' @param alpha Initial value of alpha, the background rate of observing the rare genotype under gamma = 1
#' @param log_beta Initial value of log_beta, the log of the effect size of onotological similarity
#' @param phi Character vector of HPO term IDs giving the initial value of phi, the characteristic phenotype
#' @param logit_mean_f Initial value of logit_mean_f
#' @param log_alpha_plus_beta_f Initial value of log_alpha_plus_beta_f
#' @param logit_mean_g Initial value of logit_mean_g
#' @param log_alpha_plus_beta_g Initial value of log_alpha_plus_beta_g
#' @param gamma_prior_prob Prior probability of gamma = 1
#' @param alpha_star_mean Prior mean of alpha_star given gamma = 0
#' @param alpha_mean Prior mean of alpha given gamma = 1
#' @param alpha_star_sd Prior sd of alpha_star given gamma = 0
#' @param alpha_sd Prior sd of alpha given gamma = 1
#' @param log_beta_mean Prior mean of log_beta given gamma = 1
#' @param log_beta_sd Prior sd of log_beta given gamma = 1
#' @param logit_mean_f_mean Prior mean of logit_mean_f given gamma = 1
#' @param logit_mean_f_sd Prior sd of logit_mean_f given gamma = 1
#' @param log_alpha_plus_beta_f_mean Prior mean of log_alpha_plus_beta_f given gamma = 1
#' @param log_alpha_plus_beta_f_sd Prior sd of log_alpha_plus_beta_f given gamma = 1
#' @param logit_mean_g_mean Prior mean of logit_mean_g given gamma = 1
#' @param logit_mean_g_sd Prior sd of logit_mean_g given gamma = 1
#' @param log_alpha_plus_beta_g_mean Prior mean of log_alpha_plus_beta_g given gamma = 1
#' @param log_alpha_plus_beta_g_sd Prior sd of log_alpha_plus_beta_g given gamma = 1
#' @param pseudo_phi_marginal_prior Vector of HPO term IDs to be used as prior distribution on marginal probability of single term in phi given gamma = 0
#' @param alpha_star_proposal_sd Proposal sd of local jumps in MH updates of alpha_star used during inference
#' @param alpha_proposal_sd Proposal sd of local jumps in MH updates of alpha used during inference
#' @param log_beta_proposal_sd Proposal sd of local jumps in MH updates of log_beta used during inference
#' @param logit_mean_f_proposal_sd Proposal sd of local jumps in MH updates of logit_mean_f used during inference
#' @param log_alpha_plus_beta_f_proposal_sd Proposal sd of local jumps in MH updates of log_alpha_plus_beta_f used during inference
#' @param logit_mean_g_proposal_sd Proposal sd of local jumps in MH updates of logit_mean_g used during inference
#' @param log_alpha_plus_beta_g_proposal_sd Proposal sd of local jumps in MH updates of log_alpha_plus_beta_g used during inference
#' @param phi_jumps Vector of HPO term IDs to be used as jumping distribution for proposal replacements of terms in phi during inference given gamma = 1
#' @param phi_num_leaves_geometric_rate Geometric parameter for truncated geometric distribution on number of leaf terms in phi
#' @return List (by parameter) of vectors of consecutive parameter samples from MCMC inference.
#' @examples
#' \dontrun{
#' set.seed(0)
#' data(hpo.terms)
#' disease.terms <- c("HP:0005537", "HP:0000729", "HP:0001873")
#' all.terms <- Filter(x=get.ancestors(hpo.terms, 
#'	c(disease.terms, sample(hpo.terms$id, size=50))), 
#'	f=function(tm) "HP:0000001" %in% hpo.terms$ancestors[[tm]])
#' y <- c(rep(FALSE, 96), rep(TRUE, 3))
#' x <- lapply(y, function(.y) clean.terms(
#'	hpo.terms, if (!.y) sample(all.terms, size=3) else 
#'		c(sample(all.terms, size=1), disease.terms[runif(n=3) < 0.8])))
#' sim.reg.out <- sim.reg(hpo.terms, x=x, y=y)
#' mean(sim.reg.out$gamma)
#' hpo.plot.marginal.freqs(hpo.terms, 
#'	get.term.descendancy.matrix(hpo.terms, all.terms), 
#'	sim.reg.out$phi[sim.reg.out$gamma,])
#' }
#' @export
sim.reg <- function(
	hpo.terms,
	y,
	x=NULL, 
	g=rep(0, length(y)),
	its=10000,
	thin=1,
	record_sims=FALSE,

	information.content=get.term.info.content(hpo.terms, hpo.phenotypes=x),
	ttsm=get.term.term.matrix(hpo.terms, information.content, TRUE),
	row_is_column_anc=get.term.descendancy.matrix(hpo.terms, names(information.content)),
	case_ids=unlist(mapply(SIMPLIFY=FALSE, FUN=rep, 0:(length(x)-1), sapply(x, length))),
	term_ids=as.integer(match(unlist(x), colnames(row_is_column_anc)))-1,

	return_pseudo_runs=FALSE,
	tuning_its=10000,
	tuning_burn=1000,
	burn=2000,

	gamma=(runif(1) < gamma_prior_prob),
	alpha_star=rnorm(n=1, mean=alpha_star_mean, sd=alpha_star_sd),
	alpha=rnorm(n=1, mean=alpha_mean, sd=alpha_sd),
	log_beta=rnorm(n=1, mean=log_beta_mean, sd=log_beta_sd),
	phi=sample.int(n=ncol(row_is_column_anc), size=3, replace=TRUE)-1,
	logit_mean_f=rnorm(n=1, mean=logit_mean_f_mean, sd=logit_mean_f_sd),
	log_alpha_plus_beta_f=rnorm(n=1, mean=log_alpha_plus_beta_f_mean, sd=log_alpha_plus_beta_f_sd),
	logit_mean_g=rnorm(n=1, mean=logit_mean_g_mean, sd=logit_mean_g_sd),
	log_alpha_plus_beta_g=rnorm(n=1, mean=log_alpha_plus_beta_g_mean, sd=log_alpha_plus_beta_g_sd),

	gamma_prior_prob=0.05,
	alpha_star_mean=0,
	alpha_mean=0,
	alpha_star_sd=5,
	alpha_sd=5,
	log_beta_mean=2,
	log_beta_sd=1,
	logit_mean_f_mean=1,
	logit_mean_f_sd=1,
	log_alpha_plus_beta_f_mean=2,
	log_alpha_plus_beta_f_sd=1,
	logit_mean_g_mean=0,
	logit_mean_g_sd=1.5,
	log_alpha_plus_beta_g_mean=2,
	log_alpha_plus_beta_g_sd=1,
	alpha_star_proposal_sd=2,
	alpha_proposal_sd=2,
	log_beta_proposal_sd=2,
	logit_mean_f_proposal_sd=2,
	log_alpha_plus_beta_f_proposal_sd=2,
	logit_mean_g_proposal_sd=2,
	log_alpha_plus_beta_g_proposal_sd=2,

	phi_jumps=c(0:(ncol(row_is_column_anc)-1), rep(match(unlist(lapply(x[y], get.ancestors, hpo.terms=hpo.terms)), colnames(row_is_column_anc))-1, 50)),
	pseudo_phi_marginal_prior=c(0:(ncol(row_is_column_anc)-1), rep(match(unlist(lapply(x[y], get.ancestors, hpo.terms=hpo.terms)), colnames(row_is_column_anc))-1, 50)),

	phi_num_leaves_geometric_rate=1
) {
	if (is.character(term_ids)) 
		term_ids <- match(term_ids, colnames(ttsm))-1

	if (is.character(phi_jumps)) 
		phi_jumps <- match(phi_jumps, colnames(ttsm))-1

	if (is.character(pseudo_phi_marginal_prior)) 
		pseudo_phi_marginal_prior <- match(pseudo_phi_marginal_prior, colnames(ttsm))-1

	null.out <- .sim.reg(
		hpo.terms=hpo.terms,
		ttsm=ttsm,
		row_is_column_anc=row_is_column_anc,
		y=y,
		g=g,
		case_ids=case_ids,
		term_ids=term_ids,
		its=tuning_its,
		thin=thin,
		record_sims=record_sims,
		burn=tuning_burn,

		gamma=FALSE,
		alpha_star=alpha_star,
		alpha=alpha,
		log_beta=log_beta,
		phi=phi,
		logit_mean_f=logit_mean_f,
		log_alpha_plus_beta_f=log_alpha_plus_beta_f,
		logit_mean_g=logit_mean_g,
		log_alpha_plus_beta_g=log_alpha_plus_beta_g,

		gamma_prior_prob=rep(0, length(gamma_prior_prob)),
		alpha_star_mean=alpha_star_mean,
		alpha_mean=alpha_mean,
		alpha_star_sd=alpha_star_sd,
		alpha_sd=alpha_sd,
		log_beta_mean=log_beta_mean,
		log_beta_sd=log_beta_sd,
		logit_mean_f_mean=logit_mean_f_mean,
		logit_mean_f_sd=logit_mean_f_sd,
		log_alpha_plus_beta_f_mean=log_alpha_plus_beta_f_mean,
		log_alpha_plus_beta_f_sd=log_alpha_plus_beta_f_sd,
		logit_mean_g_mean=logit_mean_g_mean,
		logit_mean_g_sd=logit_mean_g_sd,
		log_alpha_plus_beta_g_mean=log_alpha_plus_beta_g_mean,
		log_alpha_plus_beta_g_sd=log_alpha_plus_beta_g_sd,
		alpha_star_proposal_sd=alpha_star_proposal_sd,
		alpha_proposal_sd=alpha_proposal_sd,
		log_beta_proposal_sd=log_beta_proposal_sd,
		logit_mean_f_proposal_sd=logit_mean_f_proposal_sd,
		log_alpha_plus_beta_f_proposal_sd=log_alpha_plus_beta_f_proposal_sd,
		logit_mean_g_proposal_sd=logit_mean_g_proposal_sd,
		log_alpha_plus_beta_g_proposal_sd=log_alpha_plus_beta_g_proposal_sd,
		phi_jumps=phi_jumps,

		phi_num_leaves_geometric_rate=phi_num_leaves_geometric_rate
	)

	pheno.out <- .sim.reg(
		hpo.terms=hpo.terms,
		ttsm=ttsm,
		row_is_column_anc=row_is_column_anc,
		y=y,
		g=g,
		case_ids=case_ids,
		term_ids=term_ids,
		its=tuning_its,
		thin=thin,
		record_sims=record_sims,
		burn=tuning_burn,

		gamma=TRUE,
		alpha_star=alpha_star,
		alpha=alpha,
		log_beta=log_beta,
		phi=phi,
		logit_mean_f=logit_mean_f,
		log_alpha_plus_beta_f=log_alpha_plus_beta_f,
		logit_mean_g=logit_mean_g,
		log_alpha_plus_beta_g=log_alpha_plus_beta_g,

		gamma_prior_prob=rep(1, length(gamma_prior_prob)),
		alpha_star_mean=alpha_star_mean,
		alpha_mean=alpha_mean,
		alpha_star_sd=alpha_star_sd,
		alpha_sd=alpha_sd,
		log_beta_mean=log_beta_mean,
		log_beta_sd=log_beta_sd,
		logit_mean_f_mean=logit_mean_f_mean,
		logit_mean_f_sd=logit_mean_f_sd,
		log_alpha_plus_beta_f_mean=log_alpha_plus_beta_f_mean,
		log_alpha_plus_beta_f_sd=log_alpha_plus_beta_f_sd,
		logit_mean_g_mean=logit_mean_g_mean,
		logit_mean_g_sd=logit_mean_g_sd,
		log_alpha_plus_beta_g_mean=log_alpha_plus_beta_g_mean,
		log_alpha_plus_beta_g_sd=log_alpha_plus_beta_g_sd,
		alpha_star_proposal_sd=alpha_star_proposal_sd,
		alpha_proposal_sd=alpha_proposal_sd,
		log_beta_proposal_sd=log_beta_proposal_sd,
		logit_mean_f_proposal_sd=logit_mean_f_proposal_sd,
		log_alpha_plus_beta_f_proposal_sd=log_alpha_plus_beta_f_proposal_sd,
		logit_mean_g_proposal_sd=logit_mean_g_proposal_sd,
		log_alpha_plus_beta_g_proposal_sd=log_alpha_plus_beta_g_proposal_sd,
		phi_jumps=phi_jumps,

		phi_num_leaves_geometric_rate=phi_num_leaves_geometric_rate
	)
	
	result <- .sim.reg(
		hpo.terms=hpo.terms,
		ttsm=ttsm,
		row_is_column_anc=row_is_column_anc,
		y=y,
		g=g,
		case_ids=case_ids,
		term_ids=term_ids,
		its=its,
		thin=thin,
		record_sims=record_sims,
		burn=burn,

		gamma=gamma,
		alpha_star=alpha_star,
		alpha=alpha,
		log_beta=log_beta,
		phi=phi,
		logit_mean_f=logit_mean_f,
		log_alpha_plus_beta_f=log_alpha_plus_beta_f,
		logit_mean_g=logit_mean_g,
		log_alpha_plus_beta_g=log_alpha_plus_beta_g,

		gamma_prior_prob=gamma_prior_prob,

		alpha_star_mean=alpha_star_mean,
		alpha_mean=alpha_mean,
		alpha_star_sd=alpha_star_sd,
		alpha_sd=alpha_sd,
		log_beta_mean=log_beta_mean,
		log_beta_sd=log_beta_sd,
		logit_mean_f_mean=logit_mean_f_mean,
		logit_mean_f_sd=logit_mean_f_sd,
		log_alpha_plus_beta_f_mean=log_alpha_plus_beta_f_mean,
		log_alpha_plus_beta_f_sd=log_alpha_plus_beta_f_sd,
		logit_mean_g_mean=logit_mean_g_mean,
		logit_mean_g_sd=logit_mean_g_sd,
		log_alpha_plus_beta_g_mean=log_alpha_plus_beta_g_mean,
		log_alpha_plus_beta_g_sd=log_alpha_plus_beta_g_sd,
		pseudo_alpha_star_mean=mean(null.out$alpha_star[1:length(null.out$alpha_star)]),
		pseudo_alpha_star_sd=sd(null.out$alpha_star[1:length(null.out$alpha_star)]),
		pseudo_alpha_mean=mean(pheno.out$alpha[1:length(pheno.out$alpha)]),
		pseudo_alpha_sd=sd(pheno.out$alpha[1:length(pheno.out$alpha)]),
		pseudo_log_beta_mean=mean(pheno.out$log_beta[1:length(pheno.out$log_beta)]),
		pseudo_log_beta_sd=sd(pheno.out$log_beta[1:length(pheno.out$log_beta)]),
		pseudo_logit_mean_f_mean=mean(pheno.out$logit_mean_f[1:length(pheno.out$logit_mean_f)]),
		pseudo_logit_mean_f_sd=sd(pheno.out$logit_mean_f[1:length(pheno.out$logit_mean_f)]),
		pseudo_log_alpha_plus_beta_f_mean=mean(pheno.out$log_alpha_plus_beta_f[1:length(pheno.out$log_alpha_plus_beta_f)]),
		pseudo_log_alpha_plus_beta_f_sd=sd(pheno.out$log_alpha_plus_beta_f[1:length(pheno.out$log_alpha_plus_beta_f)]),
		pseudo_logit_mean_g_mean=mean(pheno.out$logit_mean_g[1:length(pheno.out$logit_mean_g)]),
		pseudo_logit_mean_g_sd=sd(pheno.out$logit_mean_g[1:length(pheno.out$logit_mean_g)]),
		pseudo_log_alpha_plus_beta_g_mean=mean(pheno.out$log_alpha_plus_beta_g[1:length(pheno.out$log_alpha_plus_beta_g)]),
		pseudo_log_alpha_plus_beta_g_sd=sd(pheno.out$log_alpha_plus_beta_g[1:length(pheno.out$log_alpha_plus_beta_g)]),
		pseudo_phi_marginal_prior=pseudo_phi_marginal_prior,
		alpha_star_proposal_sd=alpha_star_proposal_sd,
		alpha_proposal_sd=alpha_proposal_sd,
		log_beta_proposal_sd=log_beta_proposal_sd,
		logit_mean_f_proposal_sd=logit_mean_f_proposal_sd,
		log_alpha_plus_beta_f_proposal_sd=log_alpha_plus_beta_f_proposal_sd,
		logit_mean_g_proposal_sd=logit_mean_g_proposal_sd,
		log_alpha_plus_beta_g_proposal_sd=log_alpha_plus_beta_g_proposal_sd,
		phi_jumps=phi_jumps,

		phi_num_leaves_geometric_rate=phi_num_leaves_geometric_rate
	)

	if (return_pseudo_runs)
		result <- c(result, list(null.out=null.out, pheno.out=pheno.out))

	result
}


