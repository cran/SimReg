#' @importFrom ontologyIndex get_term_info_content get_term_descendancy_matrix get_ancestors
#' @importFrom ontologySimilarity get_term_sim_mat
#' @importFrom stats sd runif rnorm
#' @importFrom Rcpp evalCpp
#' @useDynLib SimReg
sim_reg_gamma1 <- function(
	ontology=NULL,
	y,
	x=NULL, 
	g=rep(0, length(y)),
	its=10000,
	thin=1,
	burn=2000,
	record_sims=FALSE,

	information_content=get_term_info_content(ontology, term_sets=x),
	term_sim_mat=get_term_sim_mat(ontology, information_content, TRUE),

	term_descendancy_matrix=get_term_descendancy_matrix(ontology, names(information_content)),
	case_ids=unlist(mapply(SIMPLIFY=FALSE, FUN=rep, 0:(length(x)-1), sapply(x, length))),
	term_ids=as.integer(match(unlist(x), colnames(term_descendancy_matrix)))-1,

	gamma=(runif(1) < gamma_prior_prob)[1],
	alpha_star=rnorm(n=1, mean=alpha_star_mean, sd=alpha_star_sd),
	alpha=rnorm(n=1, mean=alpha_mean, sd=alpha_sd),
	log_beta=rnorm(n=1, mean=log_beta_mean, sd=log_beta_sd),
	phi=sample.int(n=ncol(term_descendancy_matrix), size=3, replace=TRUE)-1,
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
	phi_jumps=c(0:(ncol(term_descendancy_matrix)-1), rep(match(unlist(lapply(x[y], get_ancestors, ontology=ontology)), colnames(term_descendancy_matrix))-1, 50)),
	pseudo_phi_marginal_prior=c(0:(ncol(term_descendancy_matrix)-1), rep(match(unlist(lapply(x[y], get_ancestors, ontology=ontology)), colnames(term_descendancy_matrix))-1, 50)),
	phi_num_leaves_geometric_rate=1,
	lit_sims=setNames(rep(1, ncol(term_descendancy_matrix)), colnames(term_sim_mat))
) {
	result <- .Call(
		"R_sim_reg",
		PACKAGE="SimReg",
		its,
		thin,
		record_sims,
		term_sim_mat,
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

		lit_sims[colnames(term_sim_mat)],
		term_descendancy_matrix,
		phi_num_leaves_geometric_rate,
		FALSE,
		FALSE,
		0,
		c(0, 0),
		0,
		0,
		1
	)

	trunc <- function(y) if (class(y) == "matrix") { if (nrow(y) == its) y[(burn+1):(nrow(y)),,drop=FALSE] else y } else { y[(burn+1):(length(y))] }

	result <- c(
		list(liks=lapply(result$liks, trunc)),
		lapply(result[setdiff(names(result), c("H", "liks"))], trunc)
	)

	result$phi_vector <- apply(result$phi, 2, function(x) colnames(term_descendancy_matrix)[x+1])
	result$phi <- mapply(SIMPLIFY=FALSE, FUN="[", split(result$phi_vector, seq(nrow(result$phi_vector))), split(as_row_leaves(term_descendancy_matrix, result$phi_vector), seq(nrow(result$phi_vector))))

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
			list(phi_accept=sapply(1:(nrow(result$phi_vector)-1), function(x) !identical(result$phi_vector[x,], result$phi_vector[x+1,]))),
			lapply(
				local({
					x <- c("alpha_star", "alpha", "log_beta", "logit_mean_f", "log_alpha_plus_beta_f", "logit_mean_g", "log_alpha_plus_beta_g")
					setNames(x, paste(x, "_accept", sep=""))
				}),
				function(x) result[[x]][1:(length(result[[x]])-1)] != result[[x]][2:length(result[[x]])]
			)
		)
	)
}

#' Similarity regression
#'
#' Performins Bayesian `similarity regression' on given binary genotype \code{y} (logical vector) against ontological term sets \code{x} (list of character vectors of term IDs). This could, for example, be a \code{list} of character vectors of HPO term IDs representing case phenotypes. It returns a list of traces for the sampled parameters. Of particular interest are the estimated mean posteriors of \code{gamma} (the model selection indicator, thus giving an estimate of the probability of an association under the model assumptions - stored in the `mean_posterior_gamma' slot in the result, i.e. \code{result$mean_posterior_gamma} (which can also be calculated \code{mean(result$gamma)}), and the characteristic ontological profile phi (which can be visualised by the functions \code{\link{phi_plot}}, \code{\link{term_pair_marginals_plot}}, and \code{\link{term_marginals}}).
#'
#' @template ontology
#' @param y Logical vector of genotypes (typically 1 for rare genotype, 0 for common genotype).
#' @param x List of character vectors of terms IDs,
#' @param g Genotype log odds offset per individual.
#' @param its Number of update cycles to perform .
#' @param thin Factor by which to thin resultant chains of parameter samples.
#' @param record_sims Logical indicating whether to record trace of similarities.
#' @param information_content Numeric vector, named by HPO IDs, containing the information content of corresponding terms.
#' @param term_sim_mat The `term-term' similarity matrix, a numeric matrix whose dimensions are named by the terms so that cell i,j contains the similarity of term i to term j.
#' @param term_descendancy_matrix Logical matrix, whose dimensions are named by HPO term IDs so that cell i,j is TRUE if i is an ancestor term of j.
#' @param case_ids IDs for the N cases from 0 to N-1, indicating which case terms in \code{term_ids} belong to (automatically determined given x).
#' @param term_ids Vector of HPO term IDs belonging to cases.
#' @param return_tuning_runs Logical indicating whether to return the MCMC output of the tuning phase of the inference procedure.
#' @param tuning_its Number of update cycles to perform in the tuning phase of the inference procedure.
#' @param tuning_burn Number of update cycles to discard in tuning phase.
#' @param burn Number of update cycles to discard .
#' @param gamma Initial value of model selection indicator gamma..
#' @param alpha_star  Initial value of alpha_star, the rate of observing the rare genotype y = 1 under gamma = 0, i.e. the no association model .
#' @param alpha Initial value of alpha, the background rate of observing the rare genotype under gamma = 1.
#' @param log_beta Initial value of log_beta, the log of the effect size of onotological similarity.
#' @param phi Character vector of HPO term IDs giving the initial value of phi, the characteristic phenotype.
#' @param logit_mean_f Initial value of logit_mean_f.
#' @param log_alpha_plus_beta_f Initial value of log_alpha_plus_beta_f.
#' @param logit_mean_g Initial value of logit_mean_g.
#' @param log_alpha_plus_beta_g Initial value of log_alpha_plus_beta_g.
#' @param gamma_prior_prob Prior probability of gamma = 1.
#' @param alpha_star_mean Prior mean of alpha_star given gamma = 0.
#' @param alpha_mean Prior mean of alpha given gamma = 1.
#' @param alpha_star_sd Prior sd of alpha_star given gamma = 0.
#' @param alpha_sd Prior sd of alpha given gamma = 1.
#' @param log_beta_mean Prior mean of log_beta given gamma = 1.
#' @param log_beta_sd Prior sd of log_beta given gamma = 1.
#' @param logit_mean_f_mean Prior mean of logit_mean_f given gamma = 1.
#' @param logit_mean_f_sd Prior sd of logit_mean_f given gamma = 1.
#' @param log_alpha_plus_beta_f_mean Prior mean of log_alpha_plus_beta_f given gamma = 1.
#' @param log_alpha_plus_beta_f_sd Prior sd of log_alpha_plus_beta_f given gamma = 1.
#' @param logit_mean_g_mean Prior mean of logit_mean_g given gamma = 1.
#' @param logit_mean_g_sd Prior sd of logit_mean_g given gamma = 1.
#' @param log_alpha_plus_beta_g_mean Prior mean of log_alpha_plus_beta_g given gamma = 1.
#' @param log_alpha_plus_beta_g_sd Prior sd of log_alpha_plus_beta_g given gamma = 1.
#' @param pseudo_phi_marginal_prior Vector of HPO term IDs to be used as prior distribution on marginal probability of single term in phi given gamma = 0.
#' @param alpha_star_proposal_sd Proposal sd of local jumps in MH updates of alpha_star used during inference.
#' @param alpha_proposal_sd Proposal sd of local jumps in MH updates of alpha used during inference.
#' @param log_beta_proposal_sd Proposal sd of local jumps in MH updates of log_beta used during inference.
#' @param logit_mean_f_proposal_sd Proposal sd of local jumps in MH updates of logit_mean_f used during inference.
#' @param log_alpha_plus_beta_f_proposal_sd Proposal sd of local jumps in MH updates of log_alpha_plus_beta_f used during inference.
#' @param logit_mean_g_proposal_sd Proposal sd of local jumps in MH updates of logit_mean_g used during inference.
#' @param log_alpha_plus_beta_g_proposal_sd Proposal sd of local jumps in MH updates of log_alpha_plus_beta_g used during inference.
#' @param phi_jumps Vector of HPO term IDs to be used as jumping distribution for proposal replacements of terms in phi during inference given gamma = 1.
#' @param phi_num_leaves_geometric_rate Geometric parameter for truncated geometric distribution on number of leaf terms in phi.
#' @param lit_sims Numeric vector of similarities (greater than 0) of literature phenotype to individual terms (named by term ID).
#' @return List (by parameter) of vectors of consecutive parameter samples from MCMC inference.
#' @examples
#' \dontrun{
#' set.seed(0)
#' data(hpo)
#' disease_terms <- c("HP:0005537", "HP:0000729", "HP:0001873")
#' all_terms <- get_ancestors(hpo, 
#'	c(disease_terms, sample(hpo$id, size=50)))
#' y <- c(rep(FALSE, 96), rep(TRUE, 3))
#' x <- lapply(y, function(.y) minimal_set(
#'	hpo, if (!.y) sample(all_terms, size=3) else 
#'		c(sample(all_terms, size=1), disease_terms[runif(n=3) < 0.8])))
#' sim_reg_out <- sim_reg(ontology=hpo, x=x, y=y)
#' mean(sim_reg_out$gamma)
#' phi_plot(hpo, 
#'	sim_reg_out$phi[sim_reg_out$gamma])
#' }
#' @export
#' @importFrom stats rnorm runif
#' @importFrom ontologyIndex get_term_descendancy_matrix
sim_reg <- function(
	ontology=NULL,
	y,
	x=NULL, 
	g=rep(0, length(y)),
	its=10000,
	thin=1,
	record_sims=FALSE,

	information_content=get_term_info_content(ontology, term_sets=x),
	term_sim_mat=prune_sim_mat(ontology=ontology, get_term_sim_mat(ontology=ontology, information_content=information_content)),
	term_descendancy_matrix=get_term_descendancy_matrix(ontology, colnames(term_sim_mat)),
	case_ids=unlist(mapply(SIMPLIFY=FALSE, FUN=rep, 0:(length(x)-1), sapply(x, length))),
	term_ids=as.integer(match(unlist(x), colnames(term_descendancy_matrix)))-1,

	return_tuning_runs=FALSE,
	tuning_its=10000,
	tuning_burn=1000,
	burn=2000,

	gamma=(runif(1) < gamma_prior_prob),
	alpha_star=rnorm(n=1, mean=alpha_star_mean, sd=alpha_star_sd),
	alpha=rnorm(n=1, mean=alpha_mean, sd=alpha_sd),
	log_beta=rnorm(n=1, mean=log_beta_mean, sd=log_beta_sd),
	phi=sample.int(n=ncol(term_descendancy_matrix), size=3, replace=TRUE)-1,
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

	phi_jumps=c(0:(ncol(term_descendancy_matrix)-1), rep(match(unlist(lapply(x[y], get_ancestors, ontology=ontology)), colnames(term_descendancy_matrix))-1, 50)),
	pseudo_phi_marginal_prior=c(0:(ncol(term_descendancy_matrix)-1), rep(match(unlist(lapply(x[y], get_ancestors, ontology=ontology)), colnames(term_descendancy_matrix))-1, 50)),

	phi_num_leaves_geometric_rate=1,
	lit_sims=setNames(rep(1, ncol(term_descendancy_matrix)), colnames(term_sim_mat))
) {
	if (is.character(term_ids)) 
		term_ids <- match(term_ids, colnames(term_sim_mat))-1

	if (is.character(phi_jumps)) 
		phi_jumps <- match(phi_jumps, colnames(term_sim_mat))-1

	if (is.character(pseudo_phi_marginal_prior)) 
		pseudo_phi_marginal_prior <- match(pseudo_phi_marginal_prior, colnames(term_sim_mat))-1

	null.out <- sim_reg_gamma1(
		ontology=ontology,
		term_sim_mat=term_sim_mat,
		term_descendancy_matrix=term_descendancy_matrix,
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

		phi_num_leaves_geometric_rate=phi_num_leaves_geometric_rate,
		lit_sims=lit_sims
	)

	pheno.out <- sim_reg_gamma1(
		ontology=ontology,
		term_sim_mat=term_sim_mat,
		term_descendancy_matrix=term_descendancy_matrix,
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

		phi_num_leaves_geometric_rate=phi_num_leaves_geometric_rate,
		lit_sims=lit_sims
	)
	
	H <- (mean(null.out$liks$likelihood) - mean(pheno.out$liks$likelihood))

	result <- sim_reg_gamma1(
		ontology=ontology,
		term_sim_mat=term_sim_mat,
		term_descendancy_matrix=term_descendancy_matrix,
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

		phi_num_leaves_geometric_rate=phi_num_leaves_geometric_rate,
		lit_sims=lit_sims
	)

	if (return_tuning_runs) {
		result$null_out <- null.out
		result$pheno_out <- pheno.out
	}

	result$mean_posterior_gamma <- mean(result$gamma)

	result
}


