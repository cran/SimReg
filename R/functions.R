#' Prune between-term similarity matrix
#'
#' Prune similarity of more specific terms to 0
#'
#' @template ontology
#' @template term_sim_mat
#' @return Numeric matrix of term-term similarities
#' @export
prune_sim_mat <- function(ontology, term_sim_mat) {
	result <- term_sim_mat
	for (term in rownames(result)) 
		result[term, ] <- ifelse(colnames(result) %in% ontology$ancestors[[term]], result[term,], 0)
	result
}

#' Get leaf matrix
#'
#' Procure logical matrix from character matrix of ontological terms, indicating whether each element is a leaf in the context of the row. Typically used to map the sampled vectors of the ontological parameter phi from the \code{\link{sim_reg}} procedure.
#'
#' @param term_descendancy_matrix Logical term_descendancy_matrix, dimensions symmetrically labelled by terms, and where by a cell value of TRUE indicates that the row is the ancestor of the column term (in the sense of the DAG structure of the HPO)
#' @param terms_matrix The character matrix of HPO terms
#' @return Logical matrix the same dimensions as the terms_matrix which indicates whethere each element is a leaf or not
#' @examples
#' library(ontologyIndex)
#' data(hpo)
#' as_row_leaves(
#' 	get_term_descendancy_matrix(hpo, c("HP:0001873", "HP:0001872")),
#' 	matrix(c("HP:0001873","HP:0001872","HP:0001873","HP:0001872"),2,2)
#' )
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib SimReg
as_row_leaves <- function(term_descendancy_matrix, terms_matrix) {
	terms <- unique(as.vector(terms_matrix))
	term.num.matrix <- apply(
		terms_matrix,
		2,
		function(x) as.integer(factor(x, levels=terms))
	) - 1

	.Call(
		"leaf_matrix",
		term_descendancy_matrix[terms,terms,drop=FALSE],
		term.num.matrix,
		PACKAGE="SimReg"
	)
}

#' @importFrom Rcpp evalCpp
#' @useDynLib SimReg
.log_odds_trace <- function(
	tsm,
	row_is_column_anc,
	as_probs=FALSE,
	x=NULL,
	case_ids=unlist(mapply(SIMPLIFY=FALSE, FUN=rep, 0:(length(x)-1), sapply(x, length))),
	term_ids=as.integer(match(unlist(x), colnames(row_is_column_anc)))-1,
	g,
	phi,
	gamma,
	alpha_star,
	alpha,
	log_beta,
	logit_mean_f,
	log_alpha_plus_beta_f,
	logit_mean_g,
	log_alpha_plus_beta_g
) {
	result <- .Call(
		"R_log_odds_trace",
		PACKAGE="SimReg",
		tsm,
		row_is_column_anc,
		term_ids,
		case_ids,
		g,
		phi,
		gamma,
		alpha_star,
		alpha,
		log_beta,
		logit_mean_f,
		log_alpha_plus_beta_f,
		logit_mean_g,
		log_alpha_plus_beta_g
	)

	if (as_probs)
		result <- 1 - 1/(1+exp(result))

	result
}	

#' Get trace of log odds of observing the rare genotype (y = 1) for individual cases from the output of \code{\link{sim_reg}}
#'
#' @param term_sim_mat Numeric matrix of similarities between individual terms, typically created with \code{get_term_sim_mat}
#' @template term_descendancy_matrix
#' @template term_sets
#' @param output Output of \code{\link{sim_reg}} function, giving trace of parameter values obtained from MCMC sample
#' @param g Genotype log odds offsets per individual
#' @param as_probs Boolean value indicating whether to convert the log odds to probabilities
#' @return Numeric matrix of log odds trace per individual
#' @export
log_odds_trace <- function(
	term_sim_mat, 
	term_descendancy_matrix, 
	term_sets, 
	output, 
	g=rep(FALSE, length(term_sets)), 
	as_probs=FALSE
) {
	.log_odds_trace(
		as_probs=as_probs,
		row_is_column_anc=term_descendancy_matrix,
		x=term_sets,
		tsm=term_sim_mat,
		g=g,
		gamma=output$gamma,
		alpha=output$alpha,
		alpha_star=output$alpha_star,
		log_beta=output$log_beta,
		logit_mean_f=output$logit_mean_f,
		log_alpha_plus_beta_f=output$log_alpha_plus_beta_f,
		logit_mean_g=output$logit_mean_g,
		log_alpha_plus_beta_g=output$log_alpha_plus_beta_g,
		phi=t(apply(output$phi_vector, 1, function(x) match(x, colnames(term_descendancy_matrix))-1))
	)
}

