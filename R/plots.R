#' Replace term ID names of vector/array with full term names
#'
#' @template ontology
#' @param x Vector/array named by term IDs
#' @export
full_names <- function(ontology, x) {
	if (is.array(x)) 
		dimnames(x) <- lapply(dimnames(x), function(i) ontology$name[i])
	else
		names(x) <- ontology$name[names(x)]
	x
}

#' Get marginal frequencies of single terms in phi
#'
#' @template phi
#' @return Numeric vector of proportion of samples each term was included in phi.
#' @export
term_marginals <- function(phi) {
	x <- sort(decreasing=TRUE, table(unlist(phi)))
	x/length(phi)
}

#' Get marginal frequencies of single terms in phi
#'
#' @template term_descendancy_matrix
#' @param phi_vector_trace Delta trace output from \code{\link{sim_reg}} function, where the rows are the values of phi at each iteration
#' @return Numeric vector of proportion of samples each term was included in phi.
term_marginals_from_phi_vec <- function(term_descendancy_matrix, phi_vector_trace) 
	with(data=as.numeric(sort(decreasing=TRUE, table(phi_vector_trace[as_row_leaves(term_descendancy_matrix, phi_vector_trace)]))), expr=x/sum(x))
 
#' Get matrix of presence of term pairs
#'
#' @template phi
#' @param terms_to_use Specify the terms to appear in the plot.
#' @param symmetric Logical value determining whether output is a symmetric matrix (\code{TRUE}) or lower triangular one (\code{FALSE}).
#' @param n Maximum number of terms to include in the matrix.
#' @return Matrix of term-pair frequencies (where by cell i,j contains the frequency of inclusion in phi for term pair i,j).
#' @export
#' @importFrom stats setNames
#' @importFrom utils combn 
term_pair_marginals <- function(
	phi,
	terms_to_use=NULL,
	symmetric=FALSE,
	n=10
) {
	all.terms <- unique(unlist(phi))
	n <- min(n, length(all.terms))
	selected <- 
		if (is.null(terms_to_use))
			names(sort(decreasing=TRUE, table(unlist(phi)))[1:n])
		else
			terms_to_use

	selected.full.trace <- sapply(phi, function(x) setNames(selected %in% x, selected)) 
	selected.trace <- selected.full.trace[,apply(selected.full.trace, 2, function(y) sum(y) >= 2)]

	term.pairs <- t(combn(selected, 2))
	heatmap.mat <- matrix(0, n, n, dimnames=rep(list(selected), 2))
	heatmap.mat[term.pairs] <- apply(
		term.pairs,
		1,
		function(pair) sum(apply(
			selected.trace[pair,],
			2,
			all
		))/length(phi)
	)

	heatmap.mat <- heatmap.mat + t(heatmap.mat)

	if (!symmetric) heatmap.mat[lower.tri(heatmap.mat, diag=TRUE)] <- NA

	heatmap.mat
}

#' Create `ontology_plot' of phi
#'
#' Create plot of marginal frequency of individual terms in context of related terms in the ontology, ready for plotting to device or exporting to dot file.
#'
#' @template ontology
#' @template phi
#' @param max_terms Specify the maximum number of terms to appear in the plot.
#' @param min_frequency Threshold frequency for including terms in plot.
#' @param colour_gradient Logical indicating whether to colour terms in the plot according to their marginal frequencies (blue being the least frequent, yellow the most).
#' @param size_gradient Logical indicating whether to colour terms in the plot according to their marginal frequencies.
#' @param custom_labels Character vector of custom labels for terms (named by corresponding term IDs).
#' @param show_proportion Logical indicating whether to append the `inclusion in phi' rate to the labels in the terms.
#' @param fillcolor Vector of colours (named by HPO term IDs) for corresponding terms in plot.
#' @param ... Additional parameters to be passed to \code{onto_plot}.
#' @return Plots graph.
#' @export
#' @importFrom grDevices colorRampPalette
#' @importFrom ontologyPlot remove_links simple_labels onto_plot
#' @importFrom ontologyIndex get_ancestors
phi_plot <- function(
	ontology, 
	phi,
	max_terms=10,
	min_frequency=0,
	colour_gradient=FALSE,
	size_gradient=TRUE,
	custom_labels=c(),
	show_proportion=TRUE,
	fillcolor=NULL,
	...
) {
	marginal_freqs <- local({ x <- term_marginals(phi); y <- x[1:min(length(x), max_terms)]; y[y > min_frequency] })
	plot_term_freqs <- local({
		missing <- setdiff(remove_links(ontology, get_ancestors(ontology, names(marginal_freqs))), names(marginal_freqs))
		setNames(
			c(marginal_freqs, rep(0, length(missing))),
			c(names(marginal_freqs), missing)
		)
	})

	onto_plot(
		ontology,
		terms=names(plot_term_freqs), 
		label=
			if (show_proportion) paste(ifelse(!(names(plot_term_freqs) %in% names(custom_labels)), simple_labels(ontology, names(plot_term_freqs)), custom_labels[names(plot_term_freqs)]), paste(round(plot_term_freqs, 2), sep=""), sep="\n")
			else ifelse(!(names(plot_term_freqs) %in% names(custom_labels)), simple_labels(ontology, names(plot_term_freqs)), custom_labels[names(plot_term_freqs)])
		,
		width="if"(
			size_gradient,
			sqrt(plot_term_freqs)*(3/max(sqrt(plot_term_freqs))),
			rep(1, length(plot_term_freqs))
		),
		fillcolor=if (is.null(fillcolor)) "if"(
			colour_gradient,
			colorRampPalette(c("#0099FF", "green3", "Yellow"))(20)[local({ x <- cut(plot_term_freqs, breaks=seq(from=min(plot_term_freqs), to=max(plot_term_freqs), by=diff(range(plot_term_freqs))/20), include.lowest=TRUE, labels=FALSE); x[which(plot_term_freqs == max(plot_term_freqs))] <- 20; x })],
			rep("cyan", length(plot_term_freqs))
		) else fillcolor,
		...
	)
}
 
#' Plot SimReg output
#'
#' Plot the output of the call to the \code{sim_reg} function. Output contains plots of marginal probabilities of terms and pairs of terms being present phi and histograms giving the distributions of variables. Note, a large device is required for a successful call.
#'
#' @param x Output of call to \code{sim_reg}.
#' @template ontology
#' @param ... Non-used arguments.
#' @return Plots
#' @export
plot.sim_reg_samples <- function(x, ontology, ...) {
	plot(summary(x), ontology=ontology, ...)
}

#' Create grid of SimReg output plots
#'
#' Create a PDF containing a set of summary plots for the output of the \code{sim_reg} application. Output contains plots of marginal probabilities of terms and pairs of terms being present phi and histograms giving the distributions of variables.
#'
#' @template ontology
#' @param file_name File to write plots to 
#' @param sim_reg_out Output of call to \code{sim_reg}
#' @return Plots graph to file
#' @export
#' @importFrom grDevices dev.off pdf
sim_reg_summary <- function(ontology, file_name, sim_reg_out) {
	pdf(file_name, width=25, height=35)
	plot(sim_reg_out, ontology=ontology)
	dev.off()
}
