#' Create ontological plot of marginal probabilities of terms
#'
#' @template ontology
#' @param term_marginals Numeric vector of marginal probabilities of inclusion in \code{phi} for individual terms, named by the term IDs.
#' @param max_terms Maximum number of terms to include in plot. Note that additional terms may be included when terms have the same marginal probability, and common ancestor terms are included.
#' @param min_probability Threshold probability of inclusion in \code{phi} for triggering inclusion in plot.
#' @param ... Additional arguments to pass to \code{onto_plot}
#' @export
#' @importFrom ontologyIndex get_ancestors
#' @importFrom ontologyPlot onto_plot remove_links official_labels
plot_term_marginals <- function(
	ontology, 
	term_marginals,
	max_terms=10,
	min_probability=0.01,
	...
) {
	subjects <- names(term_marginals)[rank(-term_marginals, ties.method="min") <= max_terms & term_marginals >= min_probability]
	
	showing <- union(subjects, remove_links(ontology, get_ancestors(ontology, subjects)))

	plot_term_freqs <- ifelse(
		showing %in% names(term_marginals),
		term_marginals[showing],
		0
	)

	onto_plot(
		ontology,
		terms=showing, 
		label=paste(official_labels(ontology, showing), paste(round(plot_term_freqs, 2), sep=""), sep="\n"),
		width=sqrt(plot_term_freqs)*(3/max(sqrt(plot_term_freqs))),
		...
	)
}
