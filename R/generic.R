#' Print \code{sim_reg_summary} object
#'
#' @param x Object of class \code{sim_reg_summary}.
#' @param ... Non-used arguments.
#' @export
#' @method print sim_reg_summary
print.sim_reg_summary <- function(x, ...) {
	width <- getOption("width")
	dashed <- paste0(rep("-", width), collapse="")
	cat("Probability of association: ", round(digits=2, x[["prob_association"]]), "\n")

	cat(dashed, "\n")
	print(data.frame(
		Name=substr(x$term_names, 1, max(10, as.integer(width)-20)),
		p=round(x$term_marginals, digits=2),
		row.names=names(x$term_marginals)
	)[seq(min(length(x$term_marginals), 10)),,drop=FALSE], row.names=FALSE)

	cat(dashed, "\n")
}

#' Print \code{sim_reg_output} object
#'
#' @param x Object of class \code{sim_reg_output}.
#' @param ... Non-used arguments.
#' @export
#' @method print sim_reg_output
print.sim_reg_output <- function(x, ...) {
	stopifnot(class(x) == "sim_reg_output")
	print.sim_reg_summary(summary(x))
}

#' Get summary of \code{sim_reg_output} object
#'
#' @param object Object of class \code{sim_reg_output}.
#' @param prior Prior probability of association.
#' @param ... Non-used arguments.
#' @export
#' @importFrom stats setNames
#' @method summary sim_reg_output
summary.sim_reg_output <- function(object, prior=0.05, ...) {
	tm <- get_term_marginals(object)
	structure(
		class="sim_reg_summary",
		list(
			prob_association=prob_association(object),
			term_marginals=tm,
			term_names=object$term_names_used[names(tm)]
		)
	)
}

#' Plot summary of \code{sim_reg_output} object
#'
#' @param x Object of class \code{sim_reg_summary}.
#' @param ... Additional arguments to pass to \code{\link{plot_term_marginals}}.
#' @export
#' @method plot sim_reg_summary
#' @importFrom graphics plot
plot.sim_reg_summary <- function(x, ...) {
	plot(main=paste0("Probability of association: ", round(digits=2, x$prob_association)), plot_term_marginals(term_marginals=x$term_marginals, ...))
}

#' @export
#' @rdname plot.sim_reg_summary
#' @method plot sim_reg_output
plot.sim_reg_output <- function(x, ...) {
	plot.sim_reg_summary(summary(x), ...)
}
