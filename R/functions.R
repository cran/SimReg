#' Get term-term similarity matrix based on Lin's definition
#'
#' @template hpo.terms
#' @template information.content 
#' @template term.descendancy.matrix
#' @param prune Logical indicating whether to make similarity of more specific terms 0 in one direction
#' @return Numeric matrix of term-term similarities
#' @references Lin, D. (1989). An Information-Theoretic Definition of Similarity.
#' @export
get.term.term.matrix <- function(hpo.terms, information.content, term.descendancy.matrix=get.term.descendancy.matrix(hpo.terms, names(information.content)), prune=FALSE) {
	result <- matrix(0, nrow=length(information.content), ncol=length(information.content), dimnames=rep(times=2, list(names(information.content))))

	for (term in names(information.content)) {
		descs <- c(match(term, names(information.content)), which(term.descendancy.matrix[term,]))
		result[descs,descs] <- pmax(matrix(information.content[term], nrow=length(descs), ncol=length(descs)), result[descs,descs])
	}

	ic.sum <- outer(information.content, information.content, "+")

	result <- 2 * result / ic.sum

	result[which(arr.ind=TRUE, is.na(result))] <- 0

	result %>%
	(function(ttsm.mat) {
		if (!prune) return(ttsm.mat)
		result <- ttsm.mat
		for (term in rownames(ttsm.mat)) 
			result[term, ] <- ifelse(colnames(ttsm.mat) %in% hpo.terms$ancestors[[term]], result[term,], 0)
		result
	})
}

#' Get leaf matrix
#'
#' Procure logical matrix from character matrix of HPO terms, indicating whether each element is a leaf in the context of the row. Typically used on the output of the \code{\link{sim.reg}}.
#'
#' @param term.descendancy.matrix Logical term descendancy matrix, dimensions symmetrically labelled by terms, and where by a cell value of TRUE indicates that the row is the ancestor of the column term (in the sense of the DAG structure of the HPO)
#' @param terms.matrix The character matrix of HPO terms
#' @return Logical matrix the same dimensions as the terms.matrix which indicates whethere each element is a leaf or not
#' @examples
#' data(hpo.terms)
#' as.row.leaves(
#' 	get.term.descendancy.matrix(hpo.terms, c("HP:0001873", "HP:0001872")),
#' 	matrix(c("HP:0001873","HP:0001872","HP:0001873","HP:0001872"),2,2)
#' )
as.row.leaves <- function(term.descendancy.matrix, terms.matrix) {
	terms <- unique(as.vector(terms.matrix))
	term.num.matrix <- apply(
		terms.matrix,
		2,
		function(x) as.integer(factor(x, levels=terms))
	) - 1

	.Call(
		"leaf_matrix",
		term.descendancy.matrix[terms,terms,drop=FALSE],
		term.num.matrix,
		PACKAGE="SimReg"
	)
}

