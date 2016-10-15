promotions <- function(ontology, phi, restrict=ontology$id[!ontology$obsolete]) {
	terms <- intersect(setdiff(unique(unlist(use.names=FALSE, ontology$children[get_ancestors(ontology, phi)])), get_ancestors(ontology, phi)), restrict)
	mapply(USE.NAMES=FALSE, SIMPLIFY=FALSE, FUN=function(promotion, dontwant) intersect(restrict, setdiff(c(phi, promotion), dontwant)), terms, ontology$parents[terms])
}

demotions <- function(ontology, phi) {
	do.call(what=c, lapply(
		unique(unlist(use.names=FALSE, ontology$parents[phi])),
		function(parent) {
			phikids <- intersect(ontology$children[[parent]], phi)
			stripped <- setdiff(phi, phikids)
			if (length(phikids) == 1)
				list(union(parent, stripped))
			else
				lapply(phikids, function(leaveout) union(setdiff(phikids, leaveout), stripped))
		}
	))
}

stack_dendro <- function(d) {
	if (!is.list(d)) d
	else c(list(unlist(use.names=FALSE, d)), do.call(what=c, lapply(d, stack_dendro)))
}

min_req_signif <- function(sum_y, N, start=1, min_ratio=exp(2), p1=2, p2=1, q1=sum_y/N, q2=1) {
	pa <- start
	while (exp(lbeta(pa+p1,p2)-lbeta(p1,p2)+lbeta(sum_y-pa+q1,N-sum_y+q2)-lbeta(q1,q2)-lbeta(sum_y+q1,N-sum_y+q2)+lbeta(q1,q2)) < min_ratio & pa < sum_y) {
		pa <- pa + 1
	}
	pa	
}

dendroterms <- function(d, term_sets) {
	if (!is.list(d)) {
		structure(d, terms=term_sets[[d]])
	} else {
		children <- lapply(d, dendroterms, term_sets=term_sets)
		structure(children, terms=Reduce(f=intersect, x=lapply(children, attr, "terms")))
	}
}

dendroterms_list <- function(dt) {
	from_descs <- if (is.list(dt)) {
		do.call(c, lapply(dt, dendroterms_list))
	}
	c(list(attr(dt, "terms")), from_descs)
}

term_set_names <- function(ts) sapply(lapply(ts, sort), paste0, collapse="-")

unique_term_sets <- function(ts) {
	ts[!duplicated(term_set_names(ts))]
}

#' @importFrom ontologyIndex get_ancestors minimal_set
#' @importFrom ontologySimilarity get_sim_grid
#' @importFrom stats hclust as.dist as.dendrogram
get_phi_roots <- function(
	ontology,
	information_content,
	term_sets,
	min_intersect=2L
) {
	w_ancs <- lapply(term_sets, get_ancestors, ontology=ontology)
	grid <- get_sim_grid(ontology=ontology, information_content=information_content, term_sets=term_sets, term_sim_method="resnik")
	den <- as.dendrogram(hclust(as.dist(1-grid), method="complete"))
	clusts <- stack_dendro(den)
	dl <- dendroterms_list(dendroterms(den, w_ancs))
	unique_term_sets(lapply(dl[sapply(clusts, length) >= min_intersect], minimal_set, ontology=ontology))
}
