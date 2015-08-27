#' Get matrix of presence of term pairs
#'
#' @template term.descendancy.matrix
#' @template phi.trace
#' @param terms.to.use Specify the terms to appear in the plot
#' @param n Maximum number of terms to include in the matrix
#' @return Matrix of term-pair frequencies (where by cell i,j contains the frequency of inclusion in phi for term pair i,j)
leaf.term.pair.marginal.matrix <- function(
	term.descendancy.matrix,
	phi.trace,
	terms.to.use=NULL,
	n=10
) {
	is.leaf.matrix <- as.row.leaves(term.descendancy.matrix, phi.trace)
	all.leaves <- unlist(phi.trace[is.leaf.matrix])
	all.terms <- unique(all.leaves)
	n <- min(n, length(all.terms))
	selected <- "if"(
		is.null(terms.to.use),
		table(all.leaves) %>% 
		sort(decreasing=TRUE) %>% 
		(function(x) x[1:n]) %>% 
		names,
		terms.to.use
	)

	selected.trace <- 
		sapply(1:nrow(phi.trace), function(x) setNames(selected %in% phi.trace[x,][is.leaf.matrix[x,]], selected)) %>%
		(function(x) {
			x[,apply(x, 2, function(y) sum(y) >= 2)]
		})

	term.pairs <- combn(selected, 2) %>% t
	heatmap.mat <- matrix(0, n, n, dimnames=rep(list(selected), 2))
	heatmap.mat[term.pairs] <- apply(
		term.pairs,
		1,
		function(pair) sum(apply(
			selected.trace[pair,],
			2,
			all
		))/nrow(phi.trace)
	)
	heatmap.mat <- heatmap.mat + t(heatmap.mat)

	heatmap.mat
}

#' Plot marginal frequency of terms
#'
#' Plot graphic depicting marginal frequency of individual HPO terms in context of related terms in the ontology
#'
#' @template hpo.terms
#' @template term.descendancy.matrix
#' @template phi.trace
#' @param max.terms Specify the maximum number of terms to appear in the plot
#' @param colour.gradient Logical indicating whether to colour terms in the plot according to their marginal frequencies (blue being the least frequent, yellow the most)
#' @param size.gradient Logical indicating whether to colour terms in the plot according to their marginal frequencies
#' @param custom.labels Character vector of custom labels for terms (named by corresponding term IDs)
#' @param show.proportion Logical indicating whether to append the `inclusion in phi' rate to the labels in the terms
#' @param colours Vector of colours (named by HPO term IDs) for corresponding terms in plot
#' @param ... Additional parameters to be passed to \code{hpo.plot}
#' @return Plots graph
hpo.plot.marginal.freqs <- function(
	hpo.terms, 
	term.descendancy.matrix,
	phi.trace,
	max.terms=20,
	colour.gradient=TRUE,
	size.gradient=TRUE,
	custom.labels=c(),
	show.proportion=TRUE,
	colours=NULL,
	...
) {
	marginal.freqs <- 
		phi.trace[as.row.leaves(term.descendancy.matrix, phi.trace)] %>%
		unlist %>%
		table %>%
		sort(decreasing=TRUE) %>%
		(function(x) x[1:min(max.terms, length(x))]) %>%
		(function(x) setNames(
			as.numeric(x)/nrow(phi.trace),
			names(x)
		)) %>%
		(function(x) {
			missing <- setdiff(remove.links(hpo.terms, get.ancestors(hpo.terms, names(x))), names(x))
			setNames(
				c(x, rep(0, length(missing))),
				c(names(x), missing)
			)
		})

	hpo.plot(
		hpo.terms,
		terms=names(marginal.freqs), 
		labels=
			if (show.proportion) paste(ifelse(!(names(marginal.freqs) %in% names(custom.labels)), get.simple.node.labels(hpo.terms, names(marginal.freqs)), custom.labels[names(marginal.freqs)]), paste(round(marginal.freqs, 2), sep=""), sep="\\\n")
			else ifelse(!(names(marginal.freqs) %in% names(custom.labels)), get.simple.node.labels(hpo.terms, names(marginal.freqs)), custom.labels[names(marginal.freqs)])
		,
		sizes="if"(
			size.gradient,
			sqrt(marginal.freqs)*(3/max(sqrt(marginal.freqs))),
			rep(1, length(marginal.freqs))
		),
		colours=if (is.null(colours)) "if"(
			colour.gradient,
			colorRampPalette(c("#0099FF", "green3", "Yellow"))(20)[cut(marginal.freqs, breaks=seq(from=min(marginal.freqs), to=max(marginal.freqs), by=diff(range(marginal.freqs))/20), include.lowest=TRUE, labels=FALSE) %>% (function(x) { x[which(marginal.freqs == max(marginal.freqs))] <- 20; x })],
			rep("cyan", length(marginal.freqs))
		) else colours,
		...
	)
}

#' Plot marginal distribution of single HPO terms
#'
#' Plots the marginal frequency of single HPO terms' inclusion in the phi parameter over the course of an application of \code{\link{sim.reg}} as bar chart.
#'
#' @template hpo.terms
#' @template term.descendancy.matrix
#' @param phi.trace Delta trace output from \code{\link{sim.reg}} function, where the rows are the values of phi at each iteration
#' @param use.clinical.names Logical indicating whether to use full clinical names for the terms (TRUE) or just the codes (FALSE)
#' @param n Maximum number of terms for plot
#' @return ggplot2 plot object
single.term.marginals.plot <- function(hpo.terms, term.descendancy.matrix, phi.trace, use.clinical.names=FALSE, n=10) 
	table(phi.trace[as.row.leaves(term.descendancy.matrix, phi.trace)]) %>%
	sort(decreasing=TRUE) %>%
	(function(x) x[1:min(n, length(x))]) %>%
	(function(x) data.frame(
		Term=names(x),
		Marginal.Freq=as.numeric(x)/sum(as.numeric(x))
	)) %>% 
	(function(x) "if"(
		use.clinical.names,
		transform(x, Term=get.shortened.names(hpo.terms, as.character(x$Term))),
		x
	)) %>%
	ggplot(
		mapping=aes_string(
			x="Term",
			y="Marginal.Freq"
		)
	) + 
	geom_bar(stat="identity") + 
	theme(axis.text.x = element_text(angle=90, hjust=1))
 
#' Get marginal distribution of HPO term pairs
#'
#' Get the marginal frequency of HPO terms pairs' inclusion in the phi parameter over the course of an application of \code{\link{sim.reg}} as matrix
#'
#' @template hpo.terms
#' @template term.descendancy.matrix
#' @param phi.trace Delta trace output from \code{\link{sim.reg}} function, where the rows are the values of phi at each iteration
#' @param use.clinical.names Logical indicating whether to use full clinical names for the terms (TRUE) or just the codes (FALSE)
#' @param n Maximum number of terms for plot
#' @param terms.to.use Specify the terms to appear in the plot
#' @param limits Range of frequencies to set colours for
#' @return ggplot2 plot object
get.two.term.marginals <- function(
	hpo.terms,
	term.descendancy.matrix, 
	phi.trace, 
	use.clinical.names=FALSE,
	n=10,
	terms.to.use=NULL,
	limits=NULL
) {
	leaf.term.pair.marginal.matrix(
		term.descendancy.matrix,
		phi.trace, 
		terms.to.use=terms.to.use,
		n=n
	) %>%
	(function(x) {
		if (use.clinical.names) {
			rownames(x) <- get.shortened.names(hpo.terms, rownames(x))
			colnames(x) <- get.shortened.names(hpo.terms, colnames(x))
			x
		} else {
			x
		}
	}) %>%
	(function(x) {
		x[lower.tri(x, diag=TRUE)] <- NA
		x
	})
}

#' Plot marginal distribution of HPO term pairs
#'
#' Plots the marginal frequency of HPO terms pairs' inclusion in the phi parameter over the course of an application of \code{\link{sim.reg}} as heatmap
#'
#' @template hpo.terms
#' @template term.descendancy.matrix
#' @param phi.trace Delta trace output from \code{\link{sim.reg}} function, where the rows are the values of phi at each iteration
#' @param use.clinical.names Logical indicating whether to use full clinical names for the terms (TRUE) or just the codes (FALSE)
#' @param n Maximum number of terms for plot
#' @param terms.to.use Specify the terms to appear in the plot
#' @param limits Range of frequencies to set colours for
#' @return ggplot2 plot object
two.term.marginals.plot <- function(
	hpo.terms,
	term.descendancy.matrix, 
	phi.trace, 
	use.clinical.names=FALSE,
	n=10,
	terms.to.use=NULL,
	limits=NULL
) {
	get.two.term.marginals(
		hpo.terms,
		term.descendancy.matrix, 
		phi.trace, 
		use.clinical.names,
		n,
		terms.to.use,
		limits
	) %>%
	(function(heatmap.mat) {
		ggplot(
			data=melt(heatmap.mat, varnames=c("Term1", "Term2"), value.name="Marginal.Frequency"),
			mapping=aes_string(x="Term1", y="Term2", fill="Marginal.Frequency")
		) + 
		scale_x_discrete(expand = c(0, 0)) +
		scale_y_discrete(expand = c(0, 0)) + 
		geom_tile() +
		theme(
			axis.text.x = element_text(angle="if"(use.clinical.names, 45, 90), hjust=1),
			panel.background = element_blank()
		) + 
		scale_fill_gradientn(na.value="white", guide="colourbar", limits=limits, colours = colorRampPalette(c("Blue", "green3","Yellow"))( 18 ))
	})
}

#' Create grid of SimReg output plots
#'
#' Create a PDF containing a set of summary plots for the output of the \code{sim.reg} application. Output contains plots of marginal probabilities of terms and pairs of terms being present phi and histograms giving the distributions of variables.
#'
#' @template hpo.terms
#' @template term.descendancy.matrix
#' @param file.name File to write plots to 
#' @param sim.reg.out Output of call to \code{sim.reg}
#' @return Plots graph to file
sim.reg.summary <- function(hpo.terms, term.descendancy.matrix, file.name, sim.reg.out) {
	x <- 1
	colour.points <- function(x) ifelse(x$gamma, "red", "black")

	plot.2d.margs <- function(x, ...) {
		image(x, col=rgb(c(rep(0, 50), seq(from=0, to=1, by=0.02)), seq(from=0, to=1, by=0.01), seq(from=0.5, to=0, by=-0.005)), axes=FALSE, ...)
		axis(1, at=0:(ncol(x)-1)/(ncol(x)-1), las=2, labels=colnames(x))
		axis(2, las=2, labels=rownames(x), at=0:(nrow(x)-1)/(nrow(x)-1)) 
		color.legend(align="rb", gradient="y", rect.col=rgb(c(rep(0, 50), seq(from=0, to=1, by=0.02)), seq(from=0, to=1, by=0.01), seq(from=0.5, to=0, by=-0.005)), xl=1.02, xr=1.06, yb=0, yt=1, legend=as.character(signif(digits=2, range(x, na.rm=TRUE))))
	}

	pbeta2 <- function(x, mean, alpha.plus.beta) pbeta(q=x, shape1=mean*alpha.plus.beta, shape2=(1-mean)*alpha.plus.beta)

	pdf(file.name, width=25, height=35)

	equal.w.cols <- c(3, 9)
	its <- length(sim.reg.out$gamma)

	layout(do.call(what=cbind, lapply(1:length(equal.w.cols), function(col.ind) sort(rep(Reduce(init=1, f="+", x=c(0, equal.w.cols)[1:col.ind]) + 0:(equal.w.cols[col.ind]-1), as.integer(200/equal.w.cols[col.ind])+1)[1:200]))))

	par(mar=c(5, 4, 4, 2))

	variables <- c("gamma", "alpha_star", "alpha", "log_beta", "logit_mean_f", "log_alpha_plus_beta_f", "logit_mean_g", "log_alpha_plus_beta_g", "phi")
	numeric.vars <- Filter(x=variables, f=function(x) class(sim.reg.out[[x]]) == "numeric")
	non.char <- Filter(x=variables, f=function(x) mode(sim.reg.out[[x]]) != "character")
	no.gamma <- setdiff(variables, "gamma")

	plot.new()
	legend(pch=19, x="topright", col=c("red", "black"), legend=expression(paste(gamma, " = 1", sep=""), paste(gamma, " = 0", sep="")))
	par(mar=c(5, 4, 4, 2))
	for (i in 1:length(non.char)) {
		text(cex=2, adj=0, x=0, y=1-i/(2*length(non.char)), paste(non.char[i], ": ", round(digits=2, mean(sim.reg.out[[non.char[i]]])), sep=""))
	}

	for (i in 1:length(no.gamma)) {
		text(cex=2, adj=0, x=0, y=0.5-i/(2*length(no.gamma)), paste(no.gamma[i], " acceptance: ", round(digits=2, mean(sim.reg.out[[paste(no.gamma[i], "_accept", sep="")]])), sep=""))
	}
	par(mar=c(5, 4, 4, 2))

	hpo.plot.marginal.freqs(hpo.terms, term.descendancy.matrix, sim.reg.out$phi[sim.reg.out$gamma,], font.size=55, main=expression(paste("Marginal Posterior on Single Terms in ", phi, sep="")), max.terms=10, colour.gradient=FALSE, colours="sky blue")	
	par(mar=c(15, 15, 3, 3))

	if (max(apply(as.row.leaves(term.descendancy.matrix, sim.reg.out$phi), 1, sum)) > 1) {
		get.two.term.marginals(use.clinical.names=TRUE, hpo.terms=hpo.terms, term.descendancy.matrix=term.descendancy.matrix, phi.trace=sim.reg.out$phi[sim.reg.out$gamma,], n=10) %>% plot.2d.margs(main=expression(paste("Marginal Posterior on Term Pairs in ", phi, sep="")))
	} else {
		plot(1)
	}

	par(mar=c(5, 4, 4, 2))

	alpha_star.range <- range(sim.reg.out$alpha_star, qnorm(p=c(0.05, 1-0.05), mean=sim.reg.out$priors$alpha_star_mean, sd=sim.reg.out$priors$alpha_star_sd))
	hist(main=NULL, x=sim.reg.out$alpha_star[sim.reg.out$gamma], ylab="Density", xlab="log beta", freq=FALSE, breaks=seq(from=alpha_star.range[1], to=alpha_star.range[2], by=diff(alpha_star.range)/100), xlim=alpha_star.range)
	curve(add=TRUE, expr=dnorm(x, mean=sim.reg.out$priors$alpha_star_mean, sd=sim.reg.out$priors$alpha_star_sd), col="blue")

	alpha.range <- range(sim.reg.out$alpha, qnorm(p=c(0.05, 1-0.05), mean=sim.reg.out$priors$alpha_mean, sd=sim.reg.out$priors$alpha_sd))
	hist(main=NULL, x=sim.reg.out$alpha[sim.reg.out$gamma], ylab="Density", xlab="log beta", freq=FALSE, breaks=seq(from=alpha.range[1], to=alpha.range[2], by=diff(alpha.range)/100), xlim=alpha.range)
	curve(add=TRUE, expr=dnorm(x, mean=sim.reg.out$priors$alpha_mean, sd=sim.reg.out$priors$alpha_sd), col="blue")

	log_beta.range <- range(sim.reg.out$log_beta, qnorm(p=c(0.05, 1-0.05), mean=sim.reg.out$priors$log_beta_mean, sd=sim.reg.out$priors$log_beta_sd))
	hist(main=NULL, x=sim.reg.out$log_beta[sim.reg.out$gamma], ylab="Density", xlab="log beta", freq=FALSE, breaks=seq(from=log_beta.range[1], to=log_beta.range[2], by=diff(log_beta.range)/100), xlim=log_beta.range)
	curve(add=TRUE, expr=dnorm(x, mean=sim.reg.out$priors$log_beta_mean, sd=sim.reg.out$priors$log_beta_sd), col="blue")

	hist(main=NULL, x=1-1/(1+exp(sim.reg.out$logit_mean_f[sim.reg.out$gamma])), xlim=c(0, 1), breaks=seq(from=0, to=1, by=0.01), freq=FALSE, xlab="mean f")
	curve(add=TRUE, col="blue", expr=dnorm(log(x)-log(1-x), mean=sim.reg.out$priors$logit_mean_f_mean, sd=sim.reg.out$priors$logit_mean_f_sd)/x/(1-x))

	log_alpha_plus_beta_f.range <- range(sim.reg.out$log_alpha_plus_beta_f[sim.reg.out$gamma], qnorm(p=c(0.05, 1-0.05), mean=sim.reg.out$priors$log_alpha_plus_beta_f_mean, sd=sim.reg.out$priors$log_alpha_plus_beta_f_sd))
	hist(main=NULL, x=sim.reg.out$log_alpha_plus_beta_f[sim.reg.out$gamma], breaks=seq(from=log_alpha_plus_beta_f.range[1], to=log_alpha_plus_beta_f.range[2], by=diff(log_alpha_plus_beta_f.range)/100), freq=FALSE, xlab="log_alpha_plus_beta_f")
	curve(add=TRUE, expr=dnorm(x, mean=sim.reg.out$priors$log_alpha_plus_beta_f_mean, sd=sim.reg.out$priors$log_alpha_plus_beta_f_sd), col="blue")

	plot(main="f transformation sample", x=NULL, xlim=c(0, 1), ylim=c(0, 1))
	for (i in sample(which(sim.reg.out$gamma), size=min(sum(sim.reg.out$gamma), 100))) curve(add=T, col=rgb(0,0,1,0.3), expr=pbeta2(x, mean=1-1/(1+exp(sim.reg.out$logit_mean_f[i])), alpha.plus.beta=exp(sim.reg.out$log_alpha_plus_beta_f[i])))

	hist(main=NULL, x=1-1/(1+exp(sim.reg.out$logit_mean_g[sim.reg.out$gamma])), xlim=c(0, 1), breaks=seq(from=0, to=1, by=0.01), freq=FALSE, xlab="mean g")
	curve(add=TRUE, col="blue", expr=dnorm(log(x)-log(1-x))/x/(1-x))

	log_alpha_plus_beta_g.range <- range(sim.reg.out$log_alpha_plus_beta_g[sim.reg.out$gamma], qnorm(p=c(0.05, 1-0.05), mean=sim.reg.out$priors$log_alpha_plus_beta_g_mean, sd=sim.reg.out$priors$log_alpha_plus_beta_g_sd))
	hist(main=NULL, x=sim.reg.out$log_alpha_plus_beta_g[sim.reg.out$gamma], breaks=seq(from=log_alpha_plus_beta_g.range[1], to=log_alpha_plus_beta_g.range[2], by=diff(log_alpha_plus_beta_g.range)/100), freq=FALSE, xlab="log_alpha_plus_beta_g")
	curve(add=TRUE, expr=dnorm(x, mean=sim.reg.out$priors$log_alpha_plus_beta_g_mean, sd=sim.reg.out$priors$log_alpha_plus_beta_g_sd), col="blue")

	plot(main="g transformation sample", x=NULL, xlim=c(0, 1), ylim=c(0, 1))
	for (i in sample(which(sim.reg.out$gamma), size=min(sum(sim.reg.out$gamma), 100))) curve(add=T, col=rgb(0,0,1,0.3), expr=pbeta2(x, mean=1-1/(1+exp(sim.reg.out$logit_mean_g[i])), alpha.plus.beta=exp(sim.reg.out$log_alpha_plus_beta_g[i])))

	dev.off()
}

#' Create grid of SimReg output plots
#'
#' Create a PDF containing a comprehensive set of summary plots for the output of the \code{sim.reg} application. Output contains plots of marginal probabilities of terms and pairs of terms being present phi, traces of variables and their log likelihoods recorded during the MCMC application, and histograms giving the distributions of variables.
#'
#' @template hpo.terms
#' @template term.descendancy.matrix
#' @param file.name File to write plots to 
#' @param sim.reg.out Output of call to \code{sim.reg}
#' @return Plots graph to file
sim.reg.full.summary <- function(hpo.terms, term.descendancy.matrix, file.name, sim.reg.out) {
	x <- 1
	colour.points <- function(x) ifelse(x$gamma, "red", "black")

	plot.2d.margs <- function(x, ...) {
		image(x, col=rgb(c(rep(0, 50), seq(from=0, to=1, by=0.02)), seq(from=0, to=1, by=0.01), seq(from=0.5, to=0, by=-0.005)), axes=FALSE, ...)
		axis(1, at=0:(ncol(x)-1)/(ncol(x)-1), las=2, labels=colnames(x))
		axis(2, las=2, labels=rownames(x), at=0:(nrow(x)-1)/(nrow(x)-1)) 
		color.legend(align="rb", gradient="y", rect.col=rgb(c(rep(0, 50), seq(from=0, to=1, by=0.02)), seq(from=0, to=1, by=0.01), seq(from=0.5, to=0, by=-0.005)), xl=1.02, xr=1.06, yb=0, yt=1, legend=as.character(signif(digits=2, range(x, na.rm=TRUE))))
	}

	pbeta2 <- function(x, mean, alpha.plus.beta) pbeta(q=x, shape1=mean*alpha.plus.beta, shape2=(1-mean)*alpha.plus.beta)

	pdf(file.name, width=25, height=35)

	equal.w.cols <- c(4, 16, 11)
	its <- length(sim.reg.out$gamma)

	layout(do.call(what=cbind, lapply(1:length(equal.w.cols), function(col.ind) sort(rep(Reduce(init=1, f="+", x=c(0, equal.w.cols)[1:col.ind]) + 0:(equal.w.cols[col.ind]-1), as.integer(200/equal.w.cols[col.ind])+1)[1:200]))))

	par(mar=c(5, 4, 4, 2))

	variables <- c("gamma", "alpha_star", "alpha", "log_beta", "logit_mean_f", "log_alpha_plus_beta_f", "logit_mean_g", "log_alpha_plus_beta_g", "phi")
	numeric.vars <- Filter(x=variables, f=function(x) class(sim.reg.out[[x]]) == "numeric")
	non.char <- Filter(x=variables, f=function(x) mode(sim.reg.out[[x]]) != "character")
	no.gamma <- setdiff(variables, "gamma")

	plot.new()
	legend(pch=19, x="topright", col=c("red", "black"), legend=expression(paste(gamma, " = 1", sep=""), paste(gamma, " = 0", sep="")))
	par(mar=c(5, 4, 4, 2))
	for (i in 1:length(non.char)) {
		text(cex=2, adj=0, x=0, y=1-i/length(non.char), paste(non.char[i], ": ", round(digits=2, mean(sim.reg.out[[non.char[i]]])), sep=""))
	}

	plot.new()
	for (i in 1:length(no.gamma)) {
		text(cex=2, adj=0, x=0, y=1-i/length(no.gamma), paste(no.gamma[i], " acceptance: ", round(digits=2, mean(sim.reg.out[[paste(no.gamma[i], "_accept", sep="")]])), sep=""))
	}
	par(mar=c(5, 4, 4, 2))

	hpo.plot.marginal.freqs(hpo.terms, term.descendancy.matrix, sim.reg.out$phi[sim.reg.out$gamma,], font.size=55, main=expression(paste("Marginal Posterior on Single Terms in ", phi, sep="")), max.terms=10, colour.gradient=FALSE, colours="sky blue")	
	par(mar=c(15, 15, 3, 3))

	if (max(apply(as.row.leaves(term.descendancy.matrix, sim.reg.out$phi), 1, sum)) > 1) {
		get.two.term.marginals(use.clinical.names=TRUE, hpo.terms=hpo.terms, term.descendancy.matrix=term.descendancy.matrix, phi.trace=sim.reg.out$phi[sim.reg.out$gamma,], n=10) %>% plot.2d.margs(main=expression(paste("Marginal Posterior on Term Pairs in ", phi, sep="")))
	} else {
		plot(1)
	}

	par(mar=c(5, 4, 4, 2))

	for (var.name in numeric.vars) {
		plot(ylab=var.name, xlab="update cycle", main=paste(var.name, " trace", sep=""), y=sim.reg.out[[var.name]][seq(from=1, to=its, by=max(1, as.integer(its/2000)))], x=seq(from=1, to=its, by=max(1, as.integer(its/2000))), col=colour.points(sim.reg.out)[seq(from=1, to=its, by=max(1, as.integer(its/2000)))])
	}

	for (var.name in c("likelihood", "alpha_star", "alpha", "log_beta", "phi", "logit_mean_f", "log_alpha_plus_beta_f", "logit_mean_g", "log_alpha_plus_beta_g")) {
		plot(xlab="update cycle", ylab="log likelihood", main=paste(if (var.name == "likelihood") "total" else var.name, " log likelihood trace", sep=""), x=seq(from=1, to=its, by=max(1, as.integer(its/2000))), y=sim.reg.out$liks[[var.name]][seq(from=1, to=its, by=max(1, as.integer(its/2000)))], col=colour.points(sim.reg.out)[seq(from=1, to=its, by=max(1, as.integer(its/2000)))]) }

	alpha_star.range <- range(sim.reg.out$alpha_star, qnorm(p=c(0.05, 1-0.05), mean=sim.reg.out$priors$alpha_star_mean, sd=sim.reg.out$priors$alpha_star_sd))
	hist(main=NULL, x=sim.reg.out$alpha_star[sim.reg.out$gamma], ylab="Density", xlab="log beta", freq=FALSE, breaks=seq(from=alpha_star.range[1], to=alpha_star.range[2], by=diff(alpha_star.range)/100), xlim=alpha_star.range)
	curve(add=TRUE, expr=dnorm(x, mean=sim.reg.out$priors$alpha_star_mean, sd=sim.reg.out$priors$alpha_star_sd), col="blue")

	alpha.range <- range(sim.reg.out$alpha, qnorm(p=c(0.05, 1-0.05), mean=sim.reg.out$priors$alpha_mean, sd=sim.reg.out$priors$alpha_sd))
	hist(main=NULL, x=sim.reg.out$alpha[sim.reg.out$gamma], ylab="Density", xlab="log beta", freq=FALSE, breaks=seq(from=alpha.range[1], to=alpha.range[2], by=diff(alpha.range)/100), xlim=alpha.range)
	curve(add=TRUE, expr=dnorm(x, mean=sim.reg.out$priors$alpha_mean, sd=sim.reg.out$priors$alpha_sd), col="blue")

	log_beta.range <- range(sim.reg.out$log_beta, qnorm(p=c(0.05, 1-0.05), mean=sim.reg.out$priors$log_beta_mean, sd=sim.reg.out$priors$log_beta_sd))
	hist(main=NULL, x=sim.reg.out$log_beta[sim.reg.out$gamma], ylab="Density", xlab="log beta", freq=FALSE, breaks=seq(from=log_beta.range[1], to=log_beta.range[2], by=diff(log_beta.range)/100), xlim=log_beta.range)
	curve(add=TRUE, expr=dnorm(x, mean=sim.reg.out$priors$log_beta_mean, sd=sim.reg.out$priors$log_beta_sd), col="blue")

	hist(main=NULL, x=1-1/(1+exp(sim.reg.out$logit_mean_f[sim.reg.out$gamma])), xlim=c(0, 1), breaks=seq(from=0, to=1, by=0.01), freq=FALSE, xlab="mean f")
	curve(add=TRUE, col="blue", expr=dnorm(log(x)-log(1-x), mean=sim.reg.out$priors$logit_mean_f_mean, sd=sim.reg.out$priors$logit_mean_f_sd)/x/(1-x))

	log_alpha_plus_beta_f.range <- range(sim.reg.out$log_alpha_plus_beta_f[sim.reg.out$gamma], qnorm(p=c(0.05, 1-0.05), mean=sim.reg.out$priors$log_alpha_plus_beta_f_mean, sd=sim.reg.out$priors$log_alpha_plus_beta_f_sd))
	hist(main=NULL, x=sim.reg.out$log_alpha_plus_beta_f[sim.reg.out$gamma], breaks=seq(from=log_alpha_plus_beta_f.range[1], to=log_alpha_plus_beta_f.range[2], by=diff(log_alpha_plus_beta_f.range)/100), freq=FALSE, xlab="log_alpha_plus_beta_f")
	curve(add=TRUE, expr=dnorm(x, mean=sim.reg.out$priors$log_alpha_plus_beta_f_mean, sd=sim.reg.out$priors$log_alpha_plus_beta_f_sd), col="blue")

	plot(main="mean of f against slope of g", 1-1/(1+exp(sim.reg.out$logit_mean_f[sim.reg.out$gamma])), sim.reg.out$log_alpha_plus_beta_f[sim.reg.out$gamma], xlab="mean f", ylab="log alpha plus beta f")

	plot(main="f transformation sample", x=NULL, xlim=c(0, 1), ylim=c(0, 1))
	for (i in sample(which(sim.reg.out$gamma), size=min(sum(sim.reg.out$gamma), 100))) curve(add=T, col=rgb(0,0,1,0.3), expr=pbeta2(x, mean=1-1/(1+exp(sim.reg.out$logit_mean_f[i])), alpha.plus.beta=exp(sim.reg.out$log_alpha_plus_beta_f[i])))

	hist(main=NULL, x=1-1/(1+exp(sim.reg.out$logit_mean_g[sim.reg.out$gamma])), xlim=c(0, 1), breaks=seq(from=0, to=1, by=0.01), freq=FALSE, xlab="mean g")
	curve(add=TRUE, col="blue", expr=dnorm(log(x)-log(1-x))/x/(1-x))

	log_alpha_plus_beta_g.range <- range(sim.reg.out$log_alpha_plus_beta_g[sim.reg.out$gamma], qnorm(p=c(0.05, 1-0.05), mean=sim.reg.out$priors$log_alpha_plus_beta_g_mean, sd=sim.reg.out$priors$log_alpha_plus_beta_g_sd))
	hist(main=NULL, x=sim.reg.out$log_alpha_plus_beta_g[sim.reg.out$gamma], breaks=seq(from=log_alpha_plus_beta_g.range[1], to=log_alpha_plus_beta_g.range[2], by=diff(log_alpha_plus_beta_g.range)/100), freq=FALSE, xlab="log_alpha_plus_beta_g")
	curve(add=TRUE, expr=dnorm(x, mean=sim.reg.out$priors$log_alpha_plus_beta_g_mean, sd=sim.reg.out$priors$log_alpha_plus_beta_g_sd), col="blue")

	plot(main="mean of g against slope of g",1-1/(1+exp(sim.reg.out$logit_mean_f[sim.reg.out$gamma])), sim.reg.out$log_alpha_plus_beta_f[sim.reg.out$gamma], xlab="mean g", ylab="log alpha plus beta g")

	plot(main="g transformation sample", x=NULL, xlim=c(0, 1), ylim=c(0, 1))
	for (i in sample(which(sim.reg.out$gamma), size=min(sum(sim.reg.out$gamma), 100))) curve(add=T, col=rgb(0,0,1,0.3), expr=pbeta2(x, mean=1-1/(1+exp(sim.reg.out$logit_mean_g[i])), alpha.plus.beta=exp(sim.reg.out$log_alpha_plus_beta_g[i])))

	dev.off()
}


