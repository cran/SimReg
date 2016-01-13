cat("Loading up...\n")
suppressPackageStartupMessages(library(SimReg))
data(hpo.terms)

set.seed(0)

local({
	gamma.prior <- 0.05
	N <- 1000
	background.case.mean.terms <- 8
	y1.mean.noise.terms <- 5
	num.y1.cases <- c(1:5 * 2, 20)
	expressivity <- 1:3 / 3
	disease.template.terms <- c("HP:0005537", "HP:0000729", "HP:0001873")
	all.terms <- Filter(x=get.ancestors(hpo.terms, c(disease.template.terms, sample(hpo.terms$id, size=200))), f=function(tm) "HP:0000001" %in% hpo.terms$ancestors[[tm]])
	info.content <- sapply(
		setNames(nm=all.terms),
		function(term) {
			inc.descs <- all.terms[sapply(all.terms, function(tm) term %in% hpo.terms$ancestors[[tm]])]

			-log(1-((length(all.terms) - length(inc.descs))/length(all.terms))^(background.case.mean.terms + 1))
		}
	)

	ttsm <- get.lin.term.term.matrix(
		hpo.terms, 
		sim.reg.sim.pure$config$info.content,
		prune=TRUE
	)

	term.desc.mat <- get.term.descendancy.matrix(hpo.terms, sim.reg.sim.pure$config$all.terms)

	simulation.replicates <- 50
	get.background.case <- function(mean.terms, all.terms) {
		clean.terms(hpo.terms, sample(all.terms, size=1+rpois(lambda=mean.terms, n=1)))
	}
	sim.reg.sim.pure$config$get.gamma1.case <- list(
		independent=function(mean.noise.terms, noise.term.dist, all.terms, key.term.freqs) {
			clean.terms(
				hpo.terms, 
				c(
					names(key.term.freqs)[runif(n=length(key.term.freqs)) < key.term.freqs], 
					noise.term.dist(mean.noise.terms, all.terms=all.terms)
				)
			)	
		},
		correlated=function(mean.noise.terms, noise.term.dist, all.terms, phenos, has.pheno.prob) {
			clean.terms(
				hpo.terms, 
				c(
					phenos[runif(n=1) < has.pheno.prob], 
					noise.term.dist(mean.noise.terms, all.terms=all.terms)
				)
			)	
		}
	)
}

cat("Performing simulations...\n")
print(system.time(
sim.reg.sim.pure$results <- lapply(
	setNames(nm=sim.reg.sim.pure$config$num.y1.cases),
	function(num.y1s) lapply(
		setNames(nm=sim.reg.sim.pure$config$gamma1.term.probabilities),
		function(gamma1.term.probs) lapply(
			list(
				null=Curry(sim.reg.sim.pure$config$get.background.case, mean.terms=sim.reg.sim.pure$config$background.case.mean.terms, all.terms=sim.reg.sim.pure$config$all.terms),
				gamma1.independent=Curry(sim.reg.sim.pure$config$get.gamma1.case$independent, mean.noise.terms=5, noise.term.dist=sim.reg.sim.pure$config$get.background.case, all.terms=sim.reg.sim.pure$config$all.terms, key.term.freqs=setNames(rep(gamma1.term.probs, length(sim.reg.sim.pure$config$case.terms)), sim.reg.sim.pure$config$case.terms)),
				gamma1.correlated=Curry(sim.reg.sim.pure$config$get.gamma1.case$correlated, mean.noise.terms=5, noise.term.dist=sim.reg.sim.pure$config$get.background.case, all.terms=sim.reg.sim.pure$config$all.terms, phenos=sim.reg.sim.pure$config$case.terms, has.pheno.prob=gamma1.term.probs)
			),
			function(term.dist) mclapply(
				1:sim.reg.sim.pure$config$replicates,
				function(replicate.number) {
					y <- c(rep(FALSE, sim.reg.sim.pure$config$N - num.y1s), rep(TRUE, num.y1s))
					x <- c(replicate(n=sim.reg.sim.pure$config$N - num.y1s, expr=sim.reg.sim.pure$config$get.background.case(sim.reg.sim.pure$config$background.case.mean.terms, sim.reg.sim.pure$config$all.terms), simplify=FALSE), replicate(n=num.y1s, expr=term.dist(), simplify=FALSE))
					result <- chib5.pseudo(
						its=sim.reg.sim.pure$config$chain.length,
						pre.its=sim.reg.sim.pure$config$pre.chain.length,
						thin=sim.reg.sim.pure$config$thin,
						mica_ic=sim.reg.sim.pure$config$ttsm,
						term_ic=sim.reg.sim.pure$config$info.content,
						row_is_column_anc=sim.reg.sim.pure$config$term.desc.mat,
						preferred_jump_group=c(0:(length(sim.reg.sim.pure$config$all.terms)-1), rep(match(unlist(lapply(x[y], get.ancestors, hpo.terms=hpo.terms)), colnames(sim.reg.sim.pure$config$term.desc.mat))-1, 50)),
						y=y,
						x=x,
						m_null_prior=0.95
					)

					mean(result$m)
				},
				mc.cores=sim.reg.sim.pure$config$cores
			) %>% simplify2array
		)
	)
)
))

save(
	sim.reg.sim.pure,
	file=paste(
		"sim.reg.sim.pure-v",
		sim.reg.sim.pure$base.data$version,
		".RData",
		sep=""
	)
)
