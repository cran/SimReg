## ----echo=FALSE----------------------------------------------------------
knitr::opts_chunk$set(dev="svg", fig.width=7, fig.height=7)

## ------------------------------------------------------------------------
library(ontologyIndex)
library(SimReg)
data(hpo)
set.seed(1)

template <- c("HP:0005537", "HP:0000729", "HP:0001873")
terms <- get_ancestors(hpo, c(template, sample(hpo$id, size=50)))

## ------------------------------------------------------------------------
y <- c(rep(TRUE, 10), rep(FALSE, 90))
x <- replicate(simplify=FALSE, n=100, expr=minimal_set(hpo, sample(terms, size=5)))

## ------------------------------------------------------------------------
y
head(x)

## ------------------------------------------------------------------------
no_assoc <- sim_reg(ontology=hpo, x=x, y=y)
prob_association(no_assoc)

## ------------------------------------------------------------------------
x_assoc <- lapply(y, function(y_i) minimal_set(hpo, c(
	sample(terms, size=5), if (y_i) sample(template, size=2))))

## ------------------------------------------------------------------------
head(x_assoc)

## ------------------------------------------------------------------------
assoc <- sim_reg(ontology=hpo, x=x_assoc, y=y)
prob_association(assoc)

## ------------------------------------------------------------------------
plot_term_marginals(hpo, get_term_marginals(assoc), max_terms=8, fontsize=30)

