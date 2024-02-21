## ----echo=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(fig.width=7, fig.height=7)

## -----------------------------------------------------------------------------
library(ontologyIndex)
library(ontologySimilarity)
library(SimReg)
data(hpo)

# only use terms which are descended from HP:0000001
pa_descs <- intersection_with_descendants(hpo, "HP:0000001", hpo$id)
hpo <- lapply(hpo, "[", pa_descs)
class(hpo) <- "ontology_index"

set.seed(1)

terms <- get_ancestors(hpo, c(hpo$id[match(c("Abnormality of thrombocytes","Hearing abnormality"), 
   hpo$name)], intersection_with_descendants(hpo, "HP:0000118", sample(hpo$id[!hpo$obsolete], size=50))))

## -----------------------------------------------------------------------------
hearing_abnormality <- hpo$id[match("Hearing abnormality", hpo$name)]
genotypes <- c(rep(TRUE, 2), rep(FALSE, 98))

#give all subjects 5 random terms and add 'hearing abnormality' for those with y_i=TRUE
phenotypes <- lapply(genotypes, function(y_i) minimal_set(hpo, c(
if (y_i) hearing_abnormality else character(0), sample(terms, size=5))))

## -----------------------------------------------------------------------------
sim_reg(ontology=hpo, x=phenotypes, y=genotypes)

## -----------------------------------------------------------------------------
term_weights <- ifelse(grepl(x=hpo$name, ignore=TRUE, pattern="hearing"), 10, 1)
names(term_weights) <- hpo$id

## -----------------------------------------------------------------------------
thrombocytes <- hpo$id[match("Abnormality of thrombocytes", hpo$name)]
literature_phenotype <- c(hearing_abnormality, thrombocytes)
info <- get_term_info_content(hpo, phenotypes)

term_weights_resnik <- apply(get_term_set_to_term_sims(
	ontology=hpo, information_content=info, terms=names(info),
	term_sim_method="resnik", term_sets=list(literature_phenotype)), 2, mean)

## -----------------------------------------------------------------------------
sim_reg(
	ontology=hpo,
	x=phenotypes,
	y=genotypes,
	term_weights=term_weights_resnik
)

## ----eval=FALSE---------------------------------------------------------------
#  annotation_df <- read.table(header=FALSE, skip=1, sep="\t",
#  	file="ALL_SOURCES_TYPICAL_FEATURES_genes_to_phenotype.txt", stringsAsFactors=FALSE, quote="")
#  hpo_by_gene <- lapply(split(f=annotation_df[,2], x=annotation_df[,4]),
#  	function(trms) minimal_set(hpo, intersect(trms, hpo$id)))

