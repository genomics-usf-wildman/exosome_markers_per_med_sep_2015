library(data.table)
library(reshape2)

args <- c("gene_grouping_variables",
          "isoform_grouping_variables",
          "tissue_specific_markers")

args <- commandArgs(trailingOnly=TRUE)

load(args[1])
load(args[2])
gc()

### combine sample name and srx into the gene counts file
.narm.mean <- function(x){
    if(length(x)==0) {return(0)} else {return(mean(x,na.rm=TRUE))}
}
## this is the Tissue Specificity Index (eq 1) from Yanai et al.
tissue.specificity.index <- function(expression){
    expression <- expression[!is.na(expression)]
    expression <- expression[is.finite(expression)]
    if (length(expression)==0) {
        return(0)
    }
    if (max(expression)==0) {
        return(0)
    }
    return(sum(1-expression/max(expression))/
               (length(expression)-1))
}

##' Calculate genes or isoforms which are specific to a tissue
##'
##' Identifies genes or isoforms which are specific to certain tissue.
##' @title specific.genes.isoforms
##' @param data -- wide table of data
##' @param min.specificity -- minimum specificity score, defaults to 0.98
##' @param min.max.expression -- minimum max expression, defaults to 10
##' @return vector of length equal to rows of data which contains the
##' specific tissue that this isoform/gene describes, or NA if it does
##' not describe a specific tissue
##' @author Don Armstrong
specific.genes.isoforms <- function(grouping,min.specificity=0.98,min.max.expression=10) {
    ## calculate the specificity index
    tissue.specific.gi <- copy(grouping)
    tissue.specific.gi[,possibly.specific:=FALSE]
    tissue.specific.gi[tissue.specificity.index >= min.specificity &
                       max >= min.max.expression,possibly.specific:=TRUE]
    tissue.specific.gi[,tissue:=as.character(NA)]
    tissue.specific.gi[possibly.specific==TRUE,
                       tissue:=tissue_max]
    tissue.specific.gi <-
        tissue.specific.gi[order(tissue,-max),
                           list(gene=gene_short_name,
                                tracking_id=tracking_id,
                                tissue=tissue,
                                gi.max=max,
                                specificity=tissue.specificity.index,
                                log.specificity=tissue.specificity.index.log2,
                                possibly.specific=possibly.specific,
                                gi.which.max=tissue_max,
                                mean=mean,
                                mean.log2=mean.log2,
                                var=var,
                                var.log2=var.log2,
                                entropy=entropy,
                                entropy.log2=entropy.log2)]
    return(tissue.specific.gi)
}

tissue.specific.genes <-
    specific.genes.isoforms(gene.grouping.variables)
tissue.specific.isoforms <-
    specific.genes.isoforms(isoform.grouping.variables)

save(tissue.specific.genes,
     tissue.specific.isoforms,
     file=args[length(args)])

