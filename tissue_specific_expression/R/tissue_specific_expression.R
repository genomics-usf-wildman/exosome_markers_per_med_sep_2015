library(data.table)
library(reshape2)

args <- c("categorized_samples",
          "combined_read_counts",
          "tissue_specific_markers")

args <- commandArgs(trailingOnly=TRUE)

load(args[1])
load(args[2])

setkey(categorized.samples,"SRX")
setnames(gene.counts,"srx","SRX")
setkey(gene.counts,"SRX")
setnames(isoform.counts,"srx","SRX")
setkey(isoform.counts,"SRX")

### combine sample name and srx into the gene counts file
combined.gene.reads <-
    categorized.samples[,list(Sample_Group,SRX,GEO_Accession)][gene.counts]

### do the same for the isoforms
combined.isoform.reads <-
    categorized.samples[,list(Sample_Group,SRX,GEO_Accession)][isoform.counts]

.narm.mean <- function(x){
    mean(x,na.rm=TRUE)
}

c.gene.reads.wide <-
    dcast(combined.gene.reads,
          gene_short_name~Sample_Group,
          fun.aggregate=.narm.mean,
          value.var="FPKM"
          )

c.isoform.reads.wide <-
    dcast(combined.isoform.reads,
          tracking_id~Sample_Group,
          fun.aggregate=.narm.mean,
          value.var="FPKM")

## this is the Tissue Specificity Index (eq 1) from Yanai et al.
tissue.specificity.index <- function(expression){
    return(sum(1-expression/max(expression))/
               (length(expression)-1))
}
##' Calculate genes or isoforms which are specific to a tissue
##'
##' Identifies genes or isoforms which are specific to certain tissue.
##' @title specific.genes.isoforms
##' @param data -- wide table of data
##' @param min.specificity -- minimum specificity score, defaults to 0.99
##' @param min.max.expression -- minimum max expression, defaults to 10
##' @return vector of length equal to rows of data which contains the
##' specific tissue that this isoform/gene describes, or NA if it does
##' not describe a specific tissue
##' @author Don Armstrong
specific.genes.isoforms <- function(data,min.specificity=0.99,min.max.expression=10) {
    ## calculate the specificity index
    specificity <-
        apply(data[,-1],1,tissue.specificity.index)
    per.gi.max <-
        apply(data[,-1],1,max)
    per.gi.which.max <-
        apply(data[,-1],1,which.max)
    possibly.specific <-
        (per.gi.max >= min.max.expression &
             specificity >= min.specificity)
    possibly.specific[is.na(possibly.specific)] <-
        FALSE
    tissue.specific.gi <-
        data.frame(gene=data[,1],
                   tissue=NA,
                   gi.max=per.gi.max,
                   specificity=specificity,
                   possibly.specific=possibly.specific,
                   gi.which.max=per.gi.which.max)
    tissue.specific.gi$tissue[possibly.specific] <- 
        colnames(data)[per.gi.which.max+1][possibly.specific]
    return(tissue.specific.gi)
}    

tissue.specific.genes <- specific.genes.isoforms(c.gene.reads.wide)
tissue.specific.isoforms <- specific.genes.isoforms(c.isoform.reads.wide)


save(tissue.specific.genes,
     tissue.specific.isoforms,
     file=args[length(args)])

