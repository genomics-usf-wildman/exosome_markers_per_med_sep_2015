library("data.table")
library("reshape2")
library("yaml")

args <- c("categorized_samples",
          "combined_read_counts",
          "brain_tissues.yaml",
          "housekeeping_genes_superset",
          "brain_specific_markers")

args <- commandArgs(trailingOnly=TRUE)

### we cannot use the interesting genes subset, because we are
### generating new tissue super-sets here.
load(args[1])
load(args[2])
### housekeeping_genes
load(args[4])

hk.g.tracking <- housekeeping.genes.superset[,tracking_id]

setkey(categorized.samples,"SRX")

grouping.by.sample.formula <-
    as.formula(tracking_id~Sample_Group)
grouping.formula <-
    as.formula(tracking_id~.)

setnames(gene.counts,"srx","SRX")
setkey(gene.counts,"SRX")

### combine sample name and srx into the gene counts file
combined.reads <-
    categorized.samples[,list(Sample_Group,SRX)][gene.counts]

.narm.mean <- function(x){
    if(length(x)==0) {return(0)} else {return(mean(x,na.rm=TRUE))}
}

gene.reads.wide <-
    data.table(dcast(combined.reads,
                     gene_short_name+tracking_id~Sample_Group,
                     fun.aggregate=.narm.mean,
                     value.var="FPKM"
                     ))

### load the list of brain tissues from the yaml file
brain.tissues <-
    data.table(t(sapply(yaml.load_file(args[3])[[1]],
                        unlist)))
    

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
specific.genes.isoforms <- function(data,tissues,tissue.name,min.specificity=0.98,min.max.expression=10) {
    ## calculate the specificity index
    genes <- data[,gene_short_name]
    tracking.id <- data[,tracking_id]
    data <- data[,-c(1:2),with=FALSE]
    if (any(colnames(data)=="universal human reference")) {
        data <- data[,colnames(data)!="universal human reference",with=FALSE]
    }
    ## calculate the tissue.name to be equal to the maximum expression
    data[,eval(parse(text=paste0(tissue.name,
                         ":=apply(data[,tissues,with=FALSE],1,max)"
                                 )))]
    data[,(tissues[!tissues %in% tissue.name]):=NULL]
        
    specificity <-
        apply(data,1,tissue.specificity.index)
    per.gi.max <-
        apply(data,1,max)
    per.gi.which.max <-
        apply(data,1,which.max)
    possibly.specific <-
        (per.gi.max >= min.max.expression &
             specificity >= min.specificity)
    possibly.specific[is.na(possibly.specific)] <-
        FALSE
    tissue.specific.gi <-
        data.frame(gene=genes,
                   tracking_id=tracking.id,
                   tissue=NA,
                   gi.max=per.gi.max,
                   specificity=specificity,
                   possibly.specific=possibly.specific,
                   gi.which.max=per.gi.which.max)

    tissue.specific.gi$tissue[possibly.specific] <- 
        colnames(data)[per.gi.which.max][possibly.specific]
    return(data.table(tissue.specific.gi[!is.na(tissue.specific.gi$tissue) &
                                             tissue.specific.gi$tissue == tissue.name,]))
}

brain.tissue.specific.genes <- list()
for (tissue in colnames(brain.tissues)[-1]) {
    tissues <- brain.tissues[eval(parse(text=paste0(tissue,"==TRUE"))),name]
    brain.tissue.specific.genes[[tissue]] <-
        specific.genes.isoforms(gene.reads.wide,
                                tissues=tissues,
                                tissue.name=tissue)
}
brain.tissue.specific.genes <-
    rbindlist(brain.tissue.specific.genes)

### identify genes which are highly expressed in brain tissues even if
### they are not specific

brain.reads <-
    gene.reads.wide[,c("gene_short_name","tracking_id",brain.tissues[,name]),
                    with=FALSE]
brain.reads[,max:=apply(brain.reads[,-(1:2),with=FALSE],1,max)]

brain.highly.expressed.genes <-
    gene.reads.wide[sort(unique(as.vector(apply(gene.reads.wide[,brain.tissues[,name],
                                                                with=FALSE],
                                                2,
                                                function(x){order(-x)[1:500]})))),
                    c("gene_short_name","tracking_id",brain.tissues[,name]),with=FALSE]

brain.highly.expressed.genes[,max:=apply(brain.highly.expressed.genes[,-(1:2),with=FALSE],1,max)]


brain.highly.expressed.genes <-
    brain.highly.expressed.genes[!(tracking_id %in% hk.g.tracking),]

save(brain.tissue.specific.genes,
     brain.highly.expressed.genes,
     file=args[length(args)])

