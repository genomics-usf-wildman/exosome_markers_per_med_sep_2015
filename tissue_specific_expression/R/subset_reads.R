library(data.table)
library(reshape2)

args <- c("categorized_samples",
          "combined_read_counts",
          "interesting_gene_reads",
          "interesting_isoform_reads"
          )

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
    categorized.samples[,list(Sample_Group,SRX)][gene.counts]

### do the same for the isoforms
combined.isoform.reads <-
    categorized.samples[,list(Sample_Group,SRX)][isoform.counts]

.narm.mean <- function(x){
    if(length(x)==0) {return(0)} else {return(mean(x,na.rm=TRUE))}
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

calculate.entropy <- function(x,bins=10,method="MM"){
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    result <- 0
    try(result <- entropy(discretize(log(x+1),numBins=bins),method=method),
        silent=TRUE)
    return(result)
}

gene.entropy <-
    dcast(combined.gene.reads,
          gene_short_name~.,
          fun.aggregate=calculate.entropy,
          value.var="FPKM")
colnames(gene.entropy)[2] <- "entropy"

gene.max <-
    dcast(combined.gene.reads,
          gene_short_name~.,
          fun.aggregate=function(x){max(x,na.rm=TRUE)},
          value.var="FPKM")
colnames(gene.max)[2] <- "max"

gene.var <-
    dcast(combined.gene.reads,
          gene_short_name~.,
          fun.aggregate=function(x){var(x,na.rm=TRUE)},
          value.var="FPKM")
colnames(gene.variance)[2] <- "var"

isoform.entropy <-
    dcast(combined.isoform.reads,
          tracking_id~.,
          fun.aggregate=calculate.entropy,
          value.var="FPKM")
colnames(isoform.entropy)[2] <- "entropy"

isoform.max <-
    dcast(combined.isoform.reads,
          tracking_id~.,
          fun.aggregate=function(x){max(x,na.rm=TRUE)},
          value.var="FPKM")
colnames(isoform.max)[2] <- "max"

isoform.var <-
    dcast(combined.isoform.reads,
          tracking_id~.,
          fun.aggregate=function(x){var(x,na.rm=TRUE)},
          value.var="FPKM")
colnames(isoform.max)[2] <- "var"



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

min.entropy <- 1.75
min.tissue.specificity <- 0.98
min.expression <- 3
min.var <- 0.01

gene.tissue.specificity <-
    data.frame(gene_short_name=c.gene.reads.wide[,1],
               tissue.specificity.index=apply(c.gene.reads.wide[,-1],1,
                   tissue.specificity.index))

isoform.tissue.specificity <-
    data.frame(tracking_id=c.isoform.reads.wide[,1],
               tissue.specificity.index=apply(c.isoform.reads.wide[,-1],1,
                   tissue.specificity.index))

interesting.genes <-
    gene.entropy[(gene.entropy[,2] >= min.entropy |
                      gene.tissue.specificity[,2] >= min.tissue.specificity)
                 & gene.max[,2] >= min.expression
                 & gene.var[,2] >= min.var,
                 1]

interesting.isoforms <-
    isoform.entropy[(isoform.entropy[,2] >= min.entropy |
                      isoform.tissue.specificity[,2] >= min.tissue.specificity)
                 & isoform.max[,2] >= min.expression
                 & isoform.var[,2] >= min.var,
                 1]
setkey(combined.gene.reads,"gene_short_name")
interesting.gene.reads <-
    combined.gene.reads[interesting.genes,list(Sample_Group,SRX,gene_id,
                                               gene_short_name,
                                               FPKM,FPKM_conf_lo,FPKM_conf_hi,
                                               FPKM_status)]

setkey(combined.isoform.reads,"tracking_id")

interesting.isoform.reads <-
    combined.isoform.reads[interesting.isoforms,
                           list(Sample_Group,SRX,
                                tracking_id,
                                gene_id,
                                gene_short_name,
                                FPKM,FPKM_conf_lo,FPKM_conf_hi,
                                FPKM_status)
                           ]

save(interesting.gene.reads,
     gene.entropy,
     gene.tissue.specificity,
     file=args[3])

save(interesting.isoform.reads,
     isoform.entropy,
     isoform.tissue.specificity,
     file=args[4])

     
