library(data.table)
library(reshape2)

args <- c("categorized_samples",
          "combined_read_counts",
          "additional_interesting_genes.txt",
          "interesting_gene_reads"
          )

args <- commandArgs(trailingOnly=TRUE)

load(args[1])
load(args[2])
additional.genes <- fread(args[3])

is.genes <- TRUE
if (any(grepl("isoform",args[length(args)]))) {
    is.genes <- FALSE
}

setkey(categorized.samples,"SRX")


## we're using tracking_id for genes instead of gene_short_name
## because SNORD5 and some other genes have multiple gene ids but also
## share a common name
if (is.genes) {
    read.counts <- gene.counts
} else {
    read.counts <- isoform.counts
}
grouping.by.sample.formula <-
    as.formula(tracking_id~Sample_Group)
grouping.formula <-
    as.formula(tracking_id~.)

additional.tracking.ids <-
    read.counts[!duplicated(tracking_id),
                list(tracking_id,gene_short_name)][gene_short_name %in%
                                                       additional.genes[[1]]][[1]]

setnames(read.counts,"srx","SRX")
setkey(read.counts,"SRX")

### combine sample name and srx into the gene counts file
combined.reads <-
    categorized.samples[,list(Sample_Group,SRX)][read.counts]

.narm.mean <- function(x){
    if(length(x)==0) {return(0)} else {return(mean(x,na.rm=TRUE))}
}

c.reads.wide <-
    dcast(combined.reads,
          grouping.by.sample.formula,
          fun.aggregate=.narm.mean,
          value.var="FPKM"
          )

calculate.entropy <- function(x,bins=10,method="MM"){
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    result <- 0
    try(result <- entropy(discretize(log(x+1),numBins=bins),method=method),
        silent=TRUE)
    return(result)
}

grouping.entropy <-
    dcast(combined.reads,
          grouping.formula,
          fun.aggregate=calculate.entropy,
          value.var="FPKM")
colnames(grouping.entropy)[2] <- "entropy"

grouping.max <-
    dcast(combined.reads,
          grouping.formula,
          fun.aggregate=function(x){max(x,na.rm=TRUE)},
          value.var="FPKM")
colnames(grouping.max)[2] <- "max"

grouping.var <-
    dcast(combined.reads,
          grouping.formula,
          fun.aggregate=function(x){var(x,na.rm=TRUE)},
          value.var="FPKM")
colnames(grouping.var)[2] <- "var"

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

tissue.specificity <-
    data.frame(grouping=c.reads.wide[,1],
               tissue.specificity.index=apply(c.reads.wide[,-1],1,
                   tissue.specificity.index))

interesting.groups <-
    grouping.entropy[((grouping.entropy[,2] >= min.entropy |
                          tissue.specificity[,2] >= min.tissue.specificity)
                      & grouping.max[,2] >= min.expression
                      & grouping.var[,2] >= min.var) |
                        grouping.entropy[,1] %in% additional.tracking.ids,
                     1]
if(is.genes) {
    setkey(combined.reads,"tracking_id") 
} else {
    setkey(combined.reads,"tracking_id")
}

interesting.reads <-
    combined.reads[interesting.groups,
                   list(Sample_Group,
                        SRX,
                        tracking_id,
                        gene_id,
                        gene_short_name,
                        FPKM,
                        FPKM_conf_lo,
                        FPKM_conf_hi,
                        FPKM_status)]
if (is.genes) {
    interesting.gene.reads <- interesting.reads
    gene.entropy <- grouping.entropy
    gene.tissue.specificity <- tissue.specificity
    save(interesting.gene.reads,
         gene.entropy,
         gene.tissue.specificity,
         file=args[length(args)])
} else {
    interesting.isoform.reads <- interesting.reads
    isoform.entropy <- grouping.entropy
    isoform.tissue.specificity <- tissue.specificity
    save(interesting.isoform.reads,
         isoform.entropy,
         isoform.tissue.specificity,
         file=args[length(args)])
}

     
