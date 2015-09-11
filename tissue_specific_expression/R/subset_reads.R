library(data.table)
library(reshape2)
library(entropy)

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
grouping.entropy <- data.table(grouping.entropy)
setkey(grouping.entropy,"tracking_id")
    

grouping.max <-
    dcast(combined.reads,
          grouping.formula,
          fun.aggregate=function(x){max(x,na.rm=TRUE)},
          value.var="FPKM")
colnames(grouping.max)[2] <- "max"
grouping.max <- data.table(grouping.max)
setkey(grouping.max,"tracking_id")

grouping.var <-
    dcast(combined.reads,
          grouping.formula,
          fun.aggregate=function(x){var(x,na.rm=TRUE)},
          value.var="FPKM")
colnames(grouping.var)[2] <- "var"
grouping.var <- data.table(grouping.var)
setkey(grouping.var,"tracking_id")

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

min.entropy <- 2
min.tissue.specificity <- 0.985
min.expression <- 10
min.var <- 1

tissue.specificity <-
    data.table(data.frame(tracking_id=c.reads.wide[,1],
                          tissue.specificity.index=apply(c.reads.wide[,-1],1,
                              tissue.specificity.index)))
setkey(tissue.specificity,"tracking_id")

grouping.variables <-
    grouping.entropy[grouping.var][grouping.max][tissue.specificity]

interesting.groups <-
    grouping.variables[((entropy >= min.entropy |
                          tissue.specificity.index >= min.tissue.specificity)
                        & max >= min.expression
                        & var >= min.var) |
                           tracking_id %in% additional.tracking.ids,
                       tracking_id]
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
    gene.grouping.variables <- grouping.variables
    save(interesting.gene.reads,
         gene.grouping.variables,
         file=args[length(args)])
} else {
    interesting.isoform.reads <- interesting.reads
    isoform.grouping.variables <- grouping.variables
    save(interesting.isoform.reads,
         isoform.grouping.variables,
         file=args[length(args)])
}

     
