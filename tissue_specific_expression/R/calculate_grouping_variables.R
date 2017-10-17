library("data.table")
library("reshape2")
library("entropy")

args <- c("categorized_samples",
          "combined_read_counts",
          "interesting_gene_reads"
          )

args <- commandArgs(trailingOnly=TRUE)


## categorized_samples
load(args[1])
## combined_read_counts
load(args[2])

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
    rm(gene.counts)
} else {
    read.counts <- isoform.counts
    rm(isoform.counts)
}

setnames(read.counts,"srx","SRX")
setkey(read.counts,"SRX")

### combine sample name and srx into the gene counts file
combined.reads <-
    categorized.samples[,list(Sample_Group,SRX)
                        ][read.counts[,list(tracking_id,gene_id,
                                            gene_short_name,
                                            FPKM,
                                            FPKM_status,
                                            SRX)]]
rm(read.counts)
rm(categorized.samples)
gc()

.narm.mean <- function(x){
    if(length(x)==0) {return(0)} else {return(mean(x,na.rm=TRUE))}
}

calculate.entropy <- function(x,bins=10,method="MM"){
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    result <- 0
    try(result <- entropy(discretize(log(x+1),numBins=bins),method=method),
        silent=TRUE)
    return(result)
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

which.max.both <- function(expression,group,names=list("tissue_max"=1,"tissue_max_2nd"=2)) {
    temp <- group[order(expression,decreasing=TRUE)]
    lapply(names,function(x){temp[x]})
}

combined.reads.mean <-
    combined.reads[,list(gene_id=gene_id[1],
                         gene_short_name=gene_short_name[1],
                         mean_FPKM=.narm.mean(FPKM)),
                   by=list(tracking_id,Sample_Group)]
rm(combined.reads)

grouping.variables <-
    combined.reads.mean[Sample_Group!="universal human reference",
                        c(list(gene_id=gene_id[1],
                               gene_short_name=gene_short_name[1],
                               entropy=calculate.entropy(mean_FPKM),
                               entropy.log2=calculate.entropy(log2(mean_FPKM+1)),
                               max=max(mean_FPKM,na.rm=TRUE),
                               max_2nd=-sort.int(-mean_FPKM,partial=2)[2],
                               var=var(mean_FPKM,na.rm=TRUE),
                               var.log2=var(log2(mean_FPKM+1),na.rm=TRUE),
                               mean=mean(mean_FPKM,na.rm=TRUE),
                               mean.log2=mean(log2(mean_FPKM+1),na.rm=TRUE)),
                          which.max.both(mean_FPKM,
                                         Sample_Group,
                                         names=list(tissue_max=1,
                                                    tissue_max_2nd=2)),
                          list(tissue.specificity.index=tissue.specificity.index(mean_FPKM),
                               tissue.specificity.index.log2=tissue.specificity.index(log2(mean_FPKM+1))
                               )),
                        by=list(tracking_id)]

if (is.genes) {
    gene.grouping.variables <- grouping.variables
    save(gene.grouping.variables,
         file=args[length(args)])
} else {
    isoform.grouping.variables <- grouping.variables
    save(isoform.grouping.variables,
         file=args[length(args)])
}

     
