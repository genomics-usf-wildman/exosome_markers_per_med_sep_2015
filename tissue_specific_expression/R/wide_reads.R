library(data.table)
library(reshape2)

args <- c("interesting_gene_reads",
          "interesting_gene_reads_wide"
          )

args <- commandArgs(trailingOnly=TRUE)

load(args[1])

is.genes <- TRUE
if (any(grepl("isoform",args[length(args)]))) {
    is.genes <- FALSE
}

if(is.genes) {
    combined.reads <- interesting.gene.reads
} else {
    combined.reads <- interesting.isoform.reads
}
    
grouping.by.sample.formula <-
    as.formula(gene_short_name+tracking_id~Sample_Group)
grouping.formula <-
    as.formula(tracking_id~.)

.narm.mean <- function(x){
    if(length(x)==0) {return(0)} else {return(mean(x,na.rm=TRUE))}
}

c.reads.wide <-
    dcast(combined.reads,
          grouping.by.sample.formula,
          fun.aggregate=.narm.mean,
          value.var="FPKM"
          )

if (is.genes) {
    interesting.gene.reads.wide <- data.table(c.reads.wide)
    save(interesting.gene.reads.wide,
         file=args[length(args)])
} else {
    interesting.isoform.reads.wide <- data.table(c.reads.wide)
    save(interesting.isoform.reads.wide,
         file=args[length(args)])
}
