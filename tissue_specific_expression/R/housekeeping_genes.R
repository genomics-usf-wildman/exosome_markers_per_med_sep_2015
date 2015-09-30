library("data.table")
library("reshape2")
library("yaml")

args <- c("categorized_samples",
          "combined_read_counts",
          "housekeeping_markers")

args <- commandArgs(trailingOnly=TRUE)

### we cannot use the interesting genes subset, because we are
### generating new tissue super-sets here.
load(args[1])
load(args[2])

setkey(categorized.samples,"SRX")

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


### identify house keeping genes. These are genes which are expressed
### in pretty much every tissue. [Note that we'll have to be careful
### with some MT genes as they are sometimes excluded.]

### We're using 60% here primarily to make sure that we get MT genes
### which in some cases have been excluded from alignments.

### genes expressed with FPKM >= 50 in >= 60% of tissues 
housekeeping.genes <-
    copy(gene.reads.wide)#[apply(gene.reads.wide[,-(1:2),with=FALSE]>=50,1,mean) >= 0.60,
                   # ]
housekeeping.genes[,mean:=apply(housekeeping.genes[,-(1:2),with=FALSE],1,mean)]
housekeeping.genes[,percentage:=apply(housekeeping.genes[,-(1:2),with=FALSE] >= 50,1,mean)]
housekeeping.genes[,percentage.10:=apply(housekeeping.genes[,-(1:2),with=FALSE] >= 10,1,mean)]

housekeeping.genes <-
    housekeeping.genes[order(-mean),list(gene_short_name,tracking_id,mean,percentage,percentage.10)]

save(housekeeping.genes,
     file=args[length(args)])

