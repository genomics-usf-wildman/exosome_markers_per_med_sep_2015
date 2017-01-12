library("data.table")
library("reshape2")
library("entropy")

args <- c("categorized_samples",
          "combined_read_counts",
          "sample_reads_wide"
          )


load(args[1])
load(args[2])

setkey(categorized.samples,"SRX")
read.counts <- gene.counts

grouping.by.sample.formula <-
    as.formula(tracking_id~Sample_Group)
grouping.formula <-
    as.formula(tracking_id~.)
setnames(read.counts,"srx","SRX")
setkey(read.counts,"SRX")
combined.reads <-
    categorized.samples[,list(Sample_Group,SRX)][read.counts]
.narm.mean <- function(x){
    if(length(x)==0) {return(0)} else {return(mean(x,na.rm=TRUE))}
}
c.reads.wide <-
    data.table(dcast(combined.reads,
                     grouping.by.sample.formula,
                     fun.aggregate=.narm.mean,
                     value.var="FPKM"
                     ))

save(c.reads.wide,
     file=args[length(args)])
