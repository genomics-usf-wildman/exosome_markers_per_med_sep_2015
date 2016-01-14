library("data.table")
library("reshape2")
library("yaml")
library("parallel")

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

gene.reads.wide[,name_or_tracking:=ifelse(is.na(gene_short_name),tracking_id,gene_short_name)]
setkey(gene.reads.wide,name_or_tracking)

housekeeping.genes <- list()

housekeeping.gene.levels <- function(gene) {
    gene.expression <-
        apply(gene.reads.wide[gene,-(c(1:2,ncol(gene.reads.wide))),
                              with=FALSE],
              2,sum)
    return(data.table(
        gene=gene,
        tracking_id=gene.reads.wide[gene,tracking_id][1],
        mean=mean(gene.expression),
        sd=sd(gene.expression),
        percentage=mean(gene.expression>50),
        log.sd=sd(log(gene.expression[gene.expression>0])),
        log.mean=mean(log(gene.expression[gene.expression>0])),
        percentage.10=mean(gene.expression>10),
        sd.10=sd(gene.expression[gene.expression>10]),
        percentage.1=mean(gene.expression>1),
        sd.1=sd(gene.expression[gene.expression>1])
        ))
}

housekeeping.genes <-
    mclapply(unique(gene.reads.wide[,name_or_tracking]),
             housekeeping.gene.levels,
             mc.cores=12)

housekeeping.genes <-
    rbindlist(housekeeping.genes)

invisible(housekeeping.genes[,scaled.sd:=log.sd/log.mean])

### identify house keeping genes. These are genes which are expressed
### in pretty much every tissue. [Note that we'll have to be careful
### with some MT genes as they are sometimes excluded.]

### We're using 60% here primarily to make sure that we get MT genes
### which in some cases have been excluded from alignments.

### genes expressed with FPKM >= 50 in >= 80% of tissues and with a scaled SD of less than 0.12

housekeeping.genes <-
    housekeeping.genes[(percentage >= 0.8 & scaled.sd < 0.12) |
                       (grepl("^MT-",gene) & mean > 100)
                      ,]
setnames(housekeeping.genes,"gene","gene_short_name")

save(housekeeping.genes,
     file=args[length(args)])

