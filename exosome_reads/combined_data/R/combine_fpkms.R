library(reshape2)
library(data.table)
### set up args for debugging purposes

args <- commandArgs(trailingOnly=TRUE)

### this is the sample_info_from_name.R file
source(args[1])
args <- args[-1]

results.file <- args[length(args)]
args <- args[-length(args)]

### these are the gene FPKM files
gene.files <- grep("_genes.fpkm_tracking",args,value=TRUE)

pb <- txtProgressBar(min=1,max=length(gene.files),style=3)
i <- 0

gene.counts <- list()
for (file in gene.files) {
    sample_info <- sample_info_from_name(file)
    temp <-
        fread(file)
    temp <-
        temp[FPKM > 0,]
    temp[,sample:=sample_info[["name"]]]
    temp[,reads:=sample_info[["reads"]]]
    temp[,subsample:=sample_info[["subsample"]]]
    
    gene.counts[[file]] <-
        temp[,list(gene_id=gene_id[1],
                   gene_short_name=gene_short_name[1],
                   locus=locus[1],
                   FPKM=sum(FPKM),
                   FPKM_conf_lo=sum(FPKM_conf_lo),
                   FPKM_conf_hi=sum(FPKM_conf_hi),
                   FPKM_status=FPKM_status[1]),
             by=list(tracking_id,sample,reads,subsample)]
    i <- i + 1
    setTxtProgressBar(pb,i)
}

close(pb)

gene.counts <- rbindlist(gene.counts,fill=TRUE)

save(file=results.file,
     gene.counts)
