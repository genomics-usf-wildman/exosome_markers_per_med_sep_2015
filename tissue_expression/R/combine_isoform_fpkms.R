library(reshape2)
library(data.table)
### set up args for debugging purposes
args <- c("SRX007167/SRX007167_isoforms.fpkm_tracking"
          )
args <- c(gsub("SRX007167","SRX007165",args),args)
args <- c(args,"combined_read_counts")

args <- commandArgs(trailingOnly=TRUE)

results.file <- args[length(args)]
args <- args[-length(args)]

### these are the isoform fpkm files
isoform.files <- grep("_isoforms.fpkm_tracking",args,value=TRUE)

gtex.file <- grep("GTEx_Analysis",args,value=TRUE)

.get.srx <- function(x){
    gsub("(SRX\\d+)/.+","\\1",x)
}

pb <- txtProgressBar(min=1,max=length(isoform.files),style=3)
i <- 0

isoform.counts <- list()
for (file in isoform.files) {
    srx.accession <- .get.srx(file)
    isoform.counts[[srx.accession]] <-
        fread(file)[,list(tracking_id,gene_id,gene_short_name,FPKM,FPKM_status)]
    isoform.counts[[srx.accession]][,srx:=srx.accession]
    i <- i + 1
    setTxtProgressBar(pb,i)
}
gc()

gtex.counts <- fread(paste0("zcat ",gtex.file),skip="TargetID")

gtex.counts[,TargetID:=gsub("\\.\\d+","",TargetID)]
gtex.counts[,Gene_Symbol:=gsub("\\.\\d+","",Gene_Symbol)]
gtex.counts[,Chr:=NULL]
gtex.counts[,Coord:=NULL]
setnames(gtex.counts,c("TargetID","Gene_Symbol"),c("tracking_id","gene_id"))
gtex.counts <-
    melt(gtex.counts,
         id.vars=c("tracking_id","gene_id"),
         value.name="FPKM",
         variable.name="srx")
gc()
isoform.counts[["gtex"]] <- gtex.counts

isoform.counts <- rbindlist(isoform.counts,fill=TRUE)

save(file=results.file,
     isoform.counts)
