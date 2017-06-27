library(reshape2)
library(data.table)
### set up args for debugging purposes

args <- commandArgs(trailingOnly=TRUE)

### this is the sample_info_from_name.R file
source(args[1])
args <- args[-1]

results.file <- args[length(args)]
args <- args[-length(args)]

### these are the isoform fpkm files
isoform.files <- grep("_isoforms.fpkm_tracking",args,value=TRUE)

pb <- txtProgressBar(min=1,max=length(isoform.files),style=3)
i <- 0

isoform.counts <- list()
for (file in isoform.files) {
    sample_info <- sample_info_from_name(file)
    temp <-
        fread(file)
    temp <-
        temp[FPKM > 0,]
    temp[,sample:=sample_info[["name"]]]
    temp[,reads:=sample_info[["reads"]]]
    temp[,subsample:=sample_info[["subsample"]]]
    isoform.counts[[file]] <- temp
    i <- i + 1
    setTxtProgressBar(pb,i)
}

isoform.counts <- rbindlist(isoform.counts)
isoform.counts[,class_code:=NULL]
isoform.counts[,nearest_ref_id:=NULL]
isoform.counts[,tss_id:=NULL]
isoform.counts[,length:=NULL]
isoform.counts[,coverage:=NULL]

save(file=results.file,
     isoform.counts)
