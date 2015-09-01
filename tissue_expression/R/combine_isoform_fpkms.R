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

.get.srx <- function(x){
    gsub("(SRX\\d+)/.+","\\1",x)
}

pb <- txtProgressBar(min=1,max=length(isoform.files),style=3)
i <- 0

isoform.counts <- list()
for (file in isoform.files) {
    srx.accession <- .get.srx(file)
    isoform.counts[[srx.accession]] <-
        fread(file)
    isoform.counts[[srx.accession]][,srx:=srx.accession]
    i <- i + 1
    setTxtProgressBar(pb,i)
}

isoform.counts <- rbindlist(isoform.counts)

save(file=results.file,
     isoform.counts)
