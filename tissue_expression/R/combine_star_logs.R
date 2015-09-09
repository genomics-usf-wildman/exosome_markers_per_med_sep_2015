library(reshape2)
library(data.table)
### set up args for debugging purposes
args <- c("SRX007167/SRX007167_isoforms.fpkm_tracking",
          "SRX007167/SRX007167_genes.fpkm_tracking",
          "SRX007167/SRX007167_star/Log.final.out")
args <- c(gsub("SRX007167","SRX007165",args),args)
args <- c(args,"cobmined_read_counts")

args <- commandArgs(trailingOnly=TRUE)

results.file <- args[length(args)]
args <- args[-length(args)]

### these files contain the alignment results from STAR
star.log.files <- grep("_star/Log.final.out",args,value=TRUE)

.get.srx <- function(x){
    gsub("(SRX\\d+)/.+","\\1",x)
}

pb <- txtProgressBar(min=1,max=length(star.log.files),style=3)
i <- 0

star.logs <- list()
for (file in star.log.files) {
    srx.accession <- .get.srx(file)
    star.log <- read.table(file,sep="|",fill=TRUE,stringsAsFactors=FALSE)
    colnames(star.log) <- c("field","value")
    star.log$value <- gsub("\\t","",star.log$value)
    star.log$field <- gsub("(^\\s+|\\s+$)","",star.log$field)
    star.log$srx <- srx.accession
    star.logs[[srx.accession]] <- data.table(star.log)[!grepl(":",field),]
    i <- i + 1
    setTxtProgressBar(pb,i)
}
close(pb)

star.logs <- rbindlist(star.logs)

save(file=results.file,
     star.logs)
