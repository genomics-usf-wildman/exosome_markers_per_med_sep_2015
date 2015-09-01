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

### these are the isoform fpkm files
isoform.files <- grep("_isoforms.fpkm_tracking",args,value=TRUE)
### these are the gene FPKM files
gene.files <- grep("_genes.fpkm_tracking",args,value=TRUE)
### these files contain the alignment results from STAR
star.log.files <- grep("_star/Log.final.out",args,value=TRUE)

gtex.file <- grep("GTEx_Analysis",args,value=TRUE)

.get.srx <- function(x){
    gsub("(SRX\\d+)/.+","\\1",x)
}

pb <- txtProgressBar(min=1,max=length(gene.files),style=3)
i <- 0

gene.counts <- list()
for (file in gene.files) {
    srx.accession <- .get.srx(file)
    gene.counts[[srx.accession]] <-
        fread(file)
    gene.counts[[srx.accession]][,srx:=srx.accession]
    i <- i + 1
    setTxtProgressBar(pb,i)
}

close(pb)

### deal with gtex data
gtex.counts <- fread(paste0("zcat ",gtex.file),skip="Name")

gtex.counts[,Name:=gsub("\\.\\d+","",Name)]
setnames(gtex.counts,c("Name","Description"),c("gene_id","gene_short_name"))
gtex.counts <-
    melt(gtex.counts,
         id.vars=c("gene_id","gene_short_name"),
         value.name="FPKM",
         variable.name="srx")
gene.counts[["gtex"]] <- gtex.counts

gene.counts <- rbindlist(gene.counts,fill=TRUE)

gene.counts[,class_code:=NULL]
gene.counts[,nearest_ref_id:=NULL]
gene.counts[,tss_id:=NULL]
gene.counts[,length:=NULL]
gene.counts[,coverage:=NULL]

save(file=results.file,
     gene.counts)
