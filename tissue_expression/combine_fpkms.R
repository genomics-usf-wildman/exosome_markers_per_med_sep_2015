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

if (!all.equal(length(isoform.files),
               length(gene.files),
               length(star.log.files))) {
    stop("the number of isoform files, gene files, and star log files must be equal")
}

.get.srx <- function(x){
    gsub("(SRX\\d+)/.+","\\1",x)
}

if (!all.equal(sort(.get.srx(isoform.files)),
               sort(.get.srx(gene.files)),
               sort(.get.srx(star.log.files)))) {
    stop("The SRX accessions of the gene counts, isoform counts, and star log files are not equal")
}


pb <- txtProgressBar(min=1,max=length(isoform.files)*3,style=3)
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

gene.counts <- list()
for (file in gene.files) {
    srx.accession <- .get.srx(file)
    gene.counts[[srx.accession]] <-
        fread(file)
    gene.counts[[srx.accession]][,srx:=srx.accession]
    i <- i + 1
    setTxtProgressBar(pb,i)
}

star.logs <- list()
for (file in star.log.files) {
    srx.accession <- .get.srx(file)
    star.log <- read.table(file,sep="|",fill=TRUE,stringsAsFactors=FALSE)
    colnames(star.log) <- c("field","value")
    star.log$value <- gsub("\\t","",star.log$value)
    star.log$field <- gsub("^\\s+","",star.log$field)
    star.log$srx <- srx.accession
    star.logs[[srx.accession]] <- data.table(star.log)[!grepl(":",field),]
    i <- i + 1
    setTxtProgressBar(pb,i)
}
close(pb)

gene.counts <- rbindlist(gene.counts)
isoform.counts <- rbindlist(isoform.counts)
star.logs <- rbindlist(star.logs)

save(file=results.file,
     gene.counts,
     isoform.counts,
     star.logs)
