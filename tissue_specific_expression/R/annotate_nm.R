library("biomaRt")
library("data.table")
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

args <- c("housekeeping_genes_eisenberg_2003.txt",
          "housekeeping_genes_eisenberg_2003_annotated.txt")

args <- commandArgs(trailingOnly=TRUE)

housekeeping.genes <- fread(args[1],header=FALSE)
if (ncol(housekeeping.genes) > 1) {
    setnames(housekeeping.genes,c("symbol","refseq"))
} else {
    setnames(housekeeping.genes,"refseq")
}
   

housekeeping.genes.annotated <-
    data.table(getBM(attributes=c("refseq_mrna", "ensembl_gene_id", "hgnc_symbol"),
                     filters = "refseq_mrna", values = housekeeping.genes[,refseq],
                     mart= ensembl))

write.table(housekeeping.genes.annotated,
            sep="\t",
            quote=FALSE,
            row.names=FALSE,
            file=args[length(args)])
