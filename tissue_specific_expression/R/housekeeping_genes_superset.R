library("data.table")

args <- c("housekeeping_genes",
          "housekeeping_genes_eisenberg_2013_annotated.txt",
          "housekeeping_genes_superset.txt"
          )

args <- commandArgs(trailingOnly=TRUE)

load(args[1])
housekeeping.genes.eisenberg <-
    fread(args[2])

setnames(housekeeping.genes.eisenberg,
         c("hgnc_symbol","ensembl_gene_id"),
         c("gene_short_name","tracking_id"))
housekeeping.genes.superset <-
    rbindlist(list(housekeeping.genes,housekeeping.genes.eisenberg),
              use.names=TRUE,fill=TRUE)
housekeeping.genes.superset <-
    housekeeping.genes.superset[!duplicated(tracking_id),list(gene_short_name,tracking_id)]

save(housekeeping.genes.superset,file=args[length(args)])
