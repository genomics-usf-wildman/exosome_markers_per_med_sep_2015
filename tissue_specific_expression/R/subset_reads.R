library("data.table")
library("reshape2")
library("entropy")

args <- c("categorized_samples",
          "combined_read_counts",
          "additional_interesting_genes.txt",
          "gene_grouping_variables",
          "interesting_gene_reads"
          )

args <- commandArgs(trailingOnly=TRUE)


## categorized_samples
load(args[1])
## combined_read_counts
load(args[2])
## additional_interesting_genes.txt
additional.genes <- fread(args[3])
## gene_grouping_variables
load(args[4])


is.genes <- TRUE
if (any(grepl("isoform",args[length(args)]))) {
    is.genes <- FALSE
}

setkey(categorized.samples,"SRX")


## we're using tracking_id for genes instead of gene_short_name
## because SNORD5 and some other genes have multiple gene ids but also
## share a common name
if (is.genes) {
    read.counts <- gene.counts
    grouping.variables <- gene.grouping.variables
    rm(gene.counts)
} else {
    read.counts <- isoform.counts
    grouping.variables <- isoform.grouping.variables
    rm(isoform.counts)
}

additional.tracking.ids <-
    read.counts[!duplicated(tracking_id),
                list(tracking_id,gene_short_name)][gene_short_name %in%
                                                       additional.genes[[1]]][[1]]

setnames(read.counts,"srx","SRX")
setkey(read.counts,"SRX")

min.entropy <- 2
min.tissue.specificity <- 0.985
min.expression <- 10

### combine sample name and srx into the gene counts file
combined.reads <-
    categorized.samples[,list(Sample_Group,SRX)
                        ][read.counts[,list(tracking_id,gene_id,
                                            gene_short_name,
                                            FPKM,
                                            SRX)]]
interesting.groups <-
    grouping.variables[((entropy >= min.entropy |
                         tissue.specificity.index >= min.tissue.specificity) &
                        max >= min.expression) |
                       tracking_id %in% additional.tracking.ids,
                       tracking_id]

setkey(combined.reads,"tracking_id") 

interesting.reads <-
    combined.reads[interesting.groups,
                   list(Sample_Group,
                        SRX,
                        tracking_id,
                        gene_id,
                        gene_short_name,
                        FPKM,
                        FPKM_status)]
if (is.genes) {
    interesting.gene.reads <- interesting.reads
    save(interesting.gene.reads,
         file=args[length(args)])
} else {
    interesting.isoform.reads <- interesting.reads
    save(interesting.isoform.reads,
         file=args[length(args)])
}

     
