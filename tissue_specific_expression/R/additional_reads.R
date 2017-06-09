library("data.table")
library("reshape2")

args <- c("categorized_samples",
          "interesting_gene_reads",
          "additional_interesting_genes.txt",
          "all_interesting_gene_reads"
          )

args <- commandArgs(trailingOnly=TRUE)


## categorized_samples
load(args[1])
## combined_read_counts
load(args[2])
## interesting_gene_reads
load(args[3])
## additional_interesting_genes.txt
additional.genes <- fread(args[4])
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
    interesting.reads <- interesting.gene.reads
    grouping.variables <- gene.grouping.variables
    rm(gene.counts)
} else {
    read.counts <- isoform.counts
    interesting.reads <- interesting.isoform.reads
    grouping.variables <- isoform.grouping.variables
    rm(isoform.counts)
}

additional.tracking.ids <-
    read.counts[!duplicated(tracking_id),
                list(tracking_id,gene_short_name)][gene_short_name %in%
                                                   additional.genes[[1]]][[1]]

additional.tracking.ids <-
    additional.tracking.ids[!(additional.tracking.ids %in%
                              interesting.reads[,unique(tracking_id)])]
setkey(read.counts,"tracking_id")
read.counts <-
    read.counts[additional.tracking.ids]

### combine sample name and srx into the gene counts file
combined.reads <-
    categorized.samples[,list(Sample_Group,SRX)
                        ][read.counts[,list(tracking_id,gene_id,
                                            gene_short_name,
                                            FPKM,
                                            SRX)]]

interesting.reads <-
    rbindlist(list(interesting.reads,
                   combined.reads[,
                                  list(Sample_Group,
                                       SRX,
                                       tracking_id,
                                       gene_id,
                                       gene_short_name,
                                       FPKM)]))

if (is.genes) {
    interesting.gene.reads <- interesting.reads
    save(interesting.gene.reads,
         file=args[length(args)])
} else {
    interesting.isoform.reads <- interesting.reads
    save(interesting.isoform.reads,
         file=args[length(args)])
}

     
