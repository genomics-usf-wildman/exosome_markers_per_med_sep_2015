library(data.table)
args <- c(grep("GSE\\d+_family.txt$",dir(),value=TRUE),
          "gse_families.txt")

args <- commandArgs(trailingOnly=TRUE)

gse.files <- args[-length(args)]
output.file <- args[length(args)]

gse.families <- list()
for (gse.file in gse.files) {
    gse.num <- gsub("(GSE\\d+).+","\\1",gse.file)
    gse.families[[gse.file]] <-
        fread(gse.file)
    gse.families[[gse.file]][,gse:=gse.num]
    
}

gse.samples <- rbindlist(gse.families)

gse.samples[,sra:=gsub("^.*(SR(X|A)\\d+)\\s*$","\\1",sra)]
save(gse.samples,file=output.file)
