args <- c("gse_family_xml/gse_families","gse_expression_samples")

args <- commandArgs(trailingOnly=TRUE)

load(args[1])

## ignore bisulfite stuff and other irrelevant sampls   
gse.samples <- gse.samples[!(grepl("Bisulfite-Seq|bisulf\\. sequencing|ChIP-[sS]eq|TAB-Seq|Chromatin accessibility|Exon arrays?|Genotyping|Whole genome sequencing|MethylC-seq",title)
                             |grepl("[bB]isulfite [sS]equencing",title)
                             |grepl("Digital Genomic Footprinting",title)
                             |grepl("Mouse",title)
                             |grepl("Ischemic|Idiopathic",title)
                             |grepl("^(WGBS|RRBS)",title)
                             ),]
gse.samples[,tissue:=as.character(NA)]
gse.samples[,type:=as.character(NA)]
gse.samples[,replicate:=as.character(NA)]

## GSE16256
gse.samples[gse=="GSE16256"
            & grepl(".+from the (.+?) cell line;.*",title),
            tissue:=gsub(".+from the (.+?) cell line;.*","\\1",title)]
gse.samples[gse=="GSE16256"
            & grepl(".+from the (.+?) cell line clone ([^\\;]+)\\;",title),
            tissue:=gsub(".+from the (.+?) cell line clone ([^\\;]+)\\;.*","\\1 clone \\2",title)]
gse.samples[gse=="GSE16256"
            & grepl(".+(?:from the|of) (.+?) (?:cell line differentiated into|[dD]erived) (.+?);.*",title),
            tissue:=gsub(".+(?:from the|of) (.+?) (?:cell line differentiated into|[dD]erived) (.+?);.*",
                "\\1 derived \\2",title)]
gse.samples[gse=="GSE16256"
            & grepl(".+of (?:STL\\d+ )?(.+?) [cC]ells;.*",title),
            tissue:=gsub(".+of (?:STL\\d+ )?(.+?) [cC]ells;.*","\\1",title)]

## GSE 18927
gse.samples[gse=="GSE18927"
            & grepl(".+\\(Total RNA\\) .+?\\;.+",title),
            tissue:=gsub(".+\\(Total RNA\\) (.+?);.+","\\1",title)]

## GSE 64016 H1 single cell
gse.samples <- gse.samples[!(gse=="GSE64016" &
                                 grepl("^(G1|S|G2)_Exp\\d.\\d+$",title)),]
gse.samples[gse=="GSE64016"
            & grepl("^H1_Exp\\d+.\\d+$",title),
            tissue:="H1"]
gse.samples[gse=="GSE64016"
            & grepl("^H1_Exp\\d+.\\d+$",title),
            type:="single cell"]
gse.samples[gse=="GSE64016"
            & grepl("^H1_Exp\\d+.\\d+$",title),
            replicate:=paste0(sra,
                gsub("^H1_Exp(\\d+.\\d+)$","\\1",title))]
                    
                
### GSE16368
gse.samples <- gse.samples[!(gse=="GSE16368" &
                                 grepl("Histone|MRE[-_]Seq|MeDIP|Input|H3K",title)),]
gse.samples[gse=="GSE16368",
            tissue:=gsub("(?:.+analysis of |(?:sm|)RNA-[sS]eq (?:\\(?polyA\\+\\)? )?)([^,;]+?)(?:| from [^,;]+)[;,].+","\\1",title)]
gse.samples[gse=="GSE16368",
            type:=gsub("^(.+) analysis of ([^,;]+?)(?:| from \\S+?)[;,].+","\\1",title)]
gse.samples[gse=="GSE16368"
            & grepl(".+analysis of ([^,;]+?) from ([^,;]+)[;,].+",title),
            replicate:=gsub(".+analysis of ([^,;]+?) from ([^,;]+)[;,].+","\\2",title)]

### GSE63818
### http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63818
gse.samples <-
    gse.samples[!(gse=="GSE63818" &
                      grepl("PBAT",title)),]

gse.samples[gse=="GSE63818"
            & grepl("^([MF])_(PGC|Soma)_(\\d+)W_(embryo\\d+_(?:sc|ps)\\d+)$",title),
            tissue:=gsub("^([MF])_(PGC|Soma)_(\\d+)W_(embryo\\d+_(?:sc|ps)\\d+)$","\\2 \\1 \\3W",title)]

gse.samples[gse=="GSE63818"
            & grepl("^(.+)_(\\d+)W_(embryo\\d+_HiSeq2[50]00(?:_deeper)?)$",title),
            tissue:=gsub("^(.+)_(\\d+)W_(embryo\\d+_HiSeq2[50]00(?:_deeper)?)$","\\1 \\2W",title)]
    
gse.samples[gse=="GSE63818"
            & grepl("^(PGC|Soma)_(\\d+)W_embryo\\d+_([MF])_(?:rep\\d+_)?HiSeq2[50]00(?:_deeper)?$",
                    title),
            tissue:=gsub("^(PGC|Soma)_(\\d+)W_embryo\\d+_([MF])_(?:rep\\d+_)?HiSeq2[50]00(?:_deeper)?$",
                "\\1 \\2 \\3W",title)]
    

### GSE57945
### http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57945
gse.samples <-
    gse.samples[!(gse=="GSE57945" &
                      !grepl("Not IBD",title)),]
gse.samples[gse=="GSE57945",
            tissue:="ileal"]

### GSE57345
### http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57345
gse.samples[gse=="GSE57345"
           & grepl("^[PC]\\d+ Non-failing ",title),
            tissue:=gsub("^[PC]\\d+ Non-failing ","",title)]

gse.samples[gse=="GSE60012",
            gsub
