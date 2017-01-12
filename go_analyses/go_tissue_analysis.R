library("data.table")
library(doMC)
ncore = multicore:::detectCores()
registerDoMC(cores = ncore)
library("topGO")
library("org.Hs.eg.db")

args <- c("sample_reads_wide","tissue_specific_markers","go_tissue_analyses")

args <- commandArgs(trailingOnly=TRUE)

load(args[1])
load(args[2])



do.go.analysis <- function(tissue,go.ontology) {

    all.genes <- c.reads.wide[eval(parse(text=paste0("`",tissue,"`>1"))),list(tracking_id)]
    setkey(all.genes,"tracking_id")
    all.genes[,tissue.specific:=0]
    all.genes[as.character(tissue.specific.genes[tissue==tissue ,tracking_id]),
              tissue.specific:=1]
    all.genes[,egid:=as.vector(unlist(sapply(mget(all.genes[,tracking_id],
                   org.Hs.egENSEMBL2EG,ifnotfound=NA),function(x){x[1]})))]
    all.genes <- all.genes[!is.na(egid),]

    gene.list <- factor(all.genes[,tissue.specific])
    names(gene.list) <- all.genes[,egid]
    
    godata <- new("topGOdata",
                  ontology=go.ontology,
                  allGenes=gene.list,
                  annot=annFUN.org,
                  mapping="org.Hs.eg.db")
    test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
    result.fisher <- getSigGroups(godata,test.stat)
    res.table <- 
        data.table(GenTable(godata,
                            classic=result.fisher,
                            orderBy="classic",
                            topNodes=length(result.fisher@score)))
    res.table[,classic:=as.numeric(gsub("< ","",classic))]
    setnames(res.table,"classic","p value")
    res.table[,FDR:=p.adjust(method="BH",`p value`)]
    res.table[,Ontology:=go.ontology]
    return(res.table)
}
 
results <- 
    c(lapply("placenta",#tissue.specific.genes[!is.na(tissue),unique(tissue)],
             do.go.analysis,
             go.ontology="MF"),
      lapply("placenta",#tissue.specific.genes[!is.na(tissue),unique(tissue)],
             do.go.analysis,
             go.ontology="BP")
      )

print(results)

go.results <- rbindlist(results)

go.results <- go.results[`p value` <= 0.05,]

save(file=args[length(args)],
     go.results)
               
