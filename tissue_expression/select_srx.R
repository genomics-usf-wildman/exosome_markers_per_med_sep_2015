library(data.table)
args <- c("roadmap_epigenomics.tsv","srx_to_sra.txt")

args <- commandArgs(trailingOnly=TRUE)


##' Get SRR number back from the database or the srx_to_sra.txt file
##'
##' .. content for \details{} ..
##' @title get.srrs
##' @param srx SRX entry
##' @param srr.listing srr.listing data.table (load from srx.to.sra)
##' @param db postgresql sradb database name
##' @return 
##' @author Don Armstrong
get.srrs <- function(srx,srr.listing=NULL,db="sradb") {
    srrs <- as.character(NULL)
    for (i in srx) {
        if (!is.null(srr.listing) &&
            any(i==srr.listing[,SRX])) {
            srrs[i] <- paste(sep="",collapse=",",srr.listing[SRX==i,SRR])
        } else {
            library("RPostgreSQL")
            con <- dbConnect(PostgreSQL(),dbname="sradb")
            srrs.result <-
                dbGetQuery(con,paste0("SELECT run_accession ",
                                      "FROM sra WHERE experiment_accession='",
                                      i,
                                      "'"))
            srrs[i] <- paste(sep="",collapse=",",srrs.result[,"run_accession"])
            dbDisconnect(con)
        }
    }
    return(srrs)
}

##' Get the number of reads for the SRR in this SRX
##'
##' Calls vdb-dump to determine the number of reads in the first SRR
##' of this SRX unless we've previously checked on this SRX.
##' @title 
##' @param srx 
##' @param srr.listing 
##' @param db 
##' @return 
##' @author Don Armstrong
get.nreads <- function(srx,srr.listing=NULL,db="sradb") {
    .nreads <- function(srr) {
        gsub(".+\\: ","",system(paste0("vdb-dump -R 1 -C NREADS ",srr),intern=TRUE))
    }
    nreads <- as.numeric(NULL)
    for (i in srx) {
        if (!is.null(srr.listing) &&
            any(i==srr.listing[,SRX])
            ) {
            srr <- srr.listing[SRX==i,SRR][1]
            if (is.na(srr.listing[,NREADS])) {
                nreads[i] <-
                    .nreads(srr)
            } else {
                nreads[i] <- srr.listing[SRX==i,NREADS][1]
            }
        } else {
            library("RPostgreSQL")
            con <- dbConnect(PostgreSQL(),dbname="sradb")
            srrs.result <-
                dbGetQuery(con,paste0("SELECT run_accession ",
                                      "FROM sra WHERE experiment_accession='",
                                      i,
                                      "'"))
            srr <- srrs.result[,"run_accession"][1]
            nreads[i] <-
                .nreads(srr)
            dbDisconnect(con)
        }
    }
    return(nreads)
}

roadmap.ep <- fread(args[1])
setnames(roadmap.ep,names(roadmap.ep),
         gsub("^_|_$","",
              gsub("[ _-]+","_",
                   gsub("[#]","",names(roadmap.ep)))))

srx.to.sra <- fread(args[2],sep="\t")

chosen.samples <-
    roadmap.ep[SRA_FTP!="" & Experiment=="mRNA-Seq",]
chosen.samples[,SRX:=gsub(".+/(SRX\\d+)$","\\1",SRA_FTP)]
chosen.samples[,SRR:=get.srrs(SRX,srx.to.sra)]
chosen.samples[,NREADS:=get.nreads(SRX,srx.to.sra)]

### only write this out if we updated the data; someone might not have
### a database who wants to regenerate this file
if (is.null(srx.to.sra) ||
    any(!(chosen.samples[,SRX] %in% srx.to.sra[,SRX]))) {
    write.table(chosen.samples[,list(SRX,SRR,NREADS)],
                file="srx_to_sra.txt",sep="\t",row.names=FALSE)
}

save(chosen.samples,file=args[length(args)])
