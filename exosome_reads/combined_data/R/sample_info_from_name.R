##' Given a filename, extract information about the sample
##'
##' Given a filename including path, like
##' E0573_TCTCGCGC-TCAGAGCC_L004_R1_001_genes.fpkm_tracking, return
##' information about the sample
##' @title sample_info_from_name
##' @param filename -- filename of the sample
##' @return list of sample information including sample name, read,
##'     and type
##' @author Don Armstrong
sample_info_from_name <- function(filename) {
    fn.components <-
        stringr::str_match(filename,
                           paste0("(?:^|.*\\/)",
                                  ## sample id
                                  "((?:SRX\\d+)?)_split",
                                  ## reads
                                  "_r(\\d+)",
                                  ## subsample
                                  "_s(\\d+)",
                                  ## file type
                                  "_(fastqc\\.zip|(isoforms|genes)\\.fpkm_tracking|",
                                  "star/Log.final.out|star/ReadsPerGene.out.tab|",
                                  "kallisto.txt)"
                                  ))[1,]
    sample_info <- list()
    sample_info[["name"]] <-
        fn.components[2]
    sample_info[["reads"]] <-
        fn.components[3]
    sample_info[["subsample"]] <-
        fn.components[4]

    ## everything is trimmed unless this is a fastqc and this
    ## fn.components[6] is ""
    sample_info[["trimmed"]] <-
        TRUE
    sample_info[["filetype"]] <-
        fn.components[5]
    sample_info
}
