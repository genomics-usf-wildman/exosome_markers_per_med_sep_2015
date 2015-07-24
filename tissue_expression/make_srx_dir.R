library(tools)
library(data.table)
args <- commandArgs(trailingOnly=TRUE)

load(args[1])

for (srx in chosen.samples[,SRX]) {
    if (!dir.exists(srx)) {
        dir.create(srx)
    }
    wd <- getwd()
    setwd(srx)
    if (file.exists("Makefile"))
        file.remove("Makefile")
    file.symlink("../srx_makefile.mk","Makefile")
    srr.mk <- paste0("SRX=",srx,"\n",
                     "SRRS=",gsub(","," ",chosen.samples[SRX==srx,SRR]),"\n",
                     "NREADS=",chosen.samples[SRX==srx,NREADS],"\n"
                     )
    cat(srr.mk,file="srx_info.mk.new")
    if (!file.exists("srx_info.mk") ||
        !(md5sum("srx_info.mk")==md5sum("srx_info.mk.new"))) {
        file.rename("srx_info.mk.new","srx_info.mk")
    } else {
        file.remove("srx_info.mk.new")
    }
    setwd(wd)
}
