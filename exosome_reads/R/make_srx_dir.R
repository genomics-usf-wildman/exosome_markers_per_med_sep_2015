library("tools")
library("data.table")
library("yaml")
args <- c("../tissue_expression/chosen_tissues",
          "../tissue_specific_expression/categorized_samples",
          "chosen_tissues.yaml")
args <- commandArgs(trailingOnly=TRUE)

load(args[1])
load(args[2])

chosen.tissues <- data.table(tissues=yaml.load_file(args[3])$tissues)

## SRR786752 claims to be paired, but only has one read
chosen.samples[SRX=="SRX252261",SRR:=gsub(",SRR786752","",SRR)]

for (srx in categorized.samples[Sample_Group %in% chosen.tissues[,tissues] &
                                grepl("^SRX",SRX),][,list(SRX=SRX[1]),by=Sample_Group][,SRX]) {
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
