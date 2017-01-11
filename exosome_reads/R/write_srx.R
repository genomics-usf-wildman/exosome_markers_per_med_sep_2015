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

args <- commandArgs(trailingOnly=TRUE)

### load chosen_samples
load(args[1])

srxes <-
    categorized.samples[Sample_Group %in% chosen.tissues[,tissues] &
                        grepl("^SRX",SRX),][,list(SRX=SRX[1]),by=Sample_Group][,SRX]
mk.output <-
    paste0("SRX_FILES=",paste(collapse=" ",srxes),"\n")
for(srx in srxes) {
    mk.output <- paste0(mk.output,
                        srx,"_SRR=",paste(collapse=" ",chosen.samples[SRX==srx,SRR]),"\n")
}
cat(mk.output,
    file=args[length(args)])
