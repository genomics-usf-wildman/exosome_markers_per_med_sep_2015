library(data.table)

args <- commandArgs(trailingOnly=TRUE)

### load chosen_samples
load(args[1])

mk.output <-
    paste0("SRX_FILES=",paste(collapse=" ",chosen.samples[,SRX]),"\n")
for(srx in chosen.samples[,SRX]) {
    mk.output <- paste0(mk.output,
                        srx,"_SRR=",paste(collapse=" ",chosen.samples[SRX==srx,SRR]),"\n")
}
cat(mk.output,
    file=args[length(args)])
