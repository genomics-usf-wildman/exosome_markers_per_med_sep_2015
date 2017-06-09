library("data.table")
library("caret")
library("reshape2")
library("doMC")
registerDoMC(cores=8)

args <- c("R/svm_training.R",
          "interesting_gene_reads",
          "svm_classification"
          )

args <- commandArgs(trailingOnly=TRUE)

### this file provides the caret.run object
source(args[1])
load(args[2])


gene.fpkm.wide <-
    dcast(interesting.gene.reads,
          Sample_Group+SRX~tracking_id,
          value.var="FPKM",
          fun.aggregate=function(x){x[1]})
## set NA to zero
gene.fpkm.wide[is.na(gene.fpkm.wide)] <- 0

## ignore genes with maximum expression less than 100 FPKM.
gene.max <- apply(gene.fpkm.wide,2,function(x){max(x,na.rm=TRUE)})
gene.fpkm.wide <- gene.fpkm.wide[,c(1,2,which(as.numeric(gene.max)>100))]
gene.names <- colnames(gene.fpkm.wide)[-(1:2)]

### for the time being, we're interested in uterus, placenta, or neither.
gene.fpkm.wide$ut.pla.none <- as.character(gene.fpkm.wide$Sample_Group)
gene.fpkm.wide$ut.pla.none[!(gene.fpkm.wide$ut.pla.none %in% c("uterus","placenta"))] <-
    "neither"
gene.fpkm.wide$ut.pla.none <-
    factor(gene.fpkm.wide$ut.pla.none)

for (group in c("uterus","placenta")) {
    gene.fpkm.wide[[group]] <-
        as.character(gene.fpkm.wide$Sample_Group)
    gene.fpkm.wide[[group]][!(gene.fpkm.wide[[group]] %in% group)] <-
        paste0("not_",group)
    gene.fpkm.wide[[group]] <-
        as.factor(gene.fpkm.wide[[group]])
}

### seed generated with echo $RANDOM
set.seed(10609)
in_training <-
    createDataPartition(gene.fpkm.wide$ut.pla.none,
                        p = .75,
                        list = FALSE)
genes.training <-
    gene.fpkm.wide[in_training,]

genes.testing <-
    gene.fpkm.wide[-in_training,]

trained.caret <- list()
### multi-class ROC https://www.kaggle.com/c/forest-cover-type-prediction/forums/t/10573/r-caret-using-roc-instead-of-accuracy-in-model-training
for (group in c("uterus","placenta")) {
    ## seed generated with echo $RANDOM;
    set.seed(29224)
    data.scaled <- scale(genes.training[,(gene.names)])
    data.scaled[!is.finite(data.scaled)] <- 0

    trained.caret[[group]] <-
        train(form=as.formula(paste0(group,"~.")),
              data=data.frame(genes.training[,group,drop=FALSE],
                  data.scaled),
              method=caret.run[["run_method"]],
              trControl=caret.run[["run_control"]],
              tuneLength=caret.run[["tune_length"]],
              metric="ROC"
              )
    message(paste0("finished ",group))
}
save.env <- new.env()
save.env[[caret.run[["object_name"]]]] <- trained.caret
save.env[["genes.training"]] <- genes.training
save.env[["genes.testing"]] <- genes.testing
save.env[["gene.names"]] <- gene.names

save(list=names(save.env),
     envir=save.env,
     file=args[length(args)])

