library(data.table)
library(caret)
library(reshape2)
library(doMC)
registerDoMC(cores=12)

args <- c("interesting_gene_reads",
          "svm_classification"
          )

args <- commandArgs(trailingOnly=TRUE)

load(args[1])


gene.fpkm.wide <-
    dcast(interesting.gene.reads,
          Sample_Group+SRX~gene_short_name,
          value.var="FPKM",
          fun.aggregate=function(x){x[1]})

### for the time being, we're interested in uterus, placenta, or neither.
gene.fpkm.wide$ut.pla.none <- as.character(gene.fpkm.wide$Sample_Group)
gene.fpkm.wide$ut.pla.none[!(gene.fpkm.wide$ut.pla.none %in% c("uterus","placenta"))] <-
    "neither"
gene.fpkm.wide$ut.pla.none <-
    factor(gene.fpkm.wide$ut.pla.none)

## set NA to zero
gene.fpkm.wide[is.na(gene.fpkm.wide)] <- 0

set.seed(10609)
in_training <-
    createDataPartition(gene.fpkm.wide$ut.pla.none,
                        p = .75,
                        list = FALSE)

genes.training <-
    gene.fpkm.wide[in_training,]

genes.testing <-
    gene.fpkm.wide[-in_training,]

knn.train <-
    train(form=ut.pla.none~.,
          data=genes.training[,-(1:2)],
          method="knn",
          trControl=trainControl(method="cv"),
          tuneLength=5
          )


svm.train <-
    train(form=ut.pla.none~.,
          data=genes.training[,-(1:2)],
          method="svmRadial",
          trControl=trainControl(method="cv",number=10)
          )


save(knn.train,
     svm.train,
     file=args[length(args)])

