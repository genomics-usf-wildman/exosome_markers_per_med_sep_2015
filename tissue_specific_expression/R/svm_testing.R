library(data.table)
library(caret)
library(reshape2)
library(doMC)
registerDoMC(cores=12)

args <- c("svm_classification",
          "R/svm_training.R",
          "svm_test_results"
          )

args <- commandArgs(trailingOnly=TRUE)

load(args[1])
source(args[2])

trained.object <- NULL
eval(parse(text=paste0("trained.object <- ",caret.run[["object_name"]])))

predictions.caret <- list()
### test the predictions using SVM
for (group in names(trained.object)) {
    predictions.caret[[group]] <-
        list()
    predictions.caret[[group]][["training"]] <- 
        predict(trained.object[[group]],
                newdata=genes.training,
                type="prob")
    predictions.caret[[group]][["training"]]$prediction <-
        colnames(predictions.caret[[group]])[apply(predictions.caret[[group]],
                                                   1,
                                                   which.max)]
    predictions.caret[[group]][["training"]]$actual <-
        genes.training[,group]
    predictions.caret[[group]][["testing"]] <-
        predict(trained.object[[group]],
                newdata=genes.training,
                type="prob")
    predictions.caret[[group]][["testing"]]$prediction <-
        colnames(predictions.caret[[group]])[apply(predictions.caret[[group]],
                                                   1,
                                                   which.max)]
    predictions.caret[[group]][["testing"]]$actual <-
        genes.testing[,group]
}

save.env <- new.env()
save.env[[caret.run[["predictions"]]]] <- predictions.caret

save(list=names(save.env),
     env=save.env,
     file=args[length(args)])

