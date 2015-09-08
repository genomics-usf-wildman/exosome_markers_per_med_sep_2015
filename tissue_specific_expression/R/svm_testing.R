library(data.table)
library(caret)
library(reshape2)
library(doMC)
registerDoMC(cores=12)

args <- c("svm_classification",
          "svm_test_results"
          )

args <- commandArgs(trailingOnly=TRUE)

load(args[1])


### test the predictions using SVM
gene.testing.predictions.svm.prob <-
    predict(svm.train,
            newdata=genes.testing,
            type="prob")
gene.testing.predictions.svm.prob$prediction <-
    colnames(gene.testing.predictions.svm.prob)[apply(gene.testing.predictions.svm.prob,1,which.max)]

gene.training.predictions.svm.prob <-
    predict(svm.train,
            newdata=genes.training,
            type="prob")
gene.training.predictions.svm.prob$prediction <-
    colnames(gene.training.predictions.svm.prob)[apply(gene.training.predictions.svm.prob,1,which.max)]

### test the predictions using the KNN
gene.testing.predictions.knn.prob <-
    predict(knn.train,
            newdata=genes.testing,
            type="prob")
gene.testing.predictions.knn.prob$prediction <-
    colnames(gene.testing.predictions.knn.prob)[apply(gene.testing.predictions.knn.prob,1,which.max)]

gene.training.predictions.knn.prob <-
    predict(knn.train,
            newdata=genes.training,
            type="prob")
gene.training.predictions.knn.prob$prediction <-
    colnames(gene.training.predictions.knn.prob)[apply(gene.training.predictions.knn.prob,1,which.max)]

genes.testing.actual <- genes.testing[,c("Sample_Group","ut.pla.none")]
genes.training.actual <- genes.training[,c("Sample_Group","ut.pla.none")]

save(gene.testing.predictions.svm.prob,
     gene.training.predictions.svm.prob,
     gene.testing.predictions.knn.prob,
     gene.training.predictions.knn.prob,
     genes.testing.actual,
     genes.training.actual,
     file=args[length(args)])

