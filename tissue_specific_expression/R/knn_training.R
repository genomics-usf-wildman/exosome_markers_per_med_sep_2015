caret.run <-
    list(run_name="knn",
         run_method="knn",
         run_control=trainControl(method="repeatedcv",number=10,repeats=10,classProbs=TRUE),
         tune_length=8,
         object_name="trained.caret.knn",
         predictions="predictions.caret.knn"
         )

registerDoMC(cores=2)
