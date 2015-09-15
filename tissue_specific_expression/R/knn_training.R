caret.run <-
    list(run_name="knn",
         run_method="knn",
         run_control=trainControl(method="repeatedcv",number=5,repeats=3,classProbs=TRUE),
         tune_length=3,
         object_name="trained.caret.knn",
         predictions="predictions.caret.knn"
         )

registerDoMC(cores=2)
