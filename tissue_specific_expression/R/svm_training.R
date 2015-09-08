
caret.run <-
    list(run_name="svm",
         run_method="svmRadial",
         run_control=trainControl(method="repeatedcv",number=10,repeats=10,classProbs=TRUE),
         tune_length=8,
         object_name="trained.caret.knn"
         )
