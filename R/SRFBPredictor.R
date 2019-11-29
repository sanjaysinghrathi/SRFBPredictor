SRFBPredictor <- function(training_data_exp_matrix, training_data_phenodata, testing_data_exp_matrix, testing_data_phenodata, biomarker, combat=TRUE, mlmodel="elastic-net") {
  # Predict response for a colorectal cancer patient given Chemoradiotherapy
  if (!requireNamespace("dplyr", quietly = TRUE)){
    BiocManager::install("dplyr")}
  library(dplyr)
  
  # Load training data
  training_em <- training_data_exp_matrix
  training_pData <- training_data_phenodata
  colnames(training_pData) <- c("response", colnames(training_pData)[2:dim(training_pData)[2]])
  if(dim(training_em)[2] != dim(training_pData)[1]){
    stop("Error: Number of samples in training data expression matrix and phenodata do not match")
  }
  if(all(rownames(training_pData)%in% colnames(training_em))== FALSE){
    stop("Error: Names of samples in training data expression matrix and phenodata do not match")
  }
  if(all(rownames(training_pData) == colnames(training_em))== FALSE){
    training_pData <- training_pData[colnames(training_em),]
  }

  # Load testing data
  testing_em <- testing_data_exp_matrix
  testing_pData <- testing_data_phenodata
  colnames(testing_pData) <- c("response", colnames(testing_pData)[2:dim(testing_pData)[2]])
  if(dim(testing_em)[2] != dim(testing_pData)[1]){
    stop("Error: Number of samples in testing data expression matrix and phenodata do not match")
  }
  if(all(rownames(testing_pData)%in% colnames(testing_em))== FALSE){
    stop("Error: Names of samples in testing data expression matrix and phenodata do not match")
  }
  if(all(rownames(testing_pData) == colnames(testing_em))== FALSE){
    testing_pData <- testing_pData[colnames(testing_em),]
  }

  # Check same column for testing and training phenodata
  if(dim(testing_pData)[2] != dim(training_pData)[2]){
    stop("Error: Number of columns in training and testing phenodata are not same")
  }

  # Remove batch effects using SVA ComBat
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")}
  # if(exists("combat")==FALSE){
  #   combat <- TRUE
  # }
  
  if(combat==TRUE){
    if (!requireNamespace("sva", quietly = TRUE)){
      BiocManager::install("sva")}
    
    library(sva)
    full_em <- cbind(training_em,testing_em)
    full_pData <- rbind(training_pData, testing_pData)
    batch <- c(rep(1, dim(training_em)[2]),rep(2, dim(testing_em)[2]))
    mod <- model.matrix(~as.factor(response), data=full_pData)
    combat_edata3 <- ComBat(dat=full_em, batch=batch, mod = mod, par.prior = TRUE, ref.batch=1)
    training_em_new <- combat_edata3[,1:dim(training_em)[2]]
    testing_em_new <- combat_edata3[, (dim(training_em)[2]+1):dim(full_em)[2]]
  }else{
    training_em_new <- training_em
    testing_em_new <- testing_em
  }

  # PCA after Decision making genes discovery in training dataset
  if (!requireNamespace("ggfortify", quietly = TRUE)){
    BiocManager::install("ggfortify")}
  library(ggfortify)

  ee <- training_em_new
  pcDat <- prcomp(t(ee[biomarker,]), scale. = TRUE)
  plot1 <-autoplot(pcDat,
           data = training_pData,
           colour="response",
           shape="response",
           label=FALSE,
           main="PCA for training data w.r.t. response",
           size=2)

  ## cleaning
  rm(ee)
  rm(pcDat)

  # PCA after Decision making genes discovery in testing set
  ee <- testing_em_new
  pcDat <- prcomp(t(ee[biomarker,]), scale. = TRUE)
  plot2 <- autoplot(pcDat,
           data = testing_pData,
           colour="response",
           shape="response",
           label=FALSE,
           main="PCA for testing data w.r.t. response",
           size=2)

  ## cleaning
  rm(ee)
  rm(pcDat)

  #' PCA after Decision making genes discovery in traing+testing set
  ee <- cbind(training_em_new, testing_em_new)
  pcDat <- prcomp(t(ee[biomarker,]), scale. = TRUE)
  batch <- c(rep("Train", dim(training_em)[2]),rep("Test", dim(testing_em)[2]))
  plot3 <- autoplot(pcDat,
           data = cbind(rbind(training_pData, testing_pData), batch),
           colour="response",
           shape="batch",
           label=FALSE,
           main="PCA for training + testing data w.r.t. response",
           size=2)

  ## cleaning
  rm(ee)
  rm(pcDat)

  # Load Signature
  biomarker <- biomarker

  # Prepare Prediction Model

  # Prepare training set
  training_data <- as.data.frame(training_em_new[as.character(biomarker),])
  training_data <- cbind(t(training_data), data.frame(training_pData[,c(dim(training_pData)[2]:1)]))
  training_data$response <- as.factor(as.character(training_data$response))

  # Prepare testing data
  #testing_data <- as.data.frame(testing_em[as.character(biomarker),])
  testing_data <- as.data.frame(testing_em_new[as.character(biomarker),])
  testing_data <- cbind(t(testing_data), data.frame(testing_pData[,c(dim(testing_pData)[2]:1)]))
  testing_data$response <- as.factor(as.character(testing_data$response))

  # Train elastic net model on training data
  if (!requireNamespace("e1071", quietly = TRUE)){
    BiocManager::install("e1071")}
  library(e1071)
  if (!requireNamespace("caret", quietly = TRUE)){
    BiocManager::install("caret")}
  library(caret)
  
  # if(exists("mlmodel")==FALSE){
  #   mlmodel <- "elastic-net"
  # }
  
  if(mlmodel=="elastic-net"){
    set.seed(123)
    garbage <- capture.output(suppressWarnings(model <- train(response ~., data = training_data, method = "glmnet", trControl = trainControl(method = "LOOCV", allowParallel = TRUE), tuneLength = 10)))
    model
    internal <- max(model$results$Accuracy, na.rm = TRUE)
    
    predictions <- predict(model, newdata = testing_data, type = "prob")
    row.names(predictions) <- row.names(testing_data)
    predictions_prob <- data.frame(predictions)
    
    predictions <- predict(model, newdata = testing_data)
    predictions
    predictions_prob$Predicted_class <- predictions
    predictions_prob$Actual_class <- testing_data$response
    
    if(nlevels(as.factor(predictions))==nlevels(as.factor(testing_data$response))){
      confusionMatrix(predictions, testing_data$response )
      pred <- confusionMatrix(predictions, testing_data$response )
    }else{
      pred <- predictions
    }
  }
  if(mlmodel == "svmLinear"){
    trctrl <- trainControl(method = "LOOCV", allowParallel = TRUE, classProbs =  TRUE)
    grid <- expand.grid(C = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2,5))
    set.seed(2345)
    garbage <- capture.output(suppressWarnings(svm_Linear_Grid <- train(response ~., data = training_data, method = "svmLinear", trControl=trctrl, preProcess = c("center", "scale"), metric="Accuracy", tuneGrid = grid, tuneLength = 10)))
    svm_Linear_Grid
    internal <- max(svm_Linear_Grid$results$Accuracy, na.rm = TRUE)
    
    predictions <- predict(svm_Linear_Grid, newdata = testing_data, type = "prob")
    predictions
    row.names(predictions) <- row.names(testing_data)
    predictions_prob <- data.frame(predictions)
    
    predictions <- predict(svm_Linear_Grid, newdata = testing_data)
    predictions
    predictions_prob$Predicted_class <- predictions
    predictions_prob$Actual_class <- testing_data$response
    
    if(nlevels(as.factor(predictions))==nlevels(as.factor(testing_data$response))){
      confusionMatrix(predictions, testing_data$response )
      pred <- confusionMatrix(predictions, testing_data$response )
    }else{
      pred <- predictions
    }
  }
  if(mlmodel == "neuralNet"){
    numFolds <- trainControl(method = 'LOOCV', allowParallel = TRUE, verbose=TRUE , search = "grid")
    grid <- expand.grid(size=c(seq(from = 1, to = 10, by = 1)),
                        decay=c(seq(from = 0.0, to = 0.5, by = 0.1)))
    
    set.seed(567)
    garbage <- capture.output(suppressWarnings(model <- train(response ~ ., training_data, method='nnet', trace = FALSE, preProcess = c('center', 'scale'), metric="Accuracy", trControl = numFolds, linout=FALSE, tuneGrid=grid)))
    model
    internal <- max(model$results$Accuracy, na.rm = TRUE)
    
    predictions <- predict(model, newdata = testing_data, type = "prob")
    row.names(predictions) <- row.names(testing_data)
    predictions_prob <- data.frame(predictions)
    
    predictions <- predict(model, newdata = testing_data)
    predictions
    predictions_prob$Predicted_class <- predictions
    predictions_prob$Actual_class <- testing_data$response
    
    if(nlevels(as.factor(predictions))==nlevels(as.factor(testing_data$response))){
      confusionMatrix(predictions, testing_data$response )
      pred <- confusionMatrix(predictions, testing_data$response )
    }else{
      pred <- predictions
    }
  }
  if(mlmodel == "randomForest"){
    if(!requireNamespace("randomForest", quietly = TRUE)){
      BiocManager::install("randomForest")}
    library(randomForest)
    set.seed(567)
    garbage <- capture.output(suppressWarnings(bestmtry <- as.data.frame(tuneRF(training_data[,1:(dim(training_data)[2]-1)], training_data$response, stepFactor=1.5, improve=1e-5, ntree=2000))))
    print(bestmtry)
    numFolds <- trainControl(method = 'LOOCV', allowParallel = TRUE, verbose=TRUE , search = "grid")
    grid <- expand.grid(.mtry=bestmtry$mtry)
    
    set.seed(567)
    garbage <- capture.output(suppressWarnings(model <- train(response ~ ., training_data, method='rf', preProcess = c('center', 'scale'), metric="Accuracy", trControl = numFolds, trace=FALSE, linout=FALSE, tuneGrid=grid)))
    model
    internal <- max(model$results$Accuracy, na.rm = TRUE)
    
    predictions <- predict(model, newdata = testing_data, type = "prob")
    row.names(predictions) <- row.names(testing_data)
    predictions_prob <- data.frame(predictions)
    
    predictions <- predict(model, newdata = testing_data)
    predictions
    predictions_prob$Predicted_class <- predictions
    predictions_prob$Actual_class <- testing_data$response
    
    if(nlevels(as.factor(predictions))==nlevels(as.factor(testing_data$response))){
      confusionMatrix(predictions, testing_data$response )
      pred <- confusionMatrix(predictions, testing_data$response )
    }else{
      pred <- predictions
    }
  }
  if(mlmodel == "svmNonLinear"){
    numFolds <- trainControl(method = 'LOOCV', allowParallel = TRUE, verbose=TRUE , search = "grid", classProbs = TRUE)
    grid <- expand.grid(sigma = c(0,0.01, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.07,0.08, 0.09, 0.1, 0.25, 0.5, 0.75,0.9),
                        C = c(0,0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2,5))
    set.seed(3233)
    garbage <- capture.output(suppressWarnings(model <- train(response ~ ., training_data, method='svmRadial', trace = FALSE, preProcess = c('center', 'scale'), metric="Accuracy", trControl = numFolds, linout=FALSE, tuneGrid=grid)))
    model
    internal <- max(model$results$Accuracy, na.rm = TRUE)
    
    predictions <- predict(model, newdata = testing_data, type = "prob")
    row.names(predictions) <- row.names(testing_data)
    predictions_prob <- data.frame(predictions)
    
    predictions <- predict(model, newdata = testing_data)
    predictions
    predictions_prob$Predicted_class <- predictions
    predictions_prob$Actual_class <- testing_data$response
    
    if(nlevels(as.factor(predictions))==nlevels(as.factor(testing_data$response))){
      confusionMatrix(predictions, testing_data$response )
      pred <- confusionMatrix(predictions, testing_data$response )
    }else{
      pred <- predictions
    }
  }
  return(list(Train_PCA=plot1, Test_PCA= plot2, TrainTest_PCA= plot3, Model_Internal_Accuracy=internal, Model_Prediction_Probabilities=predictions_prob, Model_Prediction_Accuracy=pred))
}

