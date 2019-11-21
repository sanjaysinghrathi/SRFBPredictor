SRFBPredictor <- function(training_data_exp_matrix, training_data_phenodata, testing_data_exp_matrix, testing_data_phenodata, biomarker) {
  # Predict response for a colorectal cancer patient given Chemoradiotherapy
  
  # Load training data
  training_em <- training_data_exp_matrix
  training_pData <- training_data_phenodata
  colnames(training_pData) <- c("response", colnames(training_pData)[2:dim(training_pData)[2]])
  if(dim(training_em)[2] != dim(training_pData)[1]){
    stop("Error: Number of samples in training data expression matrix and phenodata do not match")
  }
  
  # Load testing data
  testing_em <- testing_data_exp_matrix
  testing_pData <- testing_data_phenodata
  colnames(testing_pData) <- c("response", colnames(testing_pData)[2:dim(testing_pData)[2]])
  if(dim(testing_em)[2] != dim(testing_pData)[1]){
    stop("Error: Number of samples in testing data expression matrix and phenodata do not match")
  }
  
  # Check same column for testing and training phenodata
  if(dim(testing_pData)[2] != dim(training_pData)[2]){
    stop("Error: Number of columns in training and testing phenodata are not same")
  }
  
  # Remove batch effects using SVA ComBat
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")}
  
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
  plot3 <- autoplot(pcDat,
           data = rbind(training_pData, testing_pData), 
           colour="response", 
           shape="response",
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
  set.seed(123)
  model <- train(response ~., data = training_data, method = "glmnet", trControl = trainControl(method = "LOOCV", allowParallel = TRUE), tuneLength = 10)
  model
  internal <- max(model$results$Accuracy)
  
  # Make Predictions
  if (!requireNamespace("dplyr", quietly = TRUE)){
    BiocManager::install("dplyr")}
  library(dplyr)
  predictions <- model %>% predict(testing_data)
  predictions
  
  pred <- confusionMatrix(as.factor(predictions),as.factor(testing_data$response))
  return(list(Train_PCA=plot1, Test_PCA= plot2, TrainTest_PCA= plot3, Model_Internal_Accuracy=internal, Model_Prediction_Accuracy=pred))
}

