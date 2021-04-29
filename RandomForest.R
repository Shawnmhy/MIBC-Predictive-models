### This script is used to train RF classifier
### Author: Haoyang Mi
### @Johns Hopkins University
library(e1071)
library(caret)
library(randomForest)
library(dplyr)
library(ROCR)
library(glmnet)
library(LiblineaR)
library(mRMRe)
library(MLmetrics)

setwd("~/Desktop/TMA_model")

#------------------------------------------------------------#
#-------------- Module 1: Data preparation ------------------#
#------------------------------------------------------------#



#--------------------------#

clinicopathFeature <- Data.matrix[,1:2]

#--------------------------#
source('feature_preprocessing.R')
source('~/Desktop/TMA_963/FunctionHub.R')





# Replace NA with 0
allFeatures_forModel[is.na(allFeatures_forModel)] <- 0

# remove CD8/FoxP3
allFeatures_forModel <- allFeatures_forModel[,-1]



#------------- Image contexture feature ------------# allFeatures_forModel.contexture is right after this

#------------- Image contexture feature ------------# allFeatures_forModel.contexture is right after this
allFeatures_forModel.contexture <- allFeatures_forModel[, 1:1011]
allFeatures_forModel.contexture$Outcome <- allFeatures_forModel$Outcome

#------------- shape descriptor feature ------------#

allFeatures_forModel.shape <- allFeatures_forModel[, 1012:1066]
allFeatures_forModel.shape$Outcome <- allFeatures_forModel$Outcome

#------------- clustering feature -------------------#

allFeatures_forModel.cluster <- allFeatures_forModel[, 1067:1118]
allFeatures_forModel.cluster$Outcome <- allFeatures_forModel$Outcome


#------------- clustering feature -------------------#

allFeatures_forModel.spatstat <- allFeatures_forModel[, 1119:1136]
allFeatures_forModel.spatstat$Outcome <- allFeatures_forModel$Outcome



#------------- cell classification feature -------------------#

allFeatures_forModel.cellClass <- allFeatures_forModel[, 1137:1185]
allFeatures_forModel.cellClass$Outcome <- allFeatures_forModel$Outcome
#------------- vallued clinicopathology feature -------------------#

#------------- vallued clinicopathology feature -------------------#

allFeatures_forModel.patho <- cbind(clinicopathFeature, allFeatures_forModel)






#----------------------- OBTAIN 5 fold accuracy ----------------------#


allFeatures_forModel.ShapeClusTextureClass <- cbind(allFeatures_forModel.shape, allFeatures_forModel.cluster, allFeatures_forModel.cellClass, allFeatures_forModel.contexture)
allFeatures_forModel.1234 <- allFeatures_forModel.ShapeClusTextureClass

#for(feature_count in 3:60){
allFeatures_forModel.234 <- cbind(allFeatures_forModel.shape, allFeatures_forModel.cluster, allFeatures_forModel.contexture)
allFeatures_forModel.2345 <- cbind(allFeatures_forModel.shape, allFeatures_forModel.cluster, allFeatures_forModel.spatstat, allFeatures_forModel.contexture)
allFeatures_forModel.all <- cbind(allFeatures_forModel.shape, allFeatures_forModel.cluster, allFeatures_forModel.spatstat, allFeatures_forModel.cellClass, allFeatures_forModel.contexture)

#allFeatures_12345 <- allFeatures_forModel.contexture

Data.matrix <- allFeatures_forModel.2345
iterations <- 0
  
  # cross validation using bootstrapping
  
  boostrap_length = 100 # define boostap length
  
  Partitioned.data <- createDataPartition(Data.matrix$Outcome, times = boostrap_length, p = 0.2, list = T)
  
  inner.cv <- matrix(nrow = 0, ncol = 3)
  for (id in 1: length(Partitioned.data)){
    
  #  id <- 3
    
    
    x = Partitioned.data[[id]]
    
    # in the next two lines we will separate the Training set into it's 10 pieces
    training_fold = Data.matrix[-x,] # training fold =  training set minus (-) it's sub test fold
    
    test_fold = Data.matrix[x,] # here we describe the test fold individually
    
    
    
    # mRMR data frame
    training_fold$Outcome <- as.numeric(as.character(training_fold$Outcome))
    
    data <- mRMR.data(data.frame(training_fold))
    
    # assign outcome
    
    
    #mRMR.df <- mRMR.classic('mRMRe.filter', data = data, target_indices = ncol(training_fold), feature_count = 40)
    mRMR.df <- mRMR.ensemble(data = data, target_indices = ncol(training_fold), solution_count = 1, feature_count = 40)
    
    mRMR.features <- as.numeric(unlist(mRMR.df@filters))
    
    
    names.vec <- colnames(training_fold[mRMR.features])
    
    
    
    if(length(names.vec) != 0){
      
      optimized.feature[id] <- list(names.vec)
      feature.selected.mRMR <- training_fold[names.vec]
      feature.selected.mRMR$Outcome <- training_fold$Outcome
      
      # mRMR test_fold
      test_fold.mRMR <- test_fold[names.vec]
      test_fold.mRMR$Outcome <- test_fold$Outcome
      
      feature.selected.mRMR$Outcome <- as.factor(feature.selected.mRMR$Outcome)
     # for(ntree in seq(from = 50, to = 500, length = 10)){
      #  for (mtry in seq(from = 5, to = 50, length = 5)) {
          classifier = randomForest(Outcome ~ .,
                                    data = feature.selected.mRMR,
                                    #maxnodes = maxnode,
                                    ntree = 400,
                                    mtry = 5,
                                    #nodesize = nodesize,
                                    importance = TRUE)     
          # next step in the loop, we calculate the predictions and cm and we equate the accuracy
          y_pred = predict(classifier, newdata = test_fold.mRMR[names.vec],  probability = TRUE) #probability = TRUE
          
          AUC.metric <- AUC(y_pred = y_pred, y_true = test_fold.mRMR$Outcome)
          
          F1.metric <- F1_Score(y_pred = y_pred, y_true = test_fold.mRMR$Outcome, positive = '1')
          
          
          cm = table(pred = y_pred, true=test_fold.mRMR$Outcome)
          accuracy = sum(diag(cm))/sum(cm)
          
          #inner.cv <- rbind(inner.cv, cbind(gamma, cost, accuracy, AUC.metric))
          #inner.cv <- rbind(inner.cv, cbind(ntree, mtry, accuracy, AUC.metric, F1.metric))
          inner.cv <- rbind(inner.cv, cbind(accuracy, AUC.metric, F1.metric))
          
      
      iterations <- iterations + 1
      
      print('------------------')
      print(paste('finish iteration ', iterations, sep = ''))
    }
  }
  #print(mean(optimization.folds.AUC[,3]))
 #mRMR.accuracy <- rbind(mRMR.accuracy, cbind(feature_n, mean(optimization.folds[,3])))
  
#}
  
  id <- 1
  mean(as.numeric(as.character(inner.cv[, id])))
  sd(as.numeric(as.character(inner.cv[, id])))

  
  list <- optimized.feature
  selected.f.1234 <- matrix(nrow = 0, ncol = 1)
  for(run.id in 1:boostrap_length){
    
    # get the features that are selected by this run
    if(length(list[[run.id]]) > 0){
      selected.f.1234 <- rbind(selected.f.1234, as.matrix(list[[run.id]]))
    }
    
  }
  count_ <- data.frame(table(selected.f.1234))
  count_
  
  # top 30 features
  top30Feature <- count_[order(count_$Freq, decreasing = TRUE), ]
  top30Feature <- as.character(top30Feature[1:30,1])
  
    
#print(mean(optimization.folds.ACC[,3]))
#print(mean(optimization.folds.AUC[,3]))
#print(mean(optimization.folds.F1[,3]))

#print(sd(optimization.folds.ACC[,3]))
#print(sd(optimization.folds.AUC[,3]))
#print(sd(optimization.folds.F1[,3]))






#----------------------- Baseline, Random Forest, test ----------------------#

CP <- read.csv('./Clinicopathologic Features.MOD.csv')
CP$CHEMO.CYCLES.MOD <- as.numeric(as.character(CP$CHEMO.CYCLES.MOD))
CP <- CP[CP$CHEMO.CYCLES.MOD > 2 & CP$TMA,]
colnames(CP)[1] <- 'patientID'



# Select features for baseline, including Clinical T Stage
allFeatures.base.num <- mutate_all(CP[,c(3:17,19)], function(x) as.numeric(x))

allFeatures.base <- cbind(CP$patientID, allFeatures.base.num, CP$X.Nonresponder.ypT2.3.4....1)


allFeatures.base <- allFeatures.base[complete.cases(allFeatures.base),]



# categorize AGE and Clinical T Stage

allFeatures.base[allFeatures.base$AGE.AT.OPERATION < 60, 'AGE.AT.OPERATION'] <- 0
allFeatures.base[allFeatures.base$AGE.AT.OPERATION >= 60, 'AGE.AT.OPERATION'] <- 1
allFeatures.base[allFeatures.base$CLINICAL.T.STAGE <= 2, 2] <- 0
allFeatures.base[allFeatures.base$CLINICAL.T.STAGE > 2, 2] <- 1



colnames(allFeatures.base)[1] <- 'patientID'
colnames(allFeatures.base)[18] <- 'Outcome'



allFeatures.base <- allFeatures.base[,-1]

iterations <- 0

# machine learning

optimization.folds.F1 <- matrix(nrow = 0, ncol = 3)
optimization.folds.ACC <- matrix(nrow = 0, ncol = 3)
optimization.folds.AUC <- matrix(nrow = 0, ncol = 3)
#optimized.feature <- list()    


optimized.feature <- list()



# cross validation using bootstrapping

boostrap_length = 100 # define boostap length
Partitioned.data <- createDataPartition(allFeatures.base$Outcome, times = boostrap_length, p = 0.2, list = T)


for (id in 1: length(Partitioned.data)){
  
  Partitioned.data <- createDataPartition(allFeatures.base$Outcome, times = boostrap_length, p = 0.2, list = T)
  inner.cv <- matrix(nrow = 0, ncol = 3)
  
  
  #id <- 5
  x = Partitioned.data[[id]]
  
  # in the next two lines we will separate the Training set into it's 10 pieces
  training_fold = allFeatures.base[-x,] # training fold =  training set minus (-) it's sub test fold
  training_fold$Outcome <- as.factor(training_fold$Outcome)
  test_fold = allFeatures.base[x,] # here we describe the test fold individually
  
  
  #inner.cv <- matrix(nrow = 0, ncol = 3)
  #for(maxnode in 2:10){
  #  for (nodesize in 2:35) {
      classifier = randomForest(Outcome ~ .,
                                data = training_fold,
                                #maxnodes = maxnode,
                                ntree = 400,
                                mtry = 5,
                                #nodesize = nodesize,
                                importance = TRUE)     
      # next step in the loop, we calculate the predictions and cm and we equate the accuracy
      y_pred = predict(classifier, newdata = test_fold[,-ncol(test_fold)],  probability = TRUE) #probability = TRUE
      
      AUC.metric <- AUC(y_pred = y_pred, y_true = test_fold$Outcome)
      
      F1.metric <- F1_Score(y_pred = y_pred, y_true = test_fold$Outcome, positive = '1')
      
      
      cm = table(pred = y_pred, true=test_fold$Outcome)
      accuracy = sum(diag(cm))/sum(cm)
      
      #inner.cv <- rbind(inner.cv, cbind(gamma, cost, accuracy, AUC.metric))
      inner.cv <- rbind(inner.cv, cbind(accuracy, AUC.metric, F1.metric))
      
    #}
  #}
  optimization.folds.ACC <- rbind(optimization.folds.ACC, inner.cv[which.max(inner.cv[,1]),1])
  optimization.folds.AUC <- rbind(optimization.folds.AUC, inner.cv[which.max(inner.cv[,2]),2])
  optimization.folds.F1 <- rbind(optimization.folds.F1, inner.cv[which.max(inner.cv[,3]),3])
  #print(inner.cv[which.max(inner.cv[,3]),])
  #print(Accuracy.folds)
  
  iterations <- iterations + 1
  
  print('------------------')
  print(paste('finish iteration ', iterations, sep = ''))
}

# print mean and sd
print(mean(optimization.folds.ACC[,3]))
print(mean(optimization.folds.AUC[,3]))
print(mean(optimization.folds.F1[,3]))

print(sd(optimization.folds.ACC[,3]))
print(sd(optimization.folds.AUC[,3]))
print(sd(optimization.folds.F1[,3]))


