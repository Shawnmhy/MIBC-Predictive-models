### This script is used to train SVM classifier
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

source('feature_preprocessing.R')
source('~/Desktop/TMA_963/FunctionHub.R')



allFeatures_forModel <- allFeatures_forModel[,-1]
allFeatures_forModel[is.na(allFeatures_forModel)] <- 0
#write.csv(allFeatures_forModel, 'allFeatures_forModel.csv')


#------------- Image contexture feature ------------# 
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

#------------------------------------------------------------#
#------------ Module 2: Model training and testing ----------#
#------------------------------------------------------------#



#--------------------------------------------------------------#
#----------------------- SVM  Classifier ----------------------#
#--------------------------------------------------------------#



#---------------------------------------------------------------#
# Count the number of apperance for each feature




#----------------------- Shape + Clustering +  Image texture + Cell Classification + Classifier ----------------------#

allFeatures_forModel.ShapeClusTextureClass <- cbind(allFeatures_forModel.shape, allFeatures_forModel.cluster, allFeatures_forModel.cellClass, allFeatures_forModel.contexture)
allFeatures_forModel.1234 <- allFeatures_forModel.ShapeClusTextureClass

#----------------------- Computational features Classifier ----------------------#

allFeatures_forModel.all <- cbind(allFeatures_forModel.spatstat, allFeatures_forModel.cluster, allFeatures_forModel.shape, allFeatures_forModel.cellClass, allFeatures_forModel.contexture)



#----------------------- 2 + 3 + 4 Classifier ----------------------#
allFeatures_forModel.234 <- cbind(allFeatures_forModel.cluster, allFeatures_forModel.shape, allFeatures_forModel.contexture)


#----------------------- 2 + 3 + 4 + 5 Classifier ----------------------#
allFeatures_forModel.2345 <- cbind(allFeatures_forModel.spatstat, allFeatures_forModel.cluster, allFeatures_forModel.shape, allFeatures_forModel.contexture)





mRMR.accuracy <- matrix(nrow = 0, ncol = 2)



#allFeatures_12345 <- allFeatures_forModel.contexture

#Data.matrix <- cbind(Data.matrix_withCli[,1:3], allFeatures_forModel)

#Data.matrix <- Data.matrix[,-c(1:4)]
#Data.matrix <- cbind(clinic, Data.matrix)

Data.matrix <- allFeatures_forModel.all

iterations <- 0
optimized.feature <- list()



# cross validation using bootstrapping

boostrap_length = 100 # define boostap length

Partitioned.data <- createDataPartition(Data.matrix$Outcome, times = boostrap_length, p = 0.2, list = T)

inner.cv <- matrix(nrow = 0, ncol = 5)

optimization.folds.ACC <- data.frame(matrix(nrow = 0, ncol = 0))
optimization.folds.F1 <- data.frame(matrix(nrow = 0, ncol = 0))
optimization.folds.AUC <- data.frame(matrix(nrow = 0, ncol = 0))

for (id in 1: length(Partitioned.data)){
  
  #id <- 6
  #inner.cv <- matrix(nrow = 0, ncol = 5)
  
  x = Partitioned.data[[id]]
  
  # in the next two lines we will separate the Training set into it's 10 pieces
  training_fold = Data.matrix[-x,] # training fold =  training set minus (-) it's sub test fold
  
  test_fold = Data.matrix[x,] # here we describe the test fold individually
  
  
  
  # mRMR data frame
  training_fold$Outcome <- as.numeric(as.character(training_fold$Outcome))
  data <- mRMR.data(data.frame(training_fold))
  #data <- data.frame(target = feature.selected.test$Outcome, feature.selected.test)
  
  # assign outcome
  
  
  mRMR.df <- mRMR.ensemble( data = data, target_indices = ncol(training_fold), solution = 1, feature_count = 45)
  
  mRMR.features <- as.numeric(unlist(mRMR.df@filters))
  
  
  names.vec <- colnames(training_fold[mRMR.features])
  
  
  
  if(length(names.vec) != 0){
    
    optimized.feature[id] <- list(names.vec)
    feature.selected.mRMR <- training_fold[names.vec]
    feature.selected.mRMR$Outcome <- training_fold$Outcome
    
    # mRMR test_fold
    test_fold.mRMR <- test_fold[names.vec]
    test_fold.mRMR$Outcome <- test_fold$Outcome
    
    #tprfpr_df <- matrix(nrow = 0, ncol = 2)
    for(type in c('C-classification', 'nu-classification')){
      for(kernel in c('linear', 'sigmoid', 'polynomial', 'radial')){
        
        #gamma <- 0.01
        #cost <- 1
        
        classifier = svm(formula = Outcome ~ .,
                         data = feature.selected.mRMR,
                         type = type,
                         kernel = kernel,
                         #gamma = 0.001,
                         #cost = 35, 
                         probability = TRUE)      
        # next step in the loop, we calculate the predictions and cm and we equate the accuracy
        y_pred = predict(classifier, newdata = test_fold.mRMR[names.vec],  probability = TRUE) #probability = TRUE
        
        
        
        
        AUC.metric <- AUC(y_pred = y_pred, y_true = test_fold.mRMR$Outcome)
        
        F1.metric <- F1_Score(y_pred = y_pred, y_true = test_fold.mRMR$Outcome, positive = '1')
        
        
        cm = table(pred = y_pred, true=test_fold.mRMR$Outcome)
        accuracy = sum(diag(cm))/sum(cm)
        
        #inner.cv <- rbind(inner.cv, cbind(gamma, cost, accuracy, AUC.metric))
        #inner.cv <- rbind(inner.cv, cbind(gamma, cost, accuracy, AUC.metric, F1.metric))
        inner.cv <- rbind(inner.cv, cbind(type, kernel, accuracy, AUC.metric, F1.metric))
        
      }
    }
    optimization.folds.ACC <- rbind(optimization.folds.ACC, inner.cv[which.max(inner.cv[,3]),c(1,2,3)])
    
    optimization.folds.AUC <- rbind(optimization.folds.AUC, inner.cv[which.max(inner.cv[,4]),c(1,2,4)])
    optimization.folds.F1 <- rbind(optimization.folds.F1, inner.cv[which.max(inner.cv[,5]),c(1,2,5)])
    
    #print(inner.cv[which.max(inner.cv[,4]),])
    #print(Accuracy.folds)
    
    iterations <- iterations + 1
    
    print('------------------')
    print(paste('finish iteration ', iterations, sep = ''))
  }
}
#print(mean(optimization.folds.AUC[,3]))
#mRMR.accuracy <- rbind(mRMR.accuracy, cbind(feature_count, mean(optimization.folds.ACC[,3])))  
#}

#write.csv(optimization.folds, '1234_mRMR_SVM_CHEMCYCLE.csv')


inner.cv <- data.frame(inner.cv)

# which to check? 3: ACC, 4: AUC, 5: F1 score
id <- 3
#mean(as.numeric(as.character(inner.cv[inner.cv$type == 'C-classification' & inner.cv$kernel == 'linear', id])))
mean(as.numeric(as.character(inner.cv[inner.cv$type == 'C-classification' & inner.cv$kernel == 'radial', id])))
#mean(as.numeric(as.character(inner.cv[inner.cv$type == 'C-classification' & inner.cv$kernel == 'sigmoid', id])))
mean(as.numeric(as.character(inner.cv[inner.cv$type == 'C-classification' & inner.cv$kernel == 'polynomial', id])))

#mean(as.numeric(as.character(inner.cv[inner.cv$type == 'nu-classification' & inner.cv$kernel == 'linear', id])))
mean(as.numeric(as.character(inner.cv[inner.cv$type == 'nu-classification' & inner.cv$kernel == 'radial', id])))
#mean(as.numeric(as.character(inner.cv[inner.cv$type == 'nu-classification' & inner.cv$kernel == 'sigmoid', id])))
mean(as.numeric(as.character(inner.cv[inner.cv$type == 'nu-classification' & inner.cv$kernel == 'polynomial', id])))




#sd(as.numeric(as.character(inner.cv[inner.cv$type == 'C-classification' & inner.cv$kernel == 'linear', id])))
sd(as.numeric(as.character(inner.cv[inner.cv$type == 'C-classification' & inner.cv$kernel == 'radial', id])))
#sd(as.numeric(as.character(inner.cv[inner.cv$type == 'C-classification' & inner.cv$kernel == 'sigmoid', id])))
sd(as.numeric(as.character(inner.cv[inner.cv$type == 'C-classification' & inner.cv$kernel == 'polynomial', id])))


#sd(as.numeric(as.character(inner.cv[inner.cv$type == 'nu-classification' & inner.cv$kernel == 'linear', id])))
sd(as.numeric(as.character(inner.cv[inner.cv$type == 'nu-classification' & inner.cv$kernel == 'radial', id])))
#sd(as.numeric(as.character(inner.cv[inner.cv$type == 'nu-classification' & inner.cv$kernel == 'sigmoid', id])))
sd(as.numeric(as.character(inner.cv[inner.cv$type == 'nu-classification' & inner.cv$kernel == 'polynomial', id])))
#print(feature_count)

inner.cv_nu_poly <- inner.cv[inner.cv$type == 'nu-classification' & inner.cv$kernel == 'polynomial', ]
write.csv(inner.cv_nu_poly, './CHEM.CYCLES/234_mRMR_SVM_CHEMCYCLE.csv')





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
#write.csv(count_, './CHEM.CYCLES/selected_features.csv')


# top 30 features
top30Feature <- count_[order(count_$Freq, decreasing = TRUE), ]
top30Feature <- as.character(top30Feature[1:30,1])

