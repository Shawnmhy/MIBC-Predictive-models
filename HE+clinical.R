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
library(descr)


setwd("~/Desktop/TMA_model")

#------------------------------------------------------------#
#-------------- Module 1: Data preparation ------------------#
#------------------------------------------------------------#

source('~/Desktop/TMA_963/FunctionHub.R')




#------------- step 1: read features by part ----------#

# read computational features from validation cohort
{
  CP <- read.csv('./Clinicopathologic Features.MOD.csv')
  CP$CHEMO.CYCLES.MOD <- as.numeric(as.character(CP$CHEMO.CYCLES.MOD))
  CP <- CP[CP$CHEMO.CYCLES.MOD > 2,]
  colnames(CP)[1] <- 'patientID'
  
  CP_select <- CP[, c(1,14,15,19, 18)]
  colnames(CP_select)[5] <- 'Outcome'
  
  
  CP_select$AGE.AT.OPERATION <- as.numeric(as.character(CP_select$AGE.AT.OPERATION))
  CP_select$CLINICAL.T.STAGE <- as.numeric(as.character(CP_select$CLINICAL.T.STAGE))
  
  
  CP_select[CP_select$AGE.AT.OPERATION < 60, 3] <- 0
  CP_select[CP_select$AGE.AT.OPERATION >= 60, 3] <- 1
  
  CP_select[CP_select$CLINICAL.T.STAGE <= 2, 4] <- 0
  CP_select[CP_select$CLINICAL.T.STAGE > 2, 4] <- 1
  
  
  
  # read HE features
  
  patch_features963 <- read.csv('~/Desktop/TMA_963/ImageFeatureMatrix-32.csv', row.names = 1)[,1:1015]
  colnames(patch_features963)[1] <- 'TMA.core'
  
  patch_features963 <- merge(patch_features963, CP, by = 'patientID')[, c(1, 2, 4:1015)]
  patch_features963 <- patch_features963[ , -which(names(patch_features963) %in% c("Density.of.nucleus_CoV","calc_max_max","n_protrusion_Max","n_protrusion_Min", 'n_indentation_Max', 'n_indentation_Min', 'Solidity_Max'))]
  
  # read shape features
  shape_features963 <- read.csv('~/Desktop/TMA_963/TMA - Coordinates/shapeFeatures_all.csv', row.names = 1)
  colnames(shape_features963)[1] <- 'TMA.core'
  shape_features963 <- shape_features963[ , -which(names(shape_features963) %in% c("Density.of.nucleus_CoV",'calc_max_max',"n_protrusion_Max","n_protrusion_Min", 'n_indentation_Max', 'n_indentation_Min', 'Solidity_Max'))]
  
  
  # read cluster features
  cluster_features963 <- read.csv('~/Desktop/TMA_963/clusterStat.csv', row.names = 1)
  colnames(cluster_features963)[1] <- 'TMA.core'
  
  
  
  
  #------------- TMA 1042 --------------------------------#
  
  patch_features1042 <- read.csv('~/Desktop/TMA-coordinates/ImageFeatureMatrix-32.csv', row.names = 1)[,1:1015]
  colnames(patch_features1042)[1] <- 'TMA.core'
  patch_features1042 <- merge(patch_features1042, CP, by = 'patientID')[, c(1, 2, 4:1015)]
  patch_features1042 <- patch_features1042[ , -which(names(patch_features1042) %in% c("Density.of.nucleus_CoV","calc_max_max","n_protrusion_Max","n_protrusion_Min", 'n_indentation_Max', 'n_indentation_Min', 'Solidity_Max'))]
  
  
  
  
  # read shape features
  shape_features1042 <- read.csv('~/Desktop/TMA-coordinates/shapeFeatures_all.csv', row.names = 1)
  colnames(shape_features1042)[1] <- 'TMA.core'
  shape_features1042 <- shape_features1042[ , -which(names(shape_features1042) %in% c("Density.of.nucleus_CoV","calc_max_max","n_protrusion_Max","n_protrusion_Min", 'n_indentation_Max', 'n_indentation_Min', 'Solidity_Max'))]
  
  # read cluster features
  cluster_features1042 <- read.csv('~/Desktop/TMA-coordinates/clusterStat.csv', row.names = 1)
  colnames(cluster_features1042)[1] <- 'TMA.core'
  
  
  
  class_features <- read.csv('~/Desktop/TMA_963/classification_features.csv', row.names = 1)
  
  class_features.963 <- class_features[class_features$TMA.set == '963', -2]
  colnames(class_features.963)[1] <- 'TMA.core'
  
  class_features.1042 <- class_features[class_features$TMA.set == '1042',-2]
  colnames(class_features.1042)[1] <- 'TMA.core'
  
  
  
  allFeatures963 <- Reduce(function(x, y) merge(x, y, by = 'TMA.core'), list(patch_features963, shape_features963, cluster_features963, class_features.963))
  
  
  
  allFeatures963[is.na(allFeatures963)] <- 0
  
  
  allFeatures1042 <- Reduce(function(x, y) merge(x, y, by = 'TMA.core'), list(patch_features1042, shape_features1042, cluster_features1042, class_features.1042))
  
  allFeatures1042[is.na(allFeatures1042)] <- 0
  
  
  allFeatures <- rbind(allFeatures963, allFeatures1042)
  
  
  
  
  
  
  #change the column position
  
  allFeatures_forModel <- merge(CP_select, allFeatures, by = 'patientID') %>%
    select(-c('patientID', 'TMA.core'))
  
  
  allFeatures_forModel.temp <- allFeatures_forModel
  
  
  allFeatures_forModel <- allFeatures_forModel[, -4] # indicate which column is Outcome
  
  allFeatures_forModel$Outcome <- allFeatures_forModel.temp$Outcome
  
  
  
  allFeatures_forModel[,1:1171] <- mutate_all(allFeatures_forModel[,1:1171], function(x) as.numeric(as.character(x)))
  
  Data.matrix_forNormal <- allFeatures_forModel
  
  allFeatures_forModel[,1:1171] <- data.frame(scale(allFeatures_forModel[,1:1171]))
  
  
  Data.matrix <- allFeatures_forModel#[,1170:1171]
  Data.matrix[is.na(Data.matrix)] <- 0
  
}



# Strategy I + Demo -> Strategy V
allFeatures_forModel.234_ClinicoPath <- allFeatures_forModel[,-c(1122:1170)]
Data.matrix <- allFeatures_forModel.234_ClinicoPath

optimization.folds.F1 <- matrix(nrow = 0, ncol = 3)
optimization.folds.ACC <- matrix(nrow = 0, ncol = 3)
optimization.folds.AUC <- matrix(nrow = 0, ncol = 3)

iterations <- 0
optimized.feature <- list()



# cross validation using bootstrapping

boostrap_length = 100 # define boostap length

Partitioned.data <- createDataPartition(Data.matrix$Outcome, times = boostrap_length, p = 0.2, list = T)

inner.cv <- matrix(nrow = 0, ncol = 7)


# confusion matrix for each run
cm_addup <- data.frame(matrix(nrow = 2, ncol = 2))
cm_addup[2,2] <- 0
cm_addup[2,1] <- 0
cm_addup[1,2]<- 0
cm_addup[1,1]<- 0
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
  
  
  mRMR.df <- mRMR.ensemble( data = data, target_indices = ncol(training_fold), solution = 1, feature_count = 30)
  
  mRMR.features <- as.numeric(unlist(mRMR.df@filters))
  
  
  names.vec <- colnames(training_fold[mRMR.features])
  
  
  
  
  if(length(names.vec) != 0){
    
    optimized.feature[id] <- list(names.vec)
    feature.selected.mRMR <- training_fold[names.vec]
    feature.selected.mRMR$Outcome <- training_fold$Outcome
    
    # mRMR test_fold
    test_fold.mRMR <- test_fold[names.vec]
    test_fold.mRMR$Outcome <- test_fold$Outcome
    

    
    classifier = svm(formula = Outcome ~ .,
                     data = feature.selected.mRMR,
                     type = 'nu-classification',
                     kernel = 'polynomial',
                     #gamma = 0.001,
                     #cost = 10, 
                     probability = TRUE)      
    # next step in the loop, we calculate the predictions and cm and we equate the accuracy
    y_pred = predict(classifier, newdata = test_fold.mRMR[names.vec],  probability = TRUE) #probability = TRUE
    
    
    
    
    AUC.metric <- AUC(y_pred = y_pred, y_true = test_fold.mRMR$Outcome)
    
    F1.metric <- F1_Score(y_pred = y_pred, y_true = test_fold.mRMR$Outcome, positive = '1')
    
    
    cm = table(pred = y_pred, true=test_fold.mRMR$Outcome)
    
    cm_addup[1,1] <- cm_addup[1,1] + cm[1,1]
    cm_addup[1,2] <- cm_addup[1,2] + cm[1,2]
    cm_addup[2,1] <- cm_addup[2,1] + cm[2,1]
    cm_addup[2,2] <- cm_addup[2,2] + cm[2,2]
    
    accuracy = sum(diag(cm))/sum(cm)
    
    response_rate <- cm[1,1] / (cm[1,1] + cm[1,2])
    
    not_response_rate <- cm[2,1] / (cm[2,1] + cm[2,2])  
    
    #inner.cv <- rbind(inner.cv, cbind(gamma, cost, accuracy, AUC.metric))
    #inner.cv <- rbind(inner.cv, cbind(gamma, cost, accuracy, AUC.metric, F1.metric))
    inner.cv <- rbind(inner.cv, cbind(type, kernel, accuracy, AUC.metric, F1.metric, response_rate,not_response_rate))
    
    # }
    #}
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


# read metrics
inner.cv <- data.frame(inner.cv)

inner.cv[is.nan(as.numeric(as.character(inner.cv$response_rate))),6]<-0

inner.cv$response_rate <- as.numeric(as.character(inner.cv$response_rate))

mean(inner.cv$response_rate, na.rm = TRUE)
sd(inner.cv$response_rate, na.rm = TRUE)







# --------------  Validation --------------#


setwd('~/Desktop/Validation-cohort')
patch_features <- read.csv('./ImageFeatureMatrix-32.csv', row.names = 1)[,1:1015]
colnames(patch_features)[1] <- 'TMA.core'
patch_features$Response <- patch_features$Response - 1
patch_features <- patch_features[ , -which(names(patch_features) %in% c("Density.of.nucleus_CoV",'calc_max_max',"n_protrusion_Max","n_protrusion_Min", 'n_indentation_Max', 'n_indentation_Min', 'Solidity_Max'))]

# read shape features
shape_features <- read.csv('./shapeFeatures_all.csv', row.names = 1)
colnames(shape_features)[1] <- 'TMA.core'
shape_features <- shape_features[ , -which(names(shape_features) %in% c("Density.of.nucleus_CoV",'calc_max_max',"n_protrusion_Max","n_protrusion_Min", 'n_indentation_Max', 'n_indentation_Min', 'Solidity_Max'))]


# read cluster features
cluster_features <- read.csv('./clusterStat.csv', row.names = 1)
colnames(cluster_features)[1] <- 'TMA.core'

#------------- step 3: read pathology feautres ----------#


# read classification features

class_features <- read.csv('./classification_features.csv')
colnames(class_features)[1] <- 'TMA.core'

# add the valuable clinicopathology features to Feature matrix


allFeatures <- Reduce(function(x, y) merge(x, y, by = 'TMA.core'), list(patch_features, shape_features, cluster_features, class_features))




Table_patient_response <- read.csv('Table_patient_response.csv')
Table_patient_response <- Table_patient_response[complete.cases(Table_patient_response),]

Table_patient_core <- read.csv('Table_patient_core.csv')
colnames(Table_patient_core)[2] <- 'patientID'



allFeatures_forModel_valid <- merge(Table_patient_core, allFeatures, by = 'patientID')
colnames(allFeatures_forModel_valid)[1] <- 'patientID'

allFeatures_forModel_valid <- merge(Table_patient_response, allFeatures_forModel_valid, by = 'patient.id')

allFeatures_forModel_valid <- allFeatures_forModel_valid %>%  select(-c('patient.id', 'patientID', 'TMA.core', 'cn', 'ypt', 'ypn', 'response'))




allFeatures_forModel_valid$age <- as.numeric(as.character(allFeatures_forModel_valid$age))
allFeatures_forModel_valid$ct <- as.numeric(as.character(allFeatures_forModel_valid$ct))


# Categorize Clinical T Stage and Age features
allFeatures_forModel_valid[allFeatures_forModel_valid$age < 60, 1] <- 0
allFeatures_forModel_valid[allFeatures_forModel_valid$age >= 60, 1] <- 1

allFeatures_forModel_valid[allFeatures_forModel_valid$ct <= 2, 2] <- 0
allFeatures_forModel_valid[allFeatures_forModel_valid$ct > 2, 2] <- 1

colnames(allFeatures_forModel_valid)[3] <- 'Outcome'


#-------------- validation begin -----------#

allFeatures_forModel.temp <- allFeatures_forModel_valid


# remove outcome
allFeatures_forModel <- allFeatures_forModel_valid[, -3]

allFeatures_forModel$Outcome <- allFeatures_forModel.temp$Outcome



allFeatures_forModel[,1:1170] <- mutate_all(allFeatures_forModel[,1:1170], function(x) as.numeric(as.character(x)))
colnames(allFeatures_forModel)[1] <- 'AGE.AT.OPERATION'
colnames(allFeatures_forModel)[2] <- 'CLINICAL.T.STAGE'



# remove outcome, remove CD8 FoxP3
colmean <- colMeans(Data.matrix_forNormal[,2:1171])

colSTD <- apply(Data.matrix_forNormal[,2:1171], 2, sd)

allFeatures_forModel[,1:1170] <- data.frame(sweep(allFeatures_forModel[,1:1170], 2, colmean)  )
allFeatures_forModel[,1:1170] <- t(apply(allFeatures_forModel[,1:1170], 1, '/', colSTD))





allFeatures_forModel_val_234 <- allFeatures_forModel[, -c(1122:1170)]




iterations <- 0
optimized.feature <- list()



# Training fold: HE + Demo
training_fold = Data.matrix[,-c(1, 1123:1171)] # training fold =  training set minus (-) it's sub test fold

# Test (val) fold: HE + Demo
test_fold = allFeatures_forModel_val_234#[x,]


training_fold$Outcome <- as.numeric(as.character(training_fold$Outcome))
data <- mRMR.data(data.frame(training_fold))


# mRMR feature selection
mRMR.df <- mRMR.ensemble( data = data, target_indices = ncol(training_fold), solution = 1, feature_count = 30)

mRMR.features <- as.numeric(unlist(mRMR.df@filters))

# feature vector
names.vec <- colnames(training_fold[mRMR.features])


# selected feature data
optimized.feature[id] <- list(names.vec)
feature.selected.mRMR <- training_fold[names.vec]
feature.selected.mRMR$Outcome <- training_fold$Outcome

# selected feature test data
test_fold.mRMR <- test_fold[names.vec]
test_fold.mRMR$Outcome <- test_fold$Outcome



classifier = svm(formula = Outcome ~ .,
                 data = feature.selected.mRMR,
                 type = 'nu-classification',
                 kernel = 'polynomial',
                 probability = TRUE)      
# next step in the loop, we calculate the predictions and cm and we equate the accuracy

y_pred = predict(classifier, newdata = test_fold.mRMR[names.vec],  probability = TRUE) #probability = TRUE

probVal <- attr(y_pred, 'probabilities')
predObj <- prediction(probVal[,1], test_fold.mRMR$Outcome)
perf <- performance(predObj, "tpr", "fpr")

plot(perf)

#library(descr)

q1 <- quantile(probVal[,2], 0.25)
q3 <- quantile(probVal[,2], 0.75)

temp <- crosstab((probVal[,2] >= q1) + (probVal[,2] >= q3), test_fold.mRMR$Outcome)
temp



AUC.metric <- AUC(y_pred = y_pred, y_true = test_fold.mRMR$Outcome)


F1.metric <- F1_Score(y_pred = y_pred, y_true = test_fold.mRMR$Outcome, positive = '1')


cm = table(pred = y_pred, true = test_fold.mRMR$Outcome)
cm
accuracy = sum(diag(cm))/sum(cm)
accuracy
#inner.cv <- rbind(inner.cv, cbind(gamma, cost, accuracy, AUC.metric))
inner.cv <- rbind(inner.cv, cbind(type, kernel, accuracy, AUC.metric, F1.metric))


