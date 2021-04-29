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

source('~/Desktop/TMA_963/FunctionHub.R')


# read computational features for external validation cohort
{
  CP <- read.csv('./Clinicopathologic Features.MOD.csv')
  CP$CHEMO.CYCLES.MOD <- as.numeric(as.character(CP$CHEMO.CYCLES.MOD))
  CP <- CP[CP$CHEMO.CYCLES.MOD > 2 & CP$TMA,]
  colnames(CP)[1] <- 'patientID'
  #------------- TMA 963 --------------------------------#
  # read patch based feautres
  patch_features963 <- read.csv('~/Desktop/TMA_963/ImageFeatureMatrix-32.csv', row.names = 1)#[,1:1015]
  colnames(patch_features963)[1] <- 'TMA.core'
  
  patch_features963 <- merge(patch_features963, CP, by = 'patientID')[, c(1, 2, 4:1015)]
  patch_features963 <- patch_features963[ , -which(names(patch_features963) %in% c("Density.of.nucleus_CoV",'calc_max_max',"n_protrusion_Max","n_protrusion_Min", 'n_indentation_Max', 'n_indentation_Min', 'Solidity_Max'))]
  
  # read shape features
  shape_features963 <- read.csv('~/Desktop/TMA_963/TMA - Coordinates/shapeFeatures_all.csv', row.names = 1)
  colnames(shape_features963)[1] <- 'TMA.core'
  
  
  # read cluster features
  cluster_features963 <- read.csv('~/Desktop/TMA_963/clusterStat.csv')
  colnames(cluster_features963)[1] <- 'TMA.core'
  
  # read spatstat features
  spatstat_features963 <- read.csv('~/Desktop/TMA_963/SpatialStat.csv', row.names = 1)
  
  
  
  #------------- TMA 1042 --------------------------------#
  
  patch_features1042 <- read.csv('~/Desktop/TMA-coordinates/ImageFeatureMatrix-32.csv', row.names = 1)[,1:1015]
  colnames(patch_features1042)[1] <- 'TMA.core'
  
  
  patch_features1042 <- merge(patch_features1042, CP, by = 'patientID')[, c(1, 2, 4:1015)]
  patch_features1042 <- patch_features1042[ , -which(names(patch_features1042) %in% c("Density.of.nucleus_CoV","calc_max_max","n_protrusion_Max","n_protrusion_Min", 'n_indentation_Max', 'n_indentation_Min', 'Solidity_Max'))]
  
  # read shape features
  shape_features1042 <- read.csv('~/Desktop/TMA-coordinates/shapeFeatures_all.csv', row.names = 1)
  colnames(shape_features1042)[1] <- 'TMA.core'
  
  # read cluster features
  cluster_features1042 <- read.csv('~/Desktop/TMA-coordinates/clusterStat.csv', row.names = 1)
  colnames(cluster_features1042)[1] <- 'TMA.core'
  
  # read spatstat features
  spatstat_features1042 <- read.csv('~/Desktop/TMA - coregistration/SpatialStat.csv', row.names = 1)
  
  
  
  
  
  
  
  #------------- step 3: read pathology feautres ----------#
  
  # add treatment response data to feature file
  Clinicopatho_features<- read.csv('Clinicopathologic Features.csv')[,c(1,15, 19)]
  colnames(Clinicopatho_features) <- c('patientID', 'CD8-FoxP3 ratio','Outcome')
  
  
  
  #----------- step 3: sort spat features -------------#
  spat963 <- colnames(spatstat_features963)
  spat1042 <- colnames(spatstat_features1042)
  
  common_features <- intersect(spat963, spat1042)
  
  # filter uncommon features
  spatstat_features963 <- spatstat_features963[, (colnames(spatstat_features963) %in% common_features)]
  
  spatstat_features1042 <- spatstat_features1042[, (colnames(spatstat_features1042) %in% common_features)]
  
  
  
  
  
  
  
  
  #----------- step 4: combine all features -------------#
  
  
  # read outcome features
  TMA.963.meta <- read.csv('~/Desktop/TMA_963/TMA - Coordinates/TMA-HE-Meta-CellClassifier.csv')
  TMA.1042.meta <- read.csv('~/Desktop/TMA-coordinates/TMA-HE-Meta-CellClassifier.csv')
  
  TMA.963.Outcome <- TMA.963.meta[complete.cases(TMA.963.meta), c(1,4)]
  colnames(TMA.963.Outcome)[1] <- 'TMA.core'
  
  TMA.1042.Outcome <- TMA.1042.meta[TMA.1042.meta$Outcome != '', c(1,4)]
  colnames(TMA.1042.Outcome)[1] <- 'TMA.core'
  
  
  # read classification features
  
  class_features <- read.csv('~/Desktop/TMA_963/classification_features.csv', row.names = 1)
  
  class_features.963 <- class_features[class_features$TMA.set == '963', -2]
  colnames(class_features.963)[1] <- 'TMA.core'
  
  class_features.1042 <- class_features[class_features$TMA.set == '1042',-2]
  colnames(class_features.1042)[1] <- 'TMA.core'
  
  
  # add the valuable clinicopathology features to Feature matrix
  
  
  allFeatures963 <- Reduce(function(x, y) merge(x, y, by = 'TMA.core'), list(patch_features963, shape_features963, cluster_features963,spatstat_features963, class_features.963))
  
  
  
  allFeatures963[is.na(allFeatures963)] <- 0
  #write.csv(allFeatures963, 'allFeatures963_control.csv')
  
  allFeatures963 <- merge(allFeatures963, TMA.963.Outcome, by = 'TMA.core')
  
  
  
  
  allFeatures1042 <- Reduce(function(x, y) merge(x, y, by = 'TMA.core'), list(patch_features1042, shape_features1042, cluster_features1042,spatstat_features1042, class_features.1042))
  
  allFeatures1042[is.na(allFeatures1042)] <- 0
  #write.csv(allFeatures1042, 'allFeatures1042_control.csv')
  
  
  allFeatures1042 <- merge(allFeatures1042, TMA.1042.Outcome, by = 'TMA.core')
  
  #name_vector <- colnames(allFeatures963)
  #colnames(allFeatures1042) <- name_vector
  
  
  
  allFeatures <- rbind(allFeatures963, allFeatures1042)
  
}


 # read survival data
surv <- read.csv('Survival-Data.csv') %>%
  select('Case', 'AGE.AT.OPERATION', 'CD8_FOXP3', 'CLINICAL.T.STAGE')
colnames(surv)[1] <- 'patientID'

# match survival data with core
allFeatures_forModel.withBase <- merge(allFeatures, surv, by = 'patientID')



allFeatures_forModel_temp <- allFeatures_forModel.withBase[, -c(1, 2, 1194)]

# scale features
allFeatures.numbers <- mutate_all(allFeatures_forModel_temp, function(x) as.numeric(as.character(x)))

allFeatures.numbers <- data.frame(scale(allFeatures.numbers))

allFeatures_forModel <- allFeatures.numbers

allFeatures_forModel$Outcome <- as.factor(allFeatures_forModel.withBase$Outcome)

# remove features no less variance
allFeatures_forModel <- allFeatures_forModel[ , -which(names(allFeatures_forModel) %in% c("Density.of.nucleus_CoV",'calc_max_max',"n_protrusion_Max","n_protrusion_Min", 'n_indentation_Max', 'n_indentation_Min', 'Solidity_Max'))]

allFeatures_forModel[is.na(allFeatures_forModel)] <- 0





# For our best model, run Monte Carlo cross validation.
# For each core in the tesing set for each run, record the possibility of being responders

Data.matrix <- allFeatures_forModel#[, -c(1137:1185)]


iterations <- 0



Prob.matrix <- data.frame(matrix(nrow = 142, ncol = 100))

# cross validation using Monte Carlo

boostrap_length = 100 # define boostap length

Partitioned.data <- createDataPartition(Data.matrix$Outcome, times = boostrap_length, p = 0.2, list = T)

inner.cv <- matrix(nrow = 0, ncol = 5)

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
  
  
  mRMR.df <- mRMR.ensemble( data = data, target_indices = ncol(training_fold), solution = 1, feature_count = 40)
  
  mRMR.features <- as.numeric(unlist(mRMR.df@filters))
  
  
  names.vec <- colnames(training_fold[mRMR.features])
  
  
  
  if(length(names.vec) != 0){
    
    feature.selected.mRMR <- training_fold[names.vec]
    feature.selected.mRMR$Outcome <- training_fold$Outcome
    
    # mRMR test_fold
    test_fold.mRMR <- test_fold[names.vec]
    test_fold.mRMR$Outcome <- test_fold$Outcome
    
    #tprfpr_df <- matrix(nrow = 0, ncol = 2)
    
    classifier = svm(formula = Outcome ~ .,
                     data = feature.selected.mRMR,
                     type = 'nu-classification',
                     kernel = 'polynomial',
                     #gamma = 0.001,
                     #cost = 35, 
                     probability = TRUE)      
    # next step in the loop, we calculate the predictions and cm and we equate the accuracy
    y_pred = predict(classifier, newdata = test_fold.mRMR[names.vec],  probability = TRUE) #probability = TRUE
    
    #label.pred <- data.frame(y_pred)
    
    df <- attr(y_pred, 'probabilities')
    
    probs <- df[,'0']
    
    Prob.matrix[row.names(df), id] <-  probs
    
    
    AUC.metric <- AUC(y_pred = y_pred, y_true = test_fold.mRMR$Outcome)
    
    F1.metric <- F1_Score(y_pred = y_pred, y_true = test_fold.mRMR$Outcome, positive = '1')
    
    
    cm = table(pred = y_pred, true=test_fold.mRMR$Outcome)
    accuracy = sum(diag(cm))/sum(cm)
    
    inner.cv <- rbind(inner.cv, cbind('nu-classification', 'polynomial', accuracy, AUC.metric, F1.metric))
    


    #print(inner.cv[which.max(inner.cv[,4]),])
    #print(Accuracy.folds)
    
    iterations <- iterations + 1
    
    print('------------------')
    print(paste('finish iteration ', iterations, sep = ''))
  }
}





# -------------------- Survival Analysis -------------------#



for(thresh in seq(from = 0.1, to = 0.9, by = 0.05)){
  
  thresh <- 0.25
  Predicted_label_all <- data.frame(matrix(nrow = 0, ncol = 0))
  
  for(core in seq_len(142)){
    
    
    Core_data <- as.vector(unlist(Prob.matrix[core, -ncol(Prob.matrix)]))
    
    
    mean_prob <- mean(Core_data, na.rm = TRUE)
    #print(Core_data)
    Predicted_label <- 1*(mean_prob < thresh)
    
    Predicted_label_all <- rbind(Predicted_label_all, cbind(core, Predicted_label))
    
  }
  
  
  surv <- read.csv('Survival-Data.csv')
  colnames(surv)[1] <- 'patientID'
  
  Predicted_label_surv <- merge(allFeatures_forModel.withBase, surv, by = 'patientID') %>%
    select(c('patientID', 'Outcome.y', 'Time', 'Event'))
  
  #Predicted_label_surv
  
  
  Predicted_label_surv$Event <- 1 - Predicted_label_surv$Event
  Predicted_label_surv$Predicted_label <- Predicted_label_all$Predicted_label
  
  
  #crosstab(Predicted_label_surv$Outcome, Predicted_label_surv$Predicted_label)
  require(survival)
  library(survminer)
  
  
  
  fit <- survfit(Surv(Time, Event) ~ Predicted_label, data = Predicted_label_surv)
  
  
  surv_diff <- survdiff(Surv(Time, Event) ~ Predicted_label, data = Predicted_label_surv)
  #print(thresh)
  print(surv_diff)
  
  
  

  #if(summary(res.cox)$coefficients[, 5] < 0.05){
  #  print(thresh)
  #  print(summary(res.cox))
  #}
  

  survplot <- ggsurvplot(fit, data = Predicted_label_surv,
             palette = c("#75c2c9", "#dd7f73"),
             legend.title = "",
             legend.labs = c('Predicted low risk', 'Predicted high risk'),
             risk.table = TRUE, # Add confidence interval
             #risk.table.col = c("#75c2c9", "#dd7f73"),
             pval = TRUE,
             #pval.method = TRUE,
             #pval.size = 7,
             surv.scale="percent",
             font.tickslab = c(16),
             font.title = c(16),
             font.x = c(18),
             font.y = c(18),
             font.legend = c(16),
             ggtheme = theme_bw(),
             fontsize = 6,
             tables.theme = theme(axis.text = element_text(size = 16),
                                  axis.title = element_text(size = 16),
                                  title = element_text(size = 14)),
             risk.table.y.text = FALSE,
             xlim = c(0, 3400)
  ) +
    xlab('Time, (days)')
  survplot
#ggsave(print(survplot), file = paste0('Threshold_', thresh, '_survplot.png'), width = 7, height = 6, units = 'in', dpi = 300)  
}



