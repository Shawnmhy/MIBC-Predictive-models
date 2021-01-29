#####################################################
#### This script is to calculate feature matrix #####
#####################################################

library(largeVis)
library(rjson)
library(stats)
library(fractaldim)
library(ggplot2)
# polygon fit
library(sf)
library(smoothr)
# image texture
library(imager)
library(radiomics)
#calculate boundary features
library(pracma)
# redirect to working directory
setwd("~/Desktop/TMA_963")

# load function file
source('./FeatureFunction.R')
#####################################################
#-------------- define feature matrix --------------#
#####################################################


#! all unit in µm, not pixel!
######################################
#--------- TMA basic stats ----------#
######################################
# TMA core names
coreStats <- matrix(nrow = 0, ncol = 3) # names, patient ID, tissue types

######################################################
#--------- General: image texture features ----------#
######################################################

textureFeature_core <- data.frame(matrix(0, nrow = 1, ncol = 265*4))
textureFeature_total <- matrix(nrow = 0, ncol = 265*4)


# Total features count: 35


# Total features count: 37


#####################################################
#-------------- Features calculation- --------------#
#####################################################
coreStats <- matrix(nrow = 0, ncol = 3) # names, patient ID, tissue types

# read qualified cores

qualifiedCore <- read.csv('~/Desktop/TMA_963/TMA-963-Annotations/qualified_TMAcores.csv')['x']

## meta file

coreMetaFile <- read.delim('~/Desktop/TMA_963/TMA - Coordinates/TMA_HE_Meta.txt')

# calculate features for each core

for(coreName in as.character(qualifiedCore$x)){
  print(coreName)
  #coreName <- 'A-2'
  ## tiles
  tilePath <- paste('./TMA_HE_tiles/', coreName, sep ='')
  ## coordiante file
  category <- as.character(coreMetaFile[coreMetaFile$Name == coreName,]$Missing)
  
  # check availability
  if(category == 'False'){
    # assign meta stats
    patientID <- as.character(coreMetaFile[coreMetaFile$Name == coreName,]$Unique.ID)
    tissueType <- as.character(coreMetaFile[coreMetaFile$Name == coreName,]$Tissue.type)
    
    coreStats <- rbind(coreStats, cbind(coreName, patientID, tissueType))
    
    ###################################################################
    #------------------ image texutre feature begins------------------#
    ###################################################################
    
    textureFeature <- matrix(nrow = 0, ncol = 13)
    textureFeature2 <- data.frame(matrix(0, nrow = length(list.files(tilePath)), ncol = 21*12))
    for(f in 1:length(list.files(tilePath))){
      # read tile RGB image and convert to greyscale image
      tileName <- paste(tilePath, '/tiles_', f, '.jpg', sep = '')
      tileImg <- as.matrix(grayscale(load.image(tileName), method = "Luma", drop = TRUE))
      
      tileImg <- floor(tileImg*32)
      
      #plot(tileImg.plot)
      #tileImg.floor <- floor(tileImg*8)
      tst <- load.image(tileName)
      
      # 1st order stats
      textureFeature <- rbind(textureFeature, cbind(calc_features(tileImg)))
      # 2nd order stats
      count <- 0
      featureAffix <- 'glcm'
      for(angle in c(0,45,90,135)){
        for(d in c(1,3,5)){
          count <- count + 1
          hbGLCM <- glcm(tileImg, angle = angle, d= d, n_grey= 32)
          textureFeature2[f,(21*count - 20):(21*count)] <- calc_features(hbGLCM)
          # assign feature names
          names(textureFeature2)[(21*count - 20):(21*count)] <-
            c(paste(featureAffix, '_mean_',angle,'_',d,sep = ''),
              paste(featureAffix, '_variance_',angle,'_',d,sep = ''),
              paste(featureAffix, '_autoCorrelation_',angle,'_',d,sep = ''),
              paste(featureAffix, '_cProminence_',angle,'_',d,sep = ''),
              paste(featureAffix, '_cShade_',angle,'_',d,sep = ''),
              paste(featureAffix, '_cTendency_',angle,'_',d,sep = ''),
              paste(featureAffix, '_contrast_',angle,'_',d,sep = ''),
              paste(featureAffix, '_correlation_',angle,'_',d,sep = ''),
              paste(featureAffix, '_differentEntropy_',angle,'_',d,sep = ''),
              paste(featureAffix, '_dissimilarity_',angle,'_',d,sep = ''),
              paste(featureAffix, '_energy_',angle,'_',d,sep = ''),
              paste(featureAffix, '_entropy_',angle,'_',d,sep = ''),
              paste(featureAffix, '_homogeneity1_',angle,'_',d,sep = ''),
              paste(featureAffix, '_homogeneity2_',angle,'_',d,sep = ''),
              paste(featureAffix, '_IDMN_',angle,'_',d,sep = ''),
              paste(featureAffix, '_IDN_',angle,'_',d,sep = ''),
              paste(featureAffix, '_inverseVariance_',angle,'_',d,sep = ''),
              paste(featureAffix, '_maxProb_',angle,'_',d,sep = ''),
              paste(featureAffix, '_sumAverage_',angle,'_',d,sep = ''),
              paste(featureAffix, '_sumEntropy_',angle,'_',d,sep = ''),
              paste(featureAffix, '_sumVariance_',angle,'_',d,sep = '')
            )
        }
      }
      # calculate features
    }
    # combine two sets of features
    Feature_Total <- cbind(data.frame(textureFeature), textureFeature2)
    for(flag in 1:253){ # 13 1st order feature and 252 2nd order feature
      # declare a temporary matrix
      # define feature name
      metric <- switch(flag, 'calc_energy', 'calc_entropy', 'calc_kurtosis', 'calc_meanDeviation', 'calc_skewness', 'calc_uniformity',
                       'calc_mean', 'calc_median', 'calc_max', 'calc_min', 'calc_variance', 'calc_RMS', 'calc_sd', 
                       'glcm_mean_0_1', 'glcm_mean_0_3', 'glcm_mean_0_5',
                       'glcm_mean_45_1', 'glcm_mean_45_3', 'glcm_mean_45_5',
                       'glcm_mean_90_1', 'glcm_mean_90_3', 'glcm_mean_90_5',
                       'glcm_mean_135_1', 'glcm_mean_135_3', 'glcm_mean_135_5',
                       #---------------------------------------------------#
                       'glcm_variance_0_1', 'glcm_variance_0_3', 'glcm_variance_0_5',
                       'glcm_variance_45_1', 'glcm_variance_45_3', 'glcm_variance_45_5',
                       'glcm_variance_90_1', 'glcm_variance_90_3', 'glcm_variance_90_5',
                       'glcm_variance_135_1', 'glcm_variance_135_3', 'glcm_variance_135_5',
                       #---------------------------------------------------#
                       'glcm_autoCorrelation_0_1', 'glcm_autoCorrelation_0_3', 'glcm_autoCorrelation_0_5',
                       'glcm_autoCorrelation_45_1', 'glcm_autoCorrelation_45_3', 'glcm_autoCorrelation_45_5',
                       'glcm_autoCorrelation_90_1', 'glcm_autoCorrelation_90_3', 'glcm_autoCorrelation_90_5',
                       'glcm_autoCorrelation_135_1', 'glcm_autoCorrelation_135_3', 'glcm_autoCorrelation_135_5',
                       #---------------------------------------------------#
                       'glcm_cProminence_0_1', 'glcm_cProminence_0_3', 'glcm_cProminence_0_5',
                       'glcm_cProminence_45_1', 'glcm_cProminence_45_3', 'glcm_cProminence_45_5',
                       'glcm_cProminence_90_1', 'glcm_cProminence_90_3', 'glcm_cProminence_90_5',
                       'glcm_cProminence_135_1', 'glcm_cProminence_135_3', 'glcm_cProminence_135_5',
                       #---------------------------------------------------#
                       'glcm_cShade_0_1', 'glcm_cShade_0_3', 'glcm_cShade_0_5',
                       'glcm_cShade_45_1', 'glcm_cShade_45_3', 'glcm_cShade_45_5',
                       'glcm_cShade_90_1', 'glcm_cShade_90_3', 'glcm_cShade_90_5',
                       'glcm_cShade_135_1', 'glcm_cShade_135_3', 'glcm_cShade_135_5',
                       #---------------------------------------------------#
                       'glcm_cTendency_0_1', 'glcm_cTendency_0_3', 'glcm_cTendency_0_5',
                       'glcm_cTendency_45_1', 'glcm_cTendency_45_3', 'glcm_cTendency_45_5',
                       'glcm_cTendency_90_1', 'glcm_cTendency_90_3', 'glcm_cTendency_90_5',
                       'glcm_cTendency_135_1', 'glcm_cTendency_135_3', 'glcm_cTendency_135_5',
                       #---------------------------------------------------#
                       'glcm_contrast_0_1', 'glcm_contrast_0_3', 'glcm_contrast_0_5',
                       'glcm_contrast_45_1', 'glcm_contrast_45_3', 'glcm_contrast_45_5',
                       'glcm_contrast_90_1', 'glcm_contrast_90_3', 'glcm_contrast_90_5',
                       'glcm_contrast_135_1', 'glcm_contrast_135_3', 'glcm_contrast_135_5',
                       #---------------------------------------------------#
                       'glcm_correlation_0_1', 'glcm_correlation_0_3', 'glcm_correlation_0_5',
                       'glcm_correlation_45_1', 'glcm_correlation_45_3', 'glcm_correlation_45_5',
                       'glcm_correlation_90_1', 'glcm_correlation_90_3', 'glcm_correlation_90_5',
                       'glcm_correlation_135_1', 'glcm_correlation_135_3', 'glcm_correlation_135_5',
                       #---------------------------------------------------#
                       'glcm_dissimilarity_0_1', 'glcm_dissimilarity_0_3', 'glcm_dissimilarity_0_5',
                       'glcm_dissimilarity_45_1', 'glcm_dissimilarity_45_3', 'glcm_dissimilarity_45_5',
                       'glcm_dissimilarity_90_1', 'glcm_dissimilarity_90_3', 'glcm_dissimilarity_90_5',
                       'glcm_dissimilarity_135_1', 'glcm_dissimilarity_135_3', 'glcm_dissimilarity_135_5',
                       #---------------------------------------------------#
                       'glcm_energy_0_1', 'glcm_energy_0_3', 'glcm_energy_0_5',
                       'glcm_energy_45_1', 'glcm_energy_45_3', 'glcm_energy_45_5',
                       'glcm_energy_90_1', 'glcm_energy_90_3', 'glcm_energy_90_5',
                       'glcm_energy_135_1', 'glcm_energy_135_3', 'glcm_energy_135_5',
                       #---------------------------------------------------#
                       'glcm_entropy_0_1', 'glcm_entropy_0_3', 'glcm_entropy_0_5',
                       'glcm_entropy_45_1', 'glcm_entropy_45_3', 'glcm_entropy_45_5',
                       'glcm_entropy_90_1', 'glcm_entropy_90_3', 'glcm_entropy_90_5',
                       'glcm_entropy_135_1', 'glcm_entropy_135_3', 'glcm_entropy_135_5',
                       #---------------------------------------------------#
                       'glcm_homogeneity1_0_1', 'glcm_homogeneity1_0_3', 'glcm_homogeneity1_0_5',
                       'glcm_homogeneity1_45_1', 'glcm_homogeneity1_45_3', 'glcm_homogeneity1_45_5',
                       'glcm_homogeneity1_90_1', 'glcm_homogeneity1_90_3', 'glcm_homogeneity1_90_5',
                       'glcm_homogeneity1_135_1', 'glcm_homogeneity1_135_3', 'glcm_homogeneity1_135_5',
                       #---------------------------------------------------#
                       'glcm_homogeneity2_0_1', 'glcm_homogeneity2_0_3', 'glcm_homogeneity2_0_5',
                       'glcm_homogeneity2_45_1', 'glcm_homogeneity2_45_3', 'glcm_homogeneity2_45_5',
                       'glcm_homogeneity2_90_1', 'glcm_homogeneity2_90_3', 'glcm_homogeneity2_90_5',
                       'glcm_homogeneity2_135_1', 'glcm_homogeneity2_135_3', 'glcm_homogeneity2_135_5',
                       #---------------------------------------------------#
                       'glcm_IDMN_0_1', 'glcm_IDMN_0_3', 'glcm_IDMN_0_5',
                       'glcm_IDMN_45_1', 'glcm_IDMN_45_3', 'glcm_IDMN_45_5',
                       'glcm_IDMN_90_1', 'glcm_IDMN_90_3', 'glcm_IDMN_90_5',
                       'glcm_IDMN_135_1', 'glcm_IDMN_135_3', 'glcm_IDMN_135_5',
                       #---------------------------------------------------#
                       'glcm_IDN_0_1', 'glcm_IDN_0_3', 'glcm_IDN_0_5',
                       'glcm_IDN_45_1', 'glcm_IDN_45_3', 'glcm_IDN_45_5',
                       'glcm_IDN_90_1', 'glcm_IDN_90_3', 'glcm_IDN_90_5',
                       'glcm_IDN_135_1', 'glcm_IDN_135_3', 'glcm_IDN_135_5',
                       #---------------------------------------------------#
                       'glcm_inverseVariance_0_1', 'glcm_inverseVariance_0_3', 'glcm_inverseVariance_0_5',
                       'glcm_inverseVariance_45_1', 'glcm_inverseVariance_45_3', 'glcm_inverseVariance_45_5',
                       'glcm_inverseVariance_90_1', 'glcm_inverseVariance_90_3', 'glcm_inverseVariance_90_5',
                       'glcm_inverseVariance_135_1', 'glcm_inverseVariance_135_3', 'glcm_inverseVariance_135_5',
                       #---------------------------------------------------#
                       'glcm_maxProb_0_1', 'glcm_maxProb_0_3', 'glcm_maxProb_0_5',
                       'glcm_maxProb_45_1', 'glcm_maxProb_45_3', 'glcm_maxProb_45_5',
                       'glcm_maxProb_90_1', 'glcm_maxProb_90_3', 'glcm_maxProb_90_5',
                       'glcm_maxProb_135_1', 'glcm_maxProb_135_3', 'glcm_maxProb_135_5',
                       #---------------------------------------------------#
                       'glcm_sumAverage_0_1', 'glcm_sumAverage_0_3', 'glcm_sumAverage_0_5',
                       'glcm_sumAverage_45_1', 'glcm_sumAverage_45_3', 'glcm_sumAverage_45_5',
                       'glcm_sumAverage_90_1', 'glcm_sumAverage_90_3', 'glcm_sumAverage_90_5',
                       'glcm_sumAverage_135_1', 'glcm_sumAverage_135_3', 'glcm_sumAverage_135_5',
                       #---------------------------------------------------#
                       'glcm_sumEntropy_0_1', 'glcm_sumEntropy_0_3', 'glcm_sumEntropy_0_5',
                       'glcm_sumEntropy_45_1', 'glcm_sumEntropy_45_3', 'glcm_sumEntropy_45_5',
                       'glcm_sumEntropy_90_1', 'glcm_sumEntropy_90_3', 'glcm_sumEntropy_90_5',
                       'glcm_sumEntropy_135_1', 'glcm_sumEntropy_135_3', 'glcm_sumEntropy_135_5',
                       #---------------------------------------------------#
                       'glcm_sumVariance_0_1', 'glcm_sumVariance_0_3', 'glcm_sumVariance_0_5',
                       'glcm_sumVariance_45_1', 'glcm_sumVariance_45_3', 'glcm_sumVariance_45_5',
                       'glcm_sumVariance_90_1', 'glcm_sumVariance_90_3', 'glcm_sumVariance_90_5',
                       'glcm_sumVariance_135_1', 'glcm_sumVariance_135_3', 'glcm_sumVariance_135_5'
      )
      
      # calcualte statistics for 1st order features
      
      maxFeature <- max(Feature_Total[[metric]])
      minFeature <- min(Feature_Total[[metric]])
      meanFeature <- mean(Feature_Total[[metric]])
      stdFeature <- sd(Feature_Total[[metric]])
      
      
      # combine image texture features
      # note: here a trick to calculate the general term formular of arithmetic sequence is used
      textureFeature_core[1, (4*flag -3):(4*flag)] <- data.frame(cbind(maxFeature, minFeature, meanFeature, stdFeature))
      names(textureFeature_core)[(4*flag -3):(4*flag)] <- 
        c(paste(metric, '_max', sep =''),
          paste(metric, '_min', sep =''),
          paste(metric, '_mean', sep =''),
          paste(metric, '_std', sep =''))
    }
    textureFeature_total <- rbind(textureFeature_total, textureFeature_core)
    FeatureMatrix <- cbind(coreStats, textureFeature_total)
    
    #################################################################
    #------------------ image texutre feature ends------------------#
    #################################################################
    
  }
}
write.csv(FeatureMatrix, './ImageFeatureMatrix-32.csv')



FeatureMatrix.ori <- read.csv('./ImageFeatureMatrix.csv')



