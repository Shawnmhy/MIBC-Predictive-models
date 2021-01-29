######################################################
# This script is to calculate spatial stats for core #
######################################################

# @author: Haoyang Mi, Johns Hopkins University

library(R.matlab)
library(ggplot2)
library(rjson)
library(igraph)
library(jpeg)
library(grid)
library(gginnards)
library(rjson)
library(geosphere)
library(stringr)
library(pracma)
library(RANN)

setwd("~/Desktop/TMA_963")
source('FunctionHub.R')

#-----This section tells which cells are positive and negative------#

# step 1: read file -> to get classifications
IHC_allPts <- read.csv('./TMA - Coordinates/IHC_Nucleus_Centroids.csv', header = TRUE)

# (clean the dataframe)
IHC_allPts$Class <- as.character(IHC_allPts$Class)

# split the string
IHC_allPts$Image <- str_split_fixed(IHC_allPts$Image, "_", 6)[,4]

# replace CK5-6 to CK56 and uv-GATA3 to GATA3
IHC_allPts$Image <- replace(IHC_allPts$Image, which(IHC_allPts[,1] == 'uv-GATA3'), 'GATA3')

IHC_allPts$Image <- replace(IHC_allPts$Image, which(IHC_allPts[,1] == 'CK5-6'), 'CK56')

# step 2: read qualified position from file
qualifiedCore <- read.csv('./TMA-963-Annotations/qualified_TMAcores.csv')['x']



####################################################################



#################################################
# ---------- Spatial correlation -------------- #
#################################################
allAUC.i2j <- matrix(nrow = 0, ncol = 7)
allAUC.j2i <- matrix(nrow = 0, ncol = 7)
#allCount1 <- matrix(nrow = 0, ncol = 7)
#allCount2 <- matrix(nrow = 0, ncol = 7)
RATIO <- matrix(nrow = 0, ncol = 10)

# allArea
allArea <- read.csv('TMA_allArea.csv', row.names = 1)
colnames(allArea) <- c('TMA.core', 'HE', 'P16', 'P53', 'P63', 'GATA3', 'Ki67', 'CK20', 'CK56', 'Her2Neu', 'CyclinD1', 'Mean', 'CoV', "QCoD")

for(core in as.character(qualifiedCore$x)){
  #core <- 'L-13'
  Ratio <- matrix(nrow = 1, ncol = 1)
  Ratio[1,1] <- core
  # for this core, get classification files for all IHC

  for(IHC_id in 1:9){
    #IHC <- 'Ki67'
    IHC <- switch(IHC_id, 'Ki67', 'GATA3', 'CK20', 'CK56', 'CyclinD1', 'P16', 'P53', 'P63', 'Her2Neu')
    
    # area
    area <- as.numeric(allArea[allArea$TMA.core == core, ][IHC])
    
    # get classification file
    IHC_core_class <- as.character(IHC_allPts[IHC_allPts$Image == IHC & IHC_allPts$TMA.core == core, 2])
    
    # get the registered points, unit in um
    IHC_core_registered <- data.frame(readMat(paste('./TMA - Coordinates/IHC_Registered_Coords/', core, '/', IHC, '.mat', sep = '')))*0.454
    
    # assign 
    IHC_core_registered_class <- cbind(IHC_core_registered, IHC_core_class)
    #
    IHC_core_registered_class <- IHC_core_registered_class[IHC_core_registered_class$IHC_core_class == 'Positive' | IHC_core_registered_class$IHC_core_class == 'Positive: Negative' | IHC_core_registered_class$IHC_core_class == 'Positive: Positive',]
    
    assign(IHC, IHC_core_registered_class)    
    
    # calculate ratio
    #ratio <- nrow(IHC_core_registered_class)/nrow(IHC_core_registered)
    ratio <- nrow(IHC_core_registered_class)/area
    
    Ratio <- cbind(Ratio, ratio)
    
  }
  colnames(Ratio) <- c('TMA.core','Ki67_density', 'GATA3_density', 'CK20_density', 'CK56_density', 'CyclinD1_density', 'P16_density', 'P53_density', 'P63_density', 'Her2Neu_density')
  
  RATIO <- rbind(RATIO, Ratio)
  # HE section is not considered as there is only one associated IHC
  # Ki67 section
  
  Region_Ki67 <- readRDS(paste('./TMA - Coordinates/TMA_SpatStat_Rescaled/', core, '/Ki67.rds', sep = ''))
  Region_Ki67 <- lapply(Region_Ki67,FUN= function(x) x*0.454)
  
  # query: data to be centered 

  
  
  list1 <- bivarAnalysis.Kcross('Ki67', 'P53', Region_Ki67) # P53, CyclinD1
  list2 <- bivarAnalysis.Kcross('Ki67', 'CyclinD1', Region_Ki67) # P53, CyclinD1
  list3 <- bivarAnalysis.Kcross('P53', 'CyclinD1', Region_Ki67) # P53, CyclinD1
  
  # CK5-6 section CK20
  Region_CK56 <- readRDS(paste('./TMA - Coordinates/TMA_SpatStat_Rescaled/', core, '/CK5-6.rds', sep = ''))
  Region_CK56 <- lapply(Region_CK56,FUN= function(x) x*0.454)
  list4 <- bivarAnalysis.Kcross('CK20', 'CK56',Region_CK56)

  # Her2Neu section CK20
  Region_Her2Neu <- readRDS(paste('./TMA - Coordinates/TMA_SpatStat_Rescaled/', core, '/Her2Neu.rds', sep = ''))
  Region_Her2Neu <- lapply(Region_Her2Neu,FUN= function(x) x*0.454)
  
  list5 <- bivarAnalysis.Kcross('Her2Neu', 'GATA3',Region_Her2Neu)
  list6 <- bivarAnalysis.Kcross('P63', 'GATA3',Region_Her2Neu)
  list7 <- bivarAnalysis.Kcross('P63', 'Her2Neu',Region_Her2Neu)
  


  allAUC.i2j <- rbind(allAUC.i2j, cbind(list1[[1]], list2[[1]], list3[[1]], list4[[1]], list5[[1]],
                                  list6[[1]], list7[[1]]))
  # CD: common distance
  allAUC.j2i <- rbind(allAUC.j2i, cbind(list1[[2]], list2[[2]], list3[[2]], list4[[2]], list5[[2]],
                                  list6[[2]], list7[[2]]))
  #allCount1 <- rbind(allCount1, cbind(list1[[3]], list2[[3]], list3[[3]], list4[[3]], list5[[3]],
  #                            list6[[3]], list7[[3]]))
  #allCount2 <- rbind(allCount2, cbind(list1[[4]], list2[[4]], list3[[4]], list4[[4]], list5[[4]],
  #                            list6[[4]], list7[[4]]))
  print(core)
}

colnames(allAUC.i2j) <- c('Ki67_P53.AUC', 'Ki67_Cyclin.AUC', 'P53_Cyclin.AUC', 'CK20_CK56.AUC', 'Her2Neu_GATA3.AUC', 'P63_GATA3.AUC', 'P63_Her2Neu.AUC')
colnames(allAUC.j2i) <- c('P53_Ki67.AUC', 'Cyclin.Ki67.AUC', 'Cyclin_P53.AUC', 'CK56_CK20.AUC', 'GATA3_Her2Neu.AUC', 'GATA3_P63.AUC', 'Her2Neu_P63.AUC')
#colnames(allCount1) <- c('Ki67.Center.P53', 'Ki67.Center.Cyclin', 'P53.Center.Cyclin', 'CK20.Center.CK56', 'Her2.Center.GATA3', 'P63.Center.GATA3', 'P63.Center.Her2')
#colnames(allCount2) <- c('P53.Center.Ki67', 'Cyclin.Center.Ki67', 'Cyclin.Center.P53', 'CK56.Center.CK20', 'GATA3.Center.Her2', 'GATA3.Center.P63', 'Her2.Center.P63')

# summarize data
SpatStat <- data.frame(cbind(RATIO, allAUC.i2j, allAUC.j2i))

# write to file
write.csv(SpatStat, './SpatialStat.csv')

ggplot() +
  geom_polygon(aes(Region_HE[,1], Region_HE[,2]), fill = 'white', color = 'black') +
  geom_point(aes(P53_pts[,1], P53_pts[,2]), fill = 'white', color = 'green') +
  geom_point(aes(P16_pts[,1], P16_pts[,2]), fill = 'white', color = 'red') 





patch_features963 <- read.csv('~/Desktop/TMA_963/TMA - Coordinates/ImageFeatureMatrix.csv', row.names = 1)
colnames(patch_features963)[1] <- 'TMA.core'

allFeatures_patho<- read.csv('~/Desktop/TMA_model/Clinicopathologic Features.csv')[,c(1,19)]
colnames(allFeatures_patho) <- c('patientID', 'Outcome')

Comp_Patho_Features <- merge(allFeatures_patho, patch_features963, by = 'patientID')#[,1:3]


a <- merge(SpatStat, Comp_Patho_Features, by = 'TMA.core')


a[2:10] <- lapply(a[2:10], function(x) as.numeric(as.character(x)))
a$Outcome <- as.factor(a$Outcome)


r <- a[a$Outcome == '0',2]
nr <- a[a$Outcome == '1',2]

t.test(r, nr, var.equal = TRUE)

















