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

setwd("~/Desktop/TMA - coregistration")
source('FunctionHub.R')

#-----This section tells which cells are positive and negative------#

# step 1: read file
IHC_allPts <- read.csv('IHC_forDensityCount.csv', header = TRUE)[,-4]

# (clean the dataframe)
IHC_allPts$Class <- as.character(IHC_allPts$Class)
IHC_allPts <- IHC_allPts[IHC_allPts$Class != '',]

# split the string
IHC_allPts$Image <- str_split_fixed(IHC_allPts$Image, "_", 6)[,4]

# replace CK5-6 to CK56 and uv-GATA3 to GATA3
IHC_allPts$Image <- replace(IHC_allPts$Image, which(IHC_allPts[,1] == 'uv-GATA3'), 'GATA3')
IHC_allPts$Image <- replace(IHC_allPts$Image, which(IHC_allPts[,1] == 'CK5-6'), 'CK56')



# step 2: read qualified position from file
qualifiedCore <- read.csv('qualified_TMAcores.csv')['x']

# step 3: use above dat to filter dat in step 1
qIHC_allPts <- IHC_allPts[IHC_allPts$TMA.core %in% qualifiedCore$x,]

# step 4: get the correct class

Negative_id <- sapply('Negative', function(y) grep(y,qIHC_allPts$Class))

# get a all id list
all_id <- as.vector(1:nrow(qIHC_allPts))

Positive_id <- setdiff(all_id, Negative_id)

# reassign classes
qIHC_allPts$Class <- as.character(qIHC_allPts$Class)

qIHC_allPts$Class[Positive_id] <- 'Positive'

qIHC_allPts$Class[Negative_id] <- 'Negative'
####################################################################

Ratio <- matrix(nrow = 0, ncol = 3)

for(core in as.character(qualifiedCore$x)){

  # read tissue boundary file

  #core <- 'A-1'
  #--------- SpatStat Part I: cell ratio------------#

  for(IHC_id in 1:9){
    # read IHC pos files in a for loop
    IHC <- switch(IHC_id, 'Ki67', 'Cyclin', 'P16', 'P53', 'CK56', 'GATA3', 'CK20', 'Her2Neu', 'P63')
    
    # IHC coordinates
    posFile <- data.frame(readMat(paste('./IHC_Registered_Coords/', core, '/',IHC, '.mat', sep = '')))
    
    assign(IHC, posFile)
    
    # IHC ratio
    qIHC_this <- qIHC_allPts[qIHC_allPts$Image == IHC & qIHC_allPts$TMA.core == core, ]
    ratio <- nrow(qIHC_this[qIHC_this$Class == 'Positive', ])/nrow(qIHC_this)
    
    Ratio <- data.frame(rbind(Ratio, cbind(core, IHC, ratio)))
    
  }
  

}

Ratio$core <- as.character(Ratio$core)
Ratio$IHC <- as.character(Ratio$IHC)
Ratio$ratio <- as.numeric(as.character(Ratio$ratio))

# unmelt data
unmelt_Ratio <- dcast(data = Ratio, formula = core~IHC, fun.aggregate = sum, value.var = "ratio")
write.csv(unmelt_Ratio, './ratioIHC.csv')

#################################################
# ---------- Spatial correlation -------------- #
#################################################
allAUC <- matrix(nrow = 0, ncol = 12)
allCD <- matrix(nrow = 0, ncol = 12)
for(core in as.character(qualifiedCore$x)){
  
  # H&E section
  #core <- 'D-9'
  Region_HE <- readRDS(paste('./TMA_SpatStat_Rescaled/', core, '/HE.rds', sep = ''))
  
  # read pts dat

  
  list1 <- bivarAnalysis('P53', 'P16', Region_HE)

  # Ki67 section
  Region_Ki67 <- readRDS(paste('./TMA_SpatStat_Rescaled/', core, '/Ki67.rds', sep = ''))
  list2 <- bivarAnalysis('Ki67', 'Cyclin', Region_Ki67)

  
  # CK5-6 section
  Region_CK56 <- readRDS(paste('./TMA_SpatStat_Rescaled/', core, '/CK5-6.rds', sep = ''))
  list3 <- bivarAnalysis('GATA3', 'CK56',Region_CK56)
  list4 <- bivarAnalysis('Her2Neu', 'CK56',Region_CK56)
  list5 <- bivarAnalysis('CK20', 'CK56',Region_CK56)
  list6 <- bivarAnalysis('P63', 'CK56',Region_CK56)
  
  list7 <- bivarAnalysis('Her2Neu', 'GATA3',Region_CK56)
  list8 <- bivarAnalysis('CK20', 'GATA3',Region_CK56)
  list9 <- bivarAnalysis('P63', 'GATA3',Region_CK56)
  
  list10 <- bivarAnalysis('Her2Neu', 'CK20',Region_CK56)
  list11 <- bivarAnalysis('P63', 'CK20',Region_CK56)
  
  list12 <- bivarAnalysis('Her2Neu', 'P63',Region_CK56)
  
  allAUC <- rbind(allAUC, cbind(list1[[1]], list2[[1]], list3[[1]], list4[[1]], list5[[1]],
                                  list6[[1]], list7[[1]], list8[[1]], list9[[1]], list10[[1]], list11[[1]],list12[[1]]))
  # CD: common distance
  allCD <- rbind(allCD, cbind(list1[[2]], list2[[2]], list3[[2]], list4[[2]], list5[[2]],
                                  list6[[2]], list7[[2]], list8[[2]], list9[[2]], list10[[2]], list11[[2]],list12[[2]]))
  print(core)
}

colnames(allCD) <- c('P53_P16.CD', 'Ki67_Cyclin.CD', 'GATA3_CK56.CD', 'Her2Neu_CK56.CD', 'CK20_CK56.CD', 'P63_CK56.CD', 'Her2Neu_GATA3.CD',
                     'CK20_GATA3.CD', 'P63_GATA3.CD', 'Her2Neu_CK20.CD', 'P63_CK20.CD', 'Her2Neu_P63.CD')
colnames(allAUC) <-  c('P53_P16.AUC', 'Ki67_Cyclin.AUC', 'GATA3_CK56.AUC', 'Her2Neu_CK56.AUC', 'CK20_CK56.AUC', 'P63_CK56.AUC', 'Her2Neu_GATA3.AUC',
                       'CK20_GATA3.AUC', 'P63_GATA3.AUC', 'Her2Neu_CK20.AUC', 'P63_CK20.AUC', 'Her2Neu_P63.AUC')
# write to file
write.csv(allAUC, './allAUC.csv')
write.csv(allCD, './allCD.csv')

ggplot() +
  geom_polygon(aes(Region_HE[,1], Region_HE[,2]), fill = 'white', color = 'black') +
  geom_point(aes(P53_pts[,1], P53_pts[,2]), fill = 'white', color = 'green') +
  geom_point(aes(P16_pts[,1], P16_pts[,2]), fill = 'white', color = 'red') 























