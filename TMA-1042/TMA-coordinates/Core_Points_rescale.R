###########################################################
# This script is used to rescale the points for each core #
###########################################################

library(jpeg)
library(grid)
library(gginnards)
library(stringr)
library(R.matlab)

setwd("~/Desktop/TMA - coregistration")

################################################
# STEP 1: read the centroid file for each core #
################################################

Ki67_metaFile <- read.delim('./IHC_Core_Centroids/Ki67.txt', header = T,sep = '\t')
Ki67_name <- str_split_fixed(Ki67_metaFile$Image, '_', n = 6)
Ki67_metaFile$Image <- Ki67_name[,4]
#----------------------------------#
`uv-GATA3_metaFile` <- read.delim('./IHC_Core_Centroids/GATA3.txt', header = T,sep = '\t')
GATA3_name <- str_split_fixed(`uv-GATA3_metaFile`$Image, '_', n = 6)
`uv-GATA3_metaFile`$Image <- GATA3_name[,4]
#----------------------------------#
Her2Neu_metaFile <- read.delim('./IHC_Core_Centroids/Her2.txt', header = T,sep = '\t')
Her2Neu_name <- str_split_fixed(Her2Neu_metaFile$Image, '_', n = 6)
Her2Neu_metaFile$Image <- Her2Neu_name[,4]
#----------------------------------#
CK20_metaFile <- read.delim('./IHC_Core_Centroids/CK20.txt', header = T,sep = '\t')
CK20_name <- str_split_fixed(CK20_metaFile$Image, '_', n = 6)
CK20_metaFile$Image <- CK20_name[,4]
#----------------------------------#
`CK5-6_metaFile` <- read.delim('./IHC_Core_Centroids/CK56.txt', header = T,sep = '\t')
`CK5-6_name` <- str_split_fixed(`CK5-6_metaFile`$Image, '_', n = 6)
`CK5-6_metaFile`$Image <- `CK5-6_name`[,4]
#----------------------------------#
Cyclin_metaFile <- read.delim('./IHC_Core_Centroids/Cyclin.txt', header = T,sep = '\t')
Cyclin_name <- str_split_fixed(Cyclin_metaFile$Image, '_', n = 6)
Cyclin_metaFile$Image <- Cyclin_name[,4]
#----------------------------------#
P16_metaFile <- read.delim('./IHC_Core_Centroids/P16.txt', header = T,sep = '\t')
P16_name <- str_split_fixed(P16_metaFile$Image, '_', n = 6)
P16_metaFile$Image <- P16_name[,4]
#----------------------------------#
P53_metaFile <- read.delim('./IHC_Core_Centroids/P53.txt', header = T,sep = '\t')
P53_name <- str_split_fixed(P53_metaFile$Image, '_', n = 6)
P53_metaFile$Image <- P53_name[,4]
#----------------------------------#
P63_metaFile <- read.delim('./IHC_Core_Centroids/P63.txt', header = T,sep = '\t')
P63_name <- str_split_fixed(P63_metaFile$Image, '_', n = 6)
P63_metaFile$Image <- P63_name[,4]
#----------------------------------#

##########################################
# STEP 2: read the cell coordiantes file #
##########################################

IHC_ptsFile <- read.csv('IHC_allCores_Centroids.csv')

# split the image name
ImgName <- str_split_fixed(IHC_ptsFile$Image, '_', n = 6)

# assign new name
IHC_ptsFile$Image <- ImgName[,4]


###################################################
# STEP 2: write a for loop and rescale the coords #
###################################################

for(letter in LETTERS[1:12]){
  for(num in 1:10){

    # core name     
    coreName <- paste(letter, '-', num, sep = '')
    
    # create folder
    dir.create(paste("./IHC_Rescaled_Coords/", coreName, sep = ''))
    for(id in 1:9){
      
      IHC <- switch(id, 'Ki67', 'uv-GATA3', 'CK20', 'CK5-6', 'Cyclin', 'P16', 'P53', 'P63', 'Her2Neu')
      # subset DF, convert to pixel unit
      ptsFile <- IHC_ptsFile[IHC_ptsFile$TMA.core == coreName & IHC_ptsFile$Image == IHC, 4:5] 
      
      if(nrow(ptsFile) != 0){
        # get ref points
        metaFile <- eval(as.symbol(paste(IHC, '_metaFile', sep = '')))
        
        ref_x <- metaFile[metaFile$Name == coreName,]$Centroid.X.µm/0.454 - 2908/2 # unit: pixel
        ref_y <- metaFile[metaFile$Name == coreName,]$Centroid.Y.µm/0.454 - 2908/2 # unit: pixel
        
        # get the rescaled points
        rescale_x <- ptsFile$Centroid.X.µm/0.454 - ref_x
        rescale_y <- ptsFile$Centroid.Y.µm/0.454 - ref_y  
        rescale_pts <- cbind(rescale_x, rescale_y)
        # save file to Folder
        writeMat(con = paste("./IHC_Rescaled_Coords/", coreName,'/' ,IHC, '.mat', sep = ''), coords_toReg = rescale_pts)
      }
    }
  }
}   

