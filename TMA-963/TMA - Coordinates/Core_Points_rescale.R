###########################################################
# This script is used to rescale the points for each core #
###########################################################

library(jpeg)
library(grid)
library(gginnards)
library(stringr)
library(R.matlab)

setwd("~/Desktop/TMA_963/TMA - Coordinates")

################################################
# STEP 1: read the centroid file for each core #
################################################


metaFile <- read.csv('IHC_Core_Stat.csv')  

metaFile$Image <- str_split_fixed(metaFile$Image, '_', n = 6)[,4]


metaFile$Image[metaFile$Image == "uv-GATA3"] <- "GATA3"
metaFile$Image[metaFile$Image == "CyclinD1"] <- "Cyclin"


##########################################
# STEP 2: read the cell coordiantes file #
##########################################

IHC_ptsFile <- read.csv('./IHC_Nucleus_Centroids.csv')
IHC_ptsFile$Image <- str_split_fixed(IHC_ptsFile$Image, '_', n = 6)[,4]

IHC_ptsFile$Image[IHC_ptsFile$Image == "uv-GATA3"] <- "GATA3"
IHC_ptsFile$Image[IHC_ptsFile$Image == "CyclinD1"] <- "Cyclin"



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
      
      IHC <- switch(id, 'Ki67', 'GATA3', 'CK20', 'CK5-6', 'Cyclin', 'P16', 'P53', 'P63', 'Her2Neu')
      
      # subset DF, convert to pixel unit
      ptsFile <- IHC_ptsFile[IHC_ptsFile$TMA.core == coreName & IHC_ptsFile$Image == IHC, c(2,4:5)] 
      
      if(nrow(ptsFile) != 0){
        # get ref points

        ref_x <- metaFile[metaFile$Name == coreName & metaFile$Image == IHC,]$Centroid.X.µm/0.454 - 2908/2 # unit: pixel
        ref_y <- metaFile[metaFile$Name == coreName & metaFile$Image == IHC,]$Centroid.Y.µm/0.454 - 2908/2 # unit: pixel
        
        
        # get the rescaled points
        rescale_x <- ptsFile$Centroid.X.µm/0.454 - ref_x
        rescale_y <- ptsFile$Centroid.Y.µm/0.454 - ref_y  
        
        # get type
        rescale_pts <- cbind(rescale_x, rescale_y) 
          
        # save file to Folder
        writeMat(con = paste("./IHC_Rescaled_Coords/", coreName,'/' ,IHC, '.mat', sep = ''), coords_toReg = rescale_pts)
      }
    }
  }
}   

