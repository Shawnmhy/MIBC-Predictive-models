###########################################################
# This script is used to rescale the points for each core #
###########################################################

library(jpeg)
library(grid)
library(gginnards)
library(stringr)
library(R.matlab)
library(rjson)

setwd("~/Desktop/TMA - coregistration")

################################################
# STEP 1: read the centroid file for each core #
################################################

for(id in 1:10){
  
  Type <- switch(id, 'HE','Ki67', 'GATA3', 'CK20', 'CK5-6', 'Cyclin', 'P16', 'P53', 'P63', 'Her2Neu')
  
  metaFile <- read.delim(paste('./IHC_Core_Centroids/', Type, '.txt', sep = ''), header = T,sep = '\t')
  name <- str_split_fixed(metaFile$Image, '_', n = 6)
  metaFile$Image <- name[,4]

  assign(paste(Type, '_metaFile', sep = ''), metaFile)
}




###################################################
# STEP 2: write a for loop and rescale the coords #
###################################################

for(letter in LETTERS[1:12]){
  for(num in 1:13){
    
    # core name     
    #coreName <- 'A-2'
    coreName <- paste(letter, '-', num, sep = '')
    
    # create folder
    dir.create(paste("./Tissue_Boundary_Rescaled/", coreName, sep = ''))
    for(id in 1:10){
      
      Type <- switch(id, 'HE','Ki67', 'GATA3', 'CK20', 'CK5-6', 'Cyclin', 'P16', 'P53', 'P63', 'Her2Neu')
    
      # folder name
      FilePath <- paste('./TMA_IHC_TissueBoundary/',Type, '/', coreName, sep = '')
      # subset DF, convert to pixel unit

        
      if(length(dir(FilePath)) == 1){
        
        ptsFile <- data.frame(t(data.frame(fromJSON(file = paste(FilePath, '/tissue1.json', sep = ''))$geometry$coordinates[[1]])))
        
        # get ref points
        metaFile <- eval(as.symbol(paste(Type, '_metaFile', sep = '')))
        
        ref_x <- metaFile[metaFile$Name == coreName,]$Centroid.X.um/0.454 - 3030/2 # unit: pixel
        ref_y <- metaFile[metaFile$Name == coreName,]$Centroid.Y.um/0.454 - 3030/2 # unit: pixel
        
        # get the rescaled points
        rescale_x <- ptsFile$X1 - ref_x
        rescale_y <- ptsFile$X2 - ref_y  
        rescale_pts <- cbind(rescale_x, rescale_y)
        # save file to Folder
        writeMat(con = paste("./Tissue_Boundary_Rescaled/", coreName,'/' ,Type, '.mat', sep = ''), coords_toReg = rescale_pts)
      }
    }
  }
}   


###################################################
# STEP 2: write a for loop and rescale the coords #
###################################################

for(letter in LETTERS[1:12]){
  for(num in 1:13){
    
    # core name     
    #coreName <- 'A-2'
    coreName <- paste(letter, '-', num, sep = '')
    
    # create folder
    dir.create(paste("./TMA_SpatStat_Rescaled/", coreName, sep = ''))
    for(id in 1:3){
      
      Type <- switch(id, 'HE','Ki67', 'CK5-6')
      
      # folder name
      FilePath <- paste('./Reference_Countour/',Type, '/', coreName, sep = '')
      # subset DF, convert to pixel unit
      
      Region_list <- list()
      
      if(length(dir(FilePath) != 0)){
        for(n in 1:length(dir(FilePath))){
          
          ptsFile <- data.frame(t(data.frame(fromJSON(file = paste(FilePath, '/tissue', n, '.json', sep = ''))$geometry$coordinates[[1]])))
          
          # get ref points
          metaFile <- eval(as.symbol(paste(Type, '_metaFile', sep = '')))
          
          ref_x <- metaFile[metaFile$Name == coreName,]$Centroid.X.um - 0.454*3030/2 # unit: pixel
          ref_y <- metaFile[metaFile$Name == coreName,]$Centroid.Y.um - 0.454*3030/2 # unit: pixel
          
          # get the rescaled points
          rescale_x <- ptsFile$X1*0.454 - ref_x
          rescale_y <- ptsFile$X2*0.454 - ref_y  
          Region_list[[n]] <- cbind(rev(rescale_x), rev(rescale_y))
          # save file to Folder
        }
        saveRDS(Region_list, file = paste("./TMA_SpatStat_Rescaled/", coreName,'/' ,Type, '.rds', sep = ''))
      }
    }
  }
}   

