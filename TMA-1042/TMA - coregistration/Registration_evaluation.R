#########################################################
#### This script calculate the area for normal tissue####
#########################################################

library(R.matlab)
library(jpeg)
library(ggpubr)
library(gginnards)
library(ggplot2)
library(tidyverse)
library(grid)
library(sp)
library(rgeos)
library(sf)
library(largeVis)
library(concaveman)
library(RSAGA)
library(ggsignif)
# read global and local registration contours

setwd("~/Desktop/TMA - coregistration")

suffix <- '_registeredContour.mat'
label_suffix <- '.mat'
#########################
DSC_all <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(DSC_all) <- c('DSC_score', 'Marker', 'TMA.core', 'Reference')

#case <- 1

for(letter in LETTERS[1:12]){
  for(num in 1:13){
    # core name
    #coreName <- 'G-2'
    coreName <- paste(letter, '-', num, sep = '')
    
    # get the core folder
    FilePath <- paste('./Tissue_Boundary_Registered/', coreName, sep = '')
    if(length(dir(FilePath)) == 10){
      
      #------------------H&E section------------------#
      # set the reference contour
      Reference_contour <- data.frame(readMat(paste(FilePath, '/HE.mat', sep = '')))

      
      # convert to sp format
      sp_ref <- SpatialPolygons(list(Polygons(list(Polygon(as.matrix(Reference_contour[,1:2]))),1)))
      
      # buffer to make valid 
      sp_ref <- gBuffer(sp_ref, byid=TRUE, width=0)
      
      # calculate area
      ref_Area <- gArea(sp_ref) 
      
      for(marker in seq(1,2)){
        marker_name <- switch(marker, 'P53', 'P16')
        target_contour <- data.frame(readMat(paste(FilePath, '/',marker_name, label_suffix, sep = '')))
        
        #plot(target_contour)
        
        ###################################
        
        # convert to sp object for local
        sp_target <- SpatialPolygons(list(Polygons(list(Polygon(as.matrix(target_contour[,1:2]))),1)))
        sp_target <- gBuffer(sp_target, byid=TRUE, width=0)
        
        target_Area <- gArea(sp_target)
        
        
        # calculate Dice coefficietn
        intersection <- gArea(raster::intersect(sp_ref, sp_target))

        # Dice Similarity Coeeficient:
        DSC_score <- 2*intersection/(target_Area + ref_Area)

        DSC_sub <- cbind(DSC_score, marker_name, coreName, 'H&E')
        colnames(DSC_sub) <- c('DSC.score', 'Marker', 'TMA.core', 'Reference')
        DSC_all <- rbind(DSC_all, DSC_sub)
      }
      
      #------------------H&E section------------------#
      # set the reference contour
      Reference_contour <- data.frame(readMat(paste(FilePath, '/Ki67.mat', sep = '')))
      
      
      # convert to sp format
      sp_ref <- SpatialPolygons(list(Polygons(list(Polygon(as.matrix(Reference_contour[,1:2]))),1)))
      
      # buffer to make valid 
      sp_ref <- gBuffer(sp_ref, byid=TRUE, width=0)
      
      # calculate area
      ref_Area <- gArea(sp_ref) 
      
      marker_name <- 'Cyclin'
      target_contour <- data.frame(readMat(paste(FilePath, '/',marker_name, label_suffix, sep = '')))
      
      #plot(target_contour)
      
      ###################################
      
      # convert to sp object for local
      sp_target <- SpatialPolygons(list(Polygons(list(Polygon(as.matrix(target_contour[,1:2]))),1)))
      sp_target <- gBuffer(sp_target, byid=TRUE, width=0)
      
      target_Area <- gArea(sp_target)
      
      
      # calculate Dice coefficietn
      intersection <- gArea(raster::intersect(sp_ref, sp_target))
      
      # Dice Similarity Coeeficient:
      DSC_score <- 2*intersection/(target_Area + ref_Area)
      
      DSC_sub <- cbind(DSC_score, marker_name, coreName, 'Ki67')
      colnames(DSC_sub) <- c('DSC.score', 'Marker', 'TMA.core', 'Reference')
      DSC_all <- rbind(DSC_all, DSC_sub)
      
      #------------------CK5-6 section------------------#
      Reference_contour <- data.frame(readMat(paste(FilePath, '/CK5-6.mat', sep = '')))
      # convert to sp format
      sp_ref <- SpatialPolygons(list(Polygons(list(Polygon(as.matrix(Reference_contour[,1:2]))),1)))
      
      # buffer to make valid 
      sp_ref <- gBuffer(sp_ref, byid=TRUE, width=0)
      
      # calculate area
      ref_Area <- gArea(sp_ref) 
      
      for(marker in seq(1,4)){
        marker_name <- switch(marker, 'Her2Neu', 'P63', 'CK20', 'GATA3')
        target_contour <- data.frame(readMat(paste(FilePath, '/',marker_name, label_suffix, sep = '')))
        
        ###################################
        
        # convert to sp object for local
        sp_target <- SpatialPolygons(list(Polygons(list(Polygon(as.matrix(target_contour[,1:2]))),1)))
        sp_target <- gBuffer(sp_target, byid=TRUE, width=0)
        
        target_Area <- gArea(sp_target)
        
        
        # calculate Dice coefficietn
        intersection <- gArea(raster::intersect(sp_ref, sp_target))
        
        # Dice Similarity Coeeficient:
        DSC_score <- 2*intersection/(target_Area + ref_Area)
        
        DSC_sub <- cbind(DSC_score, marker_name, coreName, 'CK5-6')
        colnames(DSC_sub) <- c('DSC.score', 'Marker', 'TMA.core', 'Reference')
        DSC_all <- rbind(DSC_all, DSC_sub)
      }
      
    }
    
    

  }
}

DSC_all$DSC.score <- as.numeric(as.character(DSC_all$DSC.score))
write.csv(DSC_all, 'DSC_all.csv')

pubmean(DSC_all$DSC.score)

min(DSC_all$DSC.score)

max(DSC_all$DSC.score)

median(DSC_all$DSC.score)


# DSC plot
p <- ggplot(DSC_all, aes(x= TMA.core, y=DSC.score)) + 
  geom_boxplot() +
  theme(axis.text = element_text(angle = 90), axis.text.y = element_text(size = 10), axis.title = element_text(size = 15))

plot(p)
