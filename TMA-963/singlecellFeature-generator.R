#################################################
# This script is used to create cell classifier #
#################################################

###################################################################
#### This script is to complete the shape descriptor features #####
###################################################################
library(sf)
library(smoothr)
library(jpeg)
library(gginnards)
library(grid)
library(ggrepel)
library(rjson)
library(umap)

setwd("~/Desktop/TMA_963")

source('./FeatureFunction.R')

HE_detection <- read.csv('./TMA - Coordinates/HE_Coordinates_Features.csv')
# get rid of too large detections
HE_detection <- HE_detection[HE_detection$Area.µm.2 <= 400,]
HE_detection <- HE_detection[complete.cases(HE_detection),]


### define the zoom constant for direction visualiztion
alpha <- 5

# read qualified TMA cores
qualifiedCore <- read.csv('./TMA-963-Annotations/qualified_TMAcores.csv')['x']

# meta file
coreMetaFile <- read.delim('./TMA - Coordinates/TMA_HE_Meta.txt')

# create null matrix
for(coreName in as.character(qualifiedCore$x)){
  
  
  cellFeature_all <- matrix(nrow = 0, ncol = 11)
  orientationList <- matrix(nrow = 0, ncol = 1)
  
  # get core meta file
  #coreName <- 'L-13' # for debug

  # make directory to store feature files
  
  path = paste('./cellFeatures/', coreName, sep = '')
  dir.create(path)
  
  # read metaFile to determine if the core is missing
  category <- as.character(coreMetaFile[coreMetaFile$Name == coreName,]$Missing)
  
   
  if(category == 'False'){
    
    # get the DF subset for current core
    subDF <- HE_detection[HE_detection$TMA.core == coreName,]
    
    
    
    #############################################
    # calculate additional features begins here #
    #############################################
    
    cellShape_list <- read.csv(paste('./TMA - Coordinates/HE_Nucleus_Centroids_combined/detection_coordinate_',coreName, '.csv', sep = ''), row.names = 1)
    

    arrow_vector <- matrix(nrow = 0, ncol = 5)
    
    for(c_id in 1:nrow(cellShape_list)){
      
      if(c_id %in% subDF$Cell.ID){
        # cell boundary coordinate
        cellShape <- data.frame(cellShape_list[cellShape_list$Cell.ID == c_id, 1:2])
        

        cell_princomp <- prcomp(cellShape)
        
        # the orientation is defined in [0, 360]
        orientationAngle <- 180/pi*atan(abs(cell_princomp$rotation[2,1]/cell_princomp$rotation[1,1]))
        if(cell_princomp$rotation[2,1] > 0 & cell_princomp$rotation[1,1] < 0 ){
          orientationAngle <- 180 - orientationAngle
        } else if(cell_princomp$rotation[2,1] < 0 & cell_princomp$rotation[1,1] < 0){
          orientationAngle <- 180 + orientationAngle
        } else if (cell_princomp$rotation[2,1] < 0 & cell_princomp$rotation[1,1] > 0){
          orientationAngle <- 360 - orientationAngle
        }
        orientationList <- rbind(orientationList, orientationAngle)
        
        # save arrow vector
        center_x <- subDF[subDF$Cell.ID == c_id,]$Centroid.X.µm
        center_y <- subDF[subDF$Cell.ID == c_id,]$Centroid.Y.µm
        
        arrow_vector <- rbind(arrow_vector, cbind(c_id, center_x, center_y, 
                                                  center_x + alpha*(cell_princomp$rotation[1,1]), center_y + alpha*cell_princomp$rotation[2,1]))
        

        

    

      } else {
        c_id <- c_id + 1
      }
    }
    
    cellFeature <- data.frame(orientationList)
    colnames(cellFeature) <- c('orientation')
    
    # combine all
    cellFeature_all <- rbind(cellFeature_all, cbind(subDF[,10], subDF[,1:9], cellFeature))
    colnames(cellFeature_all)[1] <- 'cell_ID'
    
    write.csv(cellFeature_all, paste(path, '/', 'cellFeatures.csv', sep = ''))
    
    # write vector files to draw arrows
    colnames(arrow_vector) <- c('Cell.ID', 'x.start','y.start','x.end', 'y.end')
    write.csv(arrow_vector, paste(path, '/', 'arrow_vector.csv', sep = ''))
    
  }
}

