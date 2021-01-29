###################################################################
#### This script is to complete the shape descriptor features #####
###################################################################
library(sf)
library(smoothr)

setwd("~/Desktop/TMA-coordinates")

source('~/Desktop/TMA_963/FeatureFunction.R')


HE_detection <- read.table('./H&E_allCores_measurements.csv',sep = ',' , header = TRUE)
# get rid of too large detections
#HE_detection <- HE_detection[HE_detection$Area.µm.2 <= 400,]
HE_detection <- HE_detection[complete.cases(HE_detection),]

# create null matrix
shapeFeature_all <- matrix(nrow = 0, ncol = 61)

qualifiedCore <- read.csv('qualified_TMAcores.csv')['x']


# get core meta file
coreMetaFile <- read.delim('TMA_HE_Meta.txt')

for(coreName in as.character(qualifiedCore$x)){
  
  
  category <- as.character(coreMetaFile[coreMetaFile$Name == coreName,]$Missing)
  # 
  if(category == 'False'){
    
    # get the DF subset for current core
    subDF <- HE_detection[HE_detection$TMA.core == coreName,]
    
    # get the inherit stats
    # area
    areaStat <- getBaselineStat(subDF, 'Area.µm.2')
    
    # length
    lengthStat <- getBaselineStat(subDF, 'Length.µm')
    
    # circularity
    circStat <- getBaselineStat(subDF, 'Circularity')
    
    # solidity
    solidStat <- getBaselineStat(subDF, 'Solidity')
    
    # Max diameter
    maxDiameterStat <- getBaselineStat(subDF, 'Max.diameter.µm')
    
    # Min diameter
    minDiameterStat <- getBaselineStat(subDF, 'Min.diameter.µm')
    
    #############################################
    # calculate additional features begins here #
    #############################################
    
    cellShape_list <- fromJSON(file = paste('./TMA_1042_HE_Coords/detection_coordinate_',coreName, '.json', sep = ''))
    
    # loop through all cells
    convexStat <- matrix(nrow = 0, ncol = 2)
    elongation <- matrix(nrow = 0, ncol = 1)
    FD <- matrix(nrow = 0, ncol = 1)
    convexity <- matrix(nrow = 0, ncol = 1)
    curvature <- matrix(nrow = 0, ncol = 6)
    cellMOI <- matrix(nrow = 0, ncol = 1)
    #orientationList <- matrix(nrow = 0, ncol = 1)
    
    
    for(c_id in 1:length(cellShape_list)){
      
      if(c_id %in% subDF$Cell.ID){
        #c_id <- 1
        # cell boundary coordinate
        cellShape <- data.frame(t(data.frame(cellShape_list[[c_id]]$geometry$coordinates)))
        
        # convex hull stats
        chullstat <- getConvexHull(cellShape)
        convexStat <- rbind(convexStat, chullstat)
        
        # elongation
        
        cellBBox <- getMinBBox(cellShape)
        elongation <- rbind(elongation, cellBBox$width/cellBBox$height)
        
        # convexity
        
        perim <- subDF[subDF$Cell.ID == c_id,]$Length.µm # get the perimeter of this cell
        convexity <- rbind(convexity, as.numeric(chullstat[,1]/perim))
        
        # curvature
        smoothedShape <- getSmoothCellShape(cellShape, 200) # smooth the boundary
        curvature <- rbind(curvature, getCurvature(smoothedShape, 200))
        
        # cell orientation
        
        cell_princomp <- prcomp(smoothedShape)
        
        # the orientation is defined in [0, 360]
        #orientationAngle <- 180/pi*atan(abs(cell_princomp$rotation[2,1]/cell_princomp$rotation[1,1]))
        #if(cell_princomp$rotation[2,1] > 0 & cell_princomp$rotation[1,1] < 0 ){
        #  orientationAngle <- 180 - orientationAngle
        #} else if(cell_princomp$rotation[2,1] < 0 & cell_princomp$rotation[1,1] < 0){
        #  orientationAngle <- 180 + orientationAngle
        #} else if (cell_princomp$rotation[2,1] < 0 & cell_princomp$rotation[1,1] > 0){
        #  orientationAngle <- 360 - orientationAngle
        #}
        #orientationList <- rbind(orientationList, orientationAngle)
        
        # fractal dimension
        
        FD <- rbind(FD, getFractalD(smoothedShape))
        
        # moment of inertia
        cellMOI_element <- 0
        cellCentroid <- subDF[subDF$Cell.ID == c_id,]
        translate_cellShape <- sweep(cellShape,2,cbind(cellCentroid$Centroid.X.µm, cellCentroid$Centroid.Y.µm))
        for(i in 1:nrow(translate_cellShape)){
          xi <- translate_cellShape[i,1]
          yi <- translate_cellShape[i,2]
          
          if(i == nrow(translate_cellShape)){
            xi_1 <- translate_cellShape[1,1]
            yi_1 <- translate_cellShape[1,2]
          } else{
            xi_1 <- translate_cellShape[i+1,1]
            yi_1 <- translate_cellShape[i+1,2]
          }
          cellMOI_element <- cellMOI_element + (1/24)*(xi*yi_1 +2*xi*yi+ 2*xi_1*yi_1 + xi_1*yi)*(xi*yi_1 - xi_1*yi)
        }
        cellMOI <- rbind(cellMOI, cellMOI_element)
      } else {
        c_id <- c_id + 1
      }
    }
    
    shapeFeature <- data.frame(cbind(convexStat, elongation, FD, convexity, curvature, cellMOI))
    colnames(shapeFeature) <- c('chullLength', 'chullArea','elongation', 'FD', 'convexity', 'curvMean', 'curvMin', 'curvMax', 'curvStd','n_protrusion', 'n_indentation', 'MOI')
    
    # get COrE
    #colnames(orientationList) <- 'orientationAngle'
    #COrE_features <- getCOrE(orientationList, 18)
    
    # get baseline features
    
    # chull length
    chullLengthStat <- getBaselineStat(shapeFeature, 'chullLength')
    
    # chull area
    chullAreaStat <- getBaselineStat(shapeFeature, 'chullArea')
    
    # elongation
    elongationStat <- getBaselineStat(shapeFeature, 'elongation')
    
    # FD
    FDStat <- getBaselineStat(shapeFeature, 'FD')
    
    # convexity
    convexityStat <- getBaselineStat(shapeFeature, 'convexity')
    
    # curvMean
    curvMeanStat <-mean(shapeFeature$curvMean)
    
    # curvMin
    curvMinStat <-mean(shapeFeature$curvMin)
    
    # curvMax
    curvMaxStat <-mean(shapeFeature$curvMax)
    
    # curvStd
    curvStdStat <-mean(shapeFeature$curvStd)
    
    # n_protrusion
    n_protrusion_Stat <- getBaselineStat(shapeFeature, 'n_protrusion')
    
    # n_indentation
    n_indentation_Stat <- getBaselineStat(shapeFeature, 'n_indentation')
    
    
    # MOI
    MOIStat <- getBaselineStat(shapeFeature, 'MOI')
    
    # combine all
    shapeFeature_all <- rbind(shapeFeature_all, cbind(coreName, areaStat, lengthStat, circStat, 
                                                      solidStat,maxDiameterStat, minDiameterStat,
                                                      chullLengthStat, chullAreaStat, elongationStat,FDStat,
                                                      convexityStat, curvMeanStat, curvMinStat, curvMaxStat,
                                                      curvStdStat, n_protrusion_Stat, n_indentation_Stat, MOIStat)) #COrE_features removed
  }
  print(coreName)
}

row.names(shapeFeature_all) <- NULL
write.csv(shapeFeature_all, './shapeFeatures_all.csv')

