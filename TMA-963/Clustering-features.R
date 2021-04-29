######################################################
# This script is to calculate spatial stats for core #
######################################################

# @author: Haoyang Mi, Johns Hopkins University

library(R.matlab)
library(ggplot2)
library(rjson)
library(jpeg)
library(grid)
library(gginnards)
library(spatstat)
library(tripack)
library(dplyr)
library(rgeos)
library(igraph)
library(stringr)
library(concaveman)

setwd("~/Desktop/TMA-coordinates")
source('~/Desktop/TMA_963/FeatureFunction.R')
###################################################
# STEP 1: rescale H&E points #
###################################################

#-------- read meta file to get the core boundary -------------#


# read area file

TMA_allArea <- read.csv('./TMA_allArea.csv', row.names = 1)[,1:2]

# read qualified cores
qualifiedCore <- read.csv('~/Desktop/TMA-1042-Annotations/qualified_TMAcores.csv')['x']
coreClusStat <- matrix(nrow = 0, ncol = 15)

# read reference contour
for(core in as.character(qualifiedCore$x)){
  
  #core <- 'L-13'
  Coords_shapeDes <- read.csv(paste('./cellFeatures/', core,'/cellFeatures.csv', sep = ''), row.names = 1)[,c(1,3,4,11)]
  Coords_shapeDes$cell_ID <- c(1:nrow(Coords_shapeDes))
  # reference
  Region_HE <- readRDS(paste('~/Desktop/TMA - coregistration//TMA_SpatStat_Rescaled/', core, '/HE.rds', sep = ''))
  Region_HE <- lapply(Region_HE,FUN= function(x) x*0.454)
  
  # pts to ppp
  pts <- data.frame(readMat(paste('./HE_Rescaled_Coords/', core, '/HE.mat', sep = '')))*0.454
  colnames(pts) <- c('x', 'y')
  
  
  
  # replace the coordinates to rescaled coordinates
  Coords_shapeDes$Centroid.X.µm <- pts$x
  Coords_shapeDes$Centroid.Y.µm <- pts$y
  
  #
  #---------------------------------------#
  # clustering using delaunay triangulation
  
  r <- tri.mesh(pts$x, pts$y)
  
  
  # get the delaunay triangulation basic stats
  delaunayTri <- tripack::triangles(r)
  
  Num_tri <- nrow(delaunayTri)  
  # coord length of triangle 
  a <- delaunayTri[,1] # node 1 index
  b <- delaunayTri[,2] # node 2 index
  c <- delaunayTri[,3] # node 3 index
  tri_sidelength <- cbind(sqrt( ( r$x[a] - r$x[b] ) ^ 2 + ( r$y[a] - r$y[b] ) ^ 2 ), sqrt( ( r$x[a] - r$x[c] ) ^ 2 + ( r$y[a] - r$y[c] ) ^ 2 ), sqrt( ( r$x[b] - r$x[c] ) ^ 2 + ( r$y[b] - r$y[c] ) ^ 2 ))
  
  
  # perimeter of triangles
  tri_perimeter <- rowSums(tri_sidelength)
  
  # area of triangles
  s <- 0.5*tri_perimeter
  tri_area <- sqrt(s* (s - tri_sidelength[,1]) * (s - tri_sidelength[,2]) * (s - tri_sidelength[,3]))
  
  # clustering
  k <- seq_len( r$tlnew - 1 )
  i <- r$tlist[k]          
  j <- r$tlist[r$tlptr[k]]
  keep <- i > 0
  i <- abs( i[ keep ] )
  j <- abs( j[ keep ] )
  distances <- sqrt( ( r$x[i] - r$x[j] ) ^ 2 + ( r$y[i] - r$y[j] ) ^ 2 )
  threshold <- 20  # Choose the threshold manually
  i <- i[ distances < threshold ]
  j <- j[ distances < threshold ]
  #--------------------------------------------#
  
  Dual_NodeList <- Coords_shapeDes
  colnames(Dual_NodeList) <- c('cellID/nodes', 'x', 'y', 'orientation')
  
  
  Dual_EdgeList <- data.frame(cbind(i, j))
  colnames(Dual_EdgeList) <- c('col1', 'col2')
  
  Dual_EdgeList <- Dual_EdgeList %>%
    mutate(from = pmin(col1, col2), 
           to = pmax(col1, col2)) %>%
    distinct(from, to)
  
  
  colnames(Dual_EdgeList) <- c('from', 'to')
  
  # generate initial network
  g <- graph_from_data_frame(vertices = Dual_NodeList, d= Dual_EdgeList, directed = FALSE)
  
  a <- degree(g)
  
  plot(g, vertex.size = 1, vertex.label = NA)
  # clustering
  comp <- components(g)
  comp_list <- lapply(seq_along(comp$csize)[comp$csize > 1], function(x) 
    V(g)$name[comp$membership %in% x])
  #--------------------------------------------#
  
  
  # remove low density trees
  list_id <- 1
  clus_id <- 1
  
  
  
  
  Area <- matrix(nrow = 0, ncol = 1)
  posDensity <- matrix(nrow = 0, ncol = 1)
  Clus_density <- matrix(nrow = 0, ncol = 1)
  triStat <- matrix(nrow = 0, ncol = 9)
  cluster_posListAll <- matrix(nrow = 0, ncol = 5)
  
  while(list_id <= length(comp_list)){
    cluster_posList <- matrix(nrow = 0, ncol = 4)
    if(length(comp_list[[list_id]]) >= 30){
      print(list_id)
      
      for(item in comp_list[[list_id]]){
        # find all rows which contain the member
        
        # get the member pos, unit: mm^2
        item_x <- as.numeric(Dual_NodeList[Dual_NodeList$`cellID/nodes` == item,]['x'])/1000
        item_y <- as.numeric(Dual_NodeList[Dual_NodeList$`cellID/nodes` == item,]['y'])/1000
        
        cluster_posList <- rbind(cluster_posList, cbind(item_x, item_y, clus_id, item))
      }
      cluster_posList <- cluster_posList[complete.cases(cluster_posList),]
      colnames(cluster_posList) <- c('Centroid.X.mm', 'Centroid.Y.mm', 'cluster', 'cell_ID')
      
      # merge
      cluster_posList <- merge(cluster_posList, Coords_shapeDes, by= 'cell_ID')
      
      cluster_posList$Centroid.X.mm <- as.numeric(as.character(cluster_posList$Centroid.X.mm))
      cluster_posList$Centroid.Y.mm <- as.numeric(as.character(cluster_posList$Centroid.Y.mm))
      # get the cluster size
      concaveHull <- concaveman(data.matrix(cluster_posList[,2:3]), concavity = 2)
      area <- gArea(SpatialPolygons(list(Polygons(list(Polygon(concaveHull[,1:2])),1))))
      Area <- rbind(Area, area)
      
      
      
      # get the density
      posDensity <- rbind(posDensity, nrow(cluster_posList)/area)
      
      
      
      # record the original cluster index for future referral
      
      clus_id <- clus_id + 1
      list_id <- list_id + 1
      
      cluster_posListAll <- rbind(cluster_posListAll, cluster_posList)
      
    } else {
      list_id <- list_id + 1
    }
  }
  

  if(length(cluster_posListAll) != 0){
    # remove redundant columns
    cluster_posListAll <- cluster_posListAll[, c(2,3,4,7)]
    
    
    # calculate COrE features
    COrE_features <- getCOrE(cluster_posListAll, 18)
    
    COrE_feature.names <- colnames(COrE_features)
    
    
    clus_mean <- mean(posDensity) # average nucleus density per cluster
    clus_max <- max(posDensity) # average nucleus density per cluster
    clus_min <- mean(posDensity) # average nucleus density per cluster
    clus_CoV <- sd(posDensity)/clus_mean # average nucleus density per cluster
    
    
    
    # triangulation stats
    triStat <- rbind(triStat, cbind(Num_tri, mean(tri_perimeter), max(tri_perimeter), min(tri_perimeter), sd(tri_perimeter)/mean(tri_perimeter),
                                    mean(tri_area), max(tri_area), min(tri_area), sd(tri_area)/mean(tri_area)))
    
    # get the cluster density
    tissueArea <- TMA_allArea[TMA_allArea$TMA.core == core,2]
    
    
    
    # core level
    this.core <- cbind(core, clus_mean, clus_max, clus_min, clus_CoV, (clus_id - 1)/tissueArea, triStat, COrE_features)
    colnames(this.core) <- c('TMA.core', 'Density of nucleus_mean', 'Density of nucleus_max',  'Density of nucleus_min',  'Density of nucleus_CoV', 'Density of clusters', 'Number of triangles',
                             'Perimeter of triangle_mean', 'Perimeter of triangle_max', 'Perimeter of triangle_min', 'Perimeter of triangle_CoV', 'Area of triangle_mean',
                             'Area of triangle_max', 'Area of triangle_min', 'Area of triangle_CoV', COrE_feature.names)
    
    coreClusStat <- rbind(coreClusStat, this.core)
    
  } else {
    
    NA_dat <- data.frame(cbind(core, t(rep(NA, 53))))
    colnames(NA_dat) <- c('TMA.core', 'Density of nucleus_mean', 'Density of nucleus_max',  'Density of nucleus_min',  'Density of nucleus_CoV', 'Density of clusters', 'Number of triangles',
                          'Perimeter of triangle_mean', 'Perimeter of triangle_max', 'Perimeter of triangle_min', 'Perimeter of triangle_CoV', 'Area of triangle_mean',
                          'Area of triangle_max', 'Area of triangle_min', 'Area of triangle_CoV', COrE_feature.names)
    coreClusStat <- rbind(coreClusStat, NA_dat)
  }
  
} 
write.csv(coreClusStat, 'clusterStat.csv')



# heatmap

write.csv(co_occurrence, 'co_occurrence.csv')
