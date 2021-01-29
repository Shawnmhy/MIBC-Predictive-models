############################################################
#### This script is to provide functions for Feature.R #####
############################################################
library(rjson)
library(rgeos)
library(reshape2)
library(fractaldim)
#########################
# get areas for polygon #
#########################

# Label: 1
getSmoothCellShape <- function(xy, smooth_pts) {
  # convet to spatial feature
  Sp <- SpatialPolygons(list(Polygons(list(Polygon(xy)),1)))
  smooth_poly <- smooth(Sp, method='spline', n= smooth_pts)
  cellShape_smooth <- data.frame(smooth_poly@polygons[[1]]@Polygons[[1]]@coords)
  # return the coordinates of smoothed boundary
  colnames(cellShape_smooth) <- c('x', 'y')
  return(cellShape_smooth = cellShape_smooth)
}

# Label: 2
getConvexHull <- function(xy){
  
  
  # get convex hull
  cellConvexHull <- gConvexHull(SpatialPolygons(list(Polygons(list(Polygon(xy)),1))))@polygons[[1]]@Polygons[[1]]@coords
  # convert to spatial feature
  Sp <- SpatialPolygons(list(Polygons(list(Polygon(cellConvexHull)),1)))
  
  chullLength <- gLength(SpatialPolygons(list(Polygons(list(Polygon(cellConvexHull)),1))))
  
  chullArea <- gArea(SpatialPolygons(list(Polygons(list(Polygon(cellConvexHull)),1))))

  return(cbind(chullLength, chullArea))
}

getBaselineStat <- function(subDF, metric){
  
  subset<- subDF[eval(metric)]
  mean_stat <- mean(subset[,1])
  max_stat <- max(subset[,1])
  min_stat <- min(subset[,1])
  
  if(mean_stat != 0){
    cov_stat <- sd(subset[,1])/mean_stat  
  } else {
    cov_stat <- 0
  }
  
  
  # to matrix
  stat <- cbind(max_stat, min_stat, mean_stat, cov_stat)
  # rename
  colnames(stat) <- c(paste(metric, '_Max',sep = ''), paste(metric, '_Min',sep = ''), paste(metric, '_Mean',sep = ''), paste(metric, '_CoV',sep = ''))
  return(stat)
}
####################
# get bounding box #
####################


# Minimum Bounding Box function
getMinBBox <- function(xy) {
  
  # convert to matrix
  xy <- as.matrix(xy)
  stopifnot(is.matrix(xy), is.numeric(xy), nrow(xy) >= 2, ncol(xy) == 2)
  
  ## rotating calipers algorithm using the convex hull
  H    <- chull(xy)      ## hull indices, vertices ordered clockwise
  n    <- length(H)      ## number of hull vertices
  hull <- xy[H, ]        ## hull vertices
  
  ## unit basis vectors for all subspaces spanned by the hull edges
  hDir  <- diff(rbind(hull, hull[1, ])) ## hull vertices are circular
  hLens <- sqrt(rowSums(hDir^2))        ## length of basis vectors
  huDir <- diag(1/hLens) %*% hDir       ## scaled to unit length
  
  ## unit basis vectors for the orthogonal subspaces
  ## rotation by 90 deg -> y' = x, x' = -y
  ouDir <- cbind(-huDir[ , 2], huDir[ , 1])
  
  ## project hull vertices on the subspaces spanned by the hull edges, and on
  ## the subspaces spanned by their orthogonal complements - in subspace coords
  projMat <- rbind(huDir, ouDir) %*% t(hull)
  
  ## range of projections and corresponding width/height of bounding rectangle
  rangeH  <- matrix(numeric(n*2), ncol=2)  ## hull edge
  rangeO  <- matrix(numeric(n*2), ncol=2)  ## orthogonal subspace
  widths  <- numeric(n)
  heights <- numeric(n)
  
  for(i in seq(along=numeric(n))) {
    rangeH[i, ] <- range(projMat[  i, ])
    
    ## the orthogonal subspace is in the 2nd half of the matrix
    rangeO[i, ] <- range(projMat[n+i, ])
    widths[i]   <- abs(diff(rangeH[i, ]))
    heights[i]  <- abs(diff(rangeO[i, ]))
  }
  
  ## extreme projections for min-area rect in subspace coordinates
  ## hull edge leading to minimum-area
  eMin  <- which.min(widths*heights)
  hProj <- rbind(   rangeH[eMin, ], 0)
  oProj <- rbind(0, rangeO[eMin, ])
  
  ## move projections to rectangle corners
  hPts <- sweep(hProj, 1, oProj[ , 1], "+")
  oPts <- sweep(hProj, 1, oProj[ , 2], "+")
  
  ## corners in standard coordinates, rows = x,y, columns = corners
  ## in combined (4x2)-matrix: reverse point order to be usable in polygon()
  ## basis formed by hull edge and orthogonal subspace
  basis <- cbind(huDir[eMin, ], ouDir[eMin, ])
  hCorn <- basis %*% hPts
  oCorn <- basis %*% oPts
  pts   <- t(cbind(hCorn, oCorn[ , c(2, 1)]))
  
  ## angle of longer edge pointing up
  dPts <- diff(pts)
  e    <- dPts[which.max(rowSums(dPts^2)), ] ## one of the longer edges
  eUp  <- e * sign(e[2])       ## rotate upwards 180 deg if necessary
  deg  <- atan2(eUp[2], eUp[1])*180 / pi     ## angle in degrees
  

  
  return(list(pts=pts, width=widths[eMin], height=heights[eMin], angle=deg))
}

#################
# get curvature #
#################

# input is smoothed boundary
getCurvature <- function(xy, ptsCount) {
  
  #xy <- smoothedShape
  # calculate curvature at each point
  dx <- diff(c(xy$x, xy$x[1])) # Distance between x coords with wrap-around
  dy <- diff(c(xy$y, xy$y[1])) # Distance between y coords with wrap-around
  ds <- sqrt(dx^2 + dy^2)                    # Segment size between points
  ddx <- dx/ds                               # Ratio of x distance to segment size
  ddy <- dy/ds                               # Ratio of y distance to segment size
  ds2 <- (ds + c(ds[-1], ds[1]))/2           # Mean segment length either side per point
  xy$Cx <- diff(c(ddx, ddx[1]))/ds2   # Change in ddx per unit length
  xy$Cy <- diff(c(ddy, ddy[1]))/ds2   # Change in ddy per unit length
  xy$K <- (ddy * xy$Cx - ddx * xy$Cy)/
    ((ddx^2 + ddy^2)^(3/2))
  
  
  xy <- xy[complete.cases(xy), ]
  
  # 
  curvMean <- mean(xy$K)
  curvMin <- min(xy$K)
  curvMax <- max(xy$K)
  curvCoV <- sd(xy$K)/curvMean
  # get number of protrusion:
  n_protrusion <- 0
  n_indentation <- 0
  for(i in 1:nrow(xy)){
    if(xy$K[i] > 2){
      # scenario 1:
      if(i >= 5 & i <= (nrow(xy) - 5)){
        neighborhood <- xy$K[(i -5):(i+5)]
      } else if(1 < i & i <= 5){
        neighborhood <- xy$K[c((nrow(xy) - 5 + i):nrow(xy),1:i, i:(i+4))]
      } else if( nrow(xy) > i & (i + 5) > nrow(xy)){
        neighborhood <- xy$K[c((i-5):i, (i+1):nrow(xy), 1:(5 - nrow(xy) + i))]
      } else if(nrow(xy) == i){
        neighborhood <- xy$K[c((nrow(xy) - 5):nrow(xy),1:6)]
      } else if(1 == i){
        neighborhood <- xy$K[c((nrow(xy)-6):nrow(xy),1:5)]
      }
      if(max(neighborhood) == xy$K[i]){
        n_protrusion <- n_protrusion + 1
      }
    }
    if(xy$K[i] < -2){
      # scenario 1:
      if(i >= 5 & i <= (nrow(xy) - 5)){
        neighborhood <- xy$K[(i -5):(i+5)]
      } else if(1 < i & i <= 5){
        neighborhood <- xy$K[c((nrow(xy) - 5 + i):nrow(xy),1:i, i:(i+4))]

      } else if( nrow(xy) > i & (i + 5) > nrow(xy)){
        neighborhood <- xy$K[c((i-5):i, (i+1):nrow(xy), 1:(5 - nrow(xy) + i))]

      } else if(nrow(xy) == i){
        neighborhood <- xy$K[c((nrow(xy) - 5):nrow(xy),1:5)]

      } else if(1 == i){
        neighborhood <- xy$K[c((nrow(xy)-5):nrow(xy),1:6)]

      }
      if(max(neighborhood) == xy$K[i]){
        n_indentation <- n_indentation + 1
      }
    }
}
  return(cbind(curvMean, curvMin, curvMax, curvCoV, n_protrusion, n_indentation))
}

getFractalD <- function(xy) {
  fd2d <- fd.estim.boxcount(cbind(xy[,1],xy[,2]),plot.log=FALSE, plot.allpoints=FALSE, nlags="auto")
  FractalD <- as.numeric(fd2d['fd']) # fractal dimension
  return(FractalD)
}

getCOrE <- function(orientationList, w){
  
  orientationList <- data.frame(orientationList)
  
  # angle diretization
  N <- 360/w
  # define a small positive number to avoid log error
  eps <- 10^(-15)
  orientationList$orientationAngle <- ceiling((orientationList$orientationAngle)/w)*w
  
  # create co-occurence matrix
  co_occurrence <- matrix(nrow = N, ncol = N)
  
  for(col in 1:N){
    for(row in 1:N){
      if(col == row){
        co_occurrence[row, col] <- choose(length(orientationList[orientationList$orientationAngle == col*w,]),2) 
      } else {
        co_occurrence[row, col] <- length(orientationList[orientationList$orientationAngle == col*w,])*
          length(orientationList[orientationList$orientationAngle == row*w,])
      }
      
    }
  }
  p_occurance <- data.frame(co_occurrence/sum(co_occurrence))
  # should be a number
  weighted_matrix <-  p_occurance
  for(col in 1:N){
    weighted_matrix[,col] <- p_occurance[,col]*col
  }
  # calculate marginal-mean
  ux <- sum(rowSums(weighted_matrix))
  uy <- ux
  
  # calculate marginal-sigma
  sigma_x <- 0
  for(col in 1:N){
    for(row in 1:N){
      # entropy
      sigma_x <- sigma_x + p_occurance[row,col]*(row - ux)^2
    }
  }
  
  sigma_y <- sigma_x
  
  # calculate some features
  COrE_entropy <- 0 # 1
  COrE_energy <- 0 # 2
  COrE_IDM <- 0 # 3, homogeneity 
  COrE_dissimilarity <- 0 # 4
  
  COrE_corrln <- 0 # 5
  
  # px and py
  px <- colSums(p_occurance)
  py <- px
  
  # calculate Homogeneity 1 and 2
  HX <- -sum(px*log(px+eps))
  HY <- HX
  HXY1 <- 0
  HXY2 <- 0
  # entropy
  for(col in 1:N){
    for(row in 1:N){
      # entropy
      COrE_entropy <- COrE_entropy-( p_occurance[row, col]*log(p_occurance[row,col]+eps))
      COrE_energy <- COrE_energy + p_occurance[row, col]^2
      COrE_IDM <- p_occurance[row,col]/(1 +(row - col)^2)
      COrE_corrln <- COrE_corrln +  p_occurance[row,col]*((row - ux)*(col - uy)/sqrt(sigma_x*sigma_y))
      COrE_dissimilarity <- COrE_dissimilarity + p_occurance[row,col]*abs(row - col)
      HXY <- COrE_entropy
      HXY1 <- HXY1  - p_occurance[row,col]*log(px[row]*py[col]+eps)
      HXY2 <- HXY2 - px[row]*py[col]*log(px[row]*py[col]+eps)
    }
  }
  COrE_energy <- sqrt(COrE_energy)
  COrE_Homogeneity1 <- (COrE_entropy - HXY1)/max(HX, HY) # 6
  COrE_Homogeneity2 <- (1 - exp(-2*(HXY2 - COrE_entropy)))^0.5 #7
  # for diff entropy
  # flag for Diff Entropy
  k = 0
  px_y <- 0
  COrE_diffEntropy <- 0 # 8
  COrE_contrast <- 0 # 9
  diffVariance_list <- matrix(nrow = 0, ncol = 1) # 10
  while(k < N){
    for(col in 1:N){
      for(row in 1:N){
        # calculate Difference entropy
        if(abs(col - row) == k){
          px_y <- px_y + 1*p_occurance[row, col]
        } else {
          px_y <- px_y # unchanged
        }
      }
    }
    COrE_diffEntropy <- COrE_diffEntropy - px_y*log(px_y+eps)
    diffVariance_list <- rbind(diffVariance_list, px_y)
    COrE_contrast <- COrE_contrast + k^2*px_y
    k <- k + 1
    px_y <- 0
  }
  
  COrE_diffVariance <- as.numeric(var(diffVariance_list))
  
  k <- 2
  pxy <- 0
  COrE_SumAverage <- 0
  COrE_SumEntropy <- 0
  COrE_SumVariance <- 0
  while(k <= 2*N){
    for(col in 1:N){
      for(row in 1:N){
        # calculate Difference entropy
        if((col+row) == k){
          pxy <- pxy + 1*p_occurance[row, col]
        } else{
          pxy <- pxy # unchanged
        }
      }
    }
    COrE_SumAverage <- COrE_SumAverage + k*pxy # 11
    COrE_SumEntropy <- COrE_SumEntropy - pxy*log(pxy+eps) # 12
    k <- k + 1
    pxy <- 0
  }
  
  pxy <- 0
  k <- 2
  while(k <= 2*N){
    for(col in 1:N){
      for(row in 1:N){
        # calculate Difference entropy
        if((col+row) == k){
          pxy <- pxy + 1*p_occurance[row, col]
        } else{
          pxy <- pxy # unchanged
        }
      }
    }
    COrE_SumVariance <- COrE_SumVariance + pxy*(k - COrE_SumEntropy)^2 #13
    k <- k + 1
    pxy <- 0
    
  }
  
  COrE_features <- data.frame(COrE_corrln, COrE_IDM, COrE_dissimilarity, COrE_entropy, COrE_energy, COrE_diffEntropy, COrE_contrast, COrE_diffVariance, COrE_Homogeneity1, COrE_Homogeneity2, COrE_SumAverage, COrE_SumEntropy, COrE_SumVariance)
  return(COrE_features)
}

getCluster <- function(clusterDat){
  # unit: in mm
  X <- as.matrix(clusterDat[,1:2])/1000
  
  # get the cluster area, concave hull
  clusConcave <- concaveman(X, concavity = 2)
  # unit: mm^2
  clusArea  <- gArea(SpatialPolygons(list(Polygons(list(Polygon(clusConcave)),1))))
  
  # get the perimeter
  clusPeri <- gLength(SpatialPolygons(list(Polygons(list(Polygon(clusConcave)),1))))
  # get the cluster cell density
  cellDen <- nrow(X)/clusArea
  
  ## fitting ellipse
  # we can get center and covariance matrix from the cluster points:
  ellFit <- cov.wt(X)
  # eigen vectors indicate orientations
  ellEigen <- eigen(ellFit$cov)
  
  # eigen values are associated with major/minor axes length (half)
  axes <- sqrt(ellEigen$values * qchisq(.95, df=2))
  
  # get the cluster circularity
  clusCirc <- 4*pi*clusArea/(clusPeri)^2
  
  # get the cluster eccentricity
  clusEccen <- sqrt(1-(axes[2]/axes[1])^2)
  
  # get the moment of inersia
  clusMOI <- sum(sweep(X,2,colMeans(X))^2)
  
  # get the centroids of cluster
  cent <- t(matrix(colMeans(X)))
  
  clusterFeatures <- cbind(cent, clusArea, clusPeri, cellDen, as.matrix(t(c(axes[1], axes[2]))), clusCirc, clusEccen, clusMOI)
  colnames(clusterFeatures) <- c('clus_x', 'clus_y', 'clus_area', 'clus_perimeter', 'cell_density', 'major_length', 'minor_length', 'clus_circ', 'clus_eccenc','clus_moi')
  return(clusterFeatures)
}







