##############################
## fitting point process model
##############################

# load spatstat package
library(spatstat) 
# test CSR, return true if reject
testRejectCSR <- function(pppPattern, alpha){
  if(pppPattern$n > 3){
    # Quatrat test
    # not used
    #n = floor(sqrt(pppPattern$n + 1 / 2))
    #qTest <- quadrat.test(pppPattern, nx = n, ny = n)
    #if (qTest$p.value < alpha){
    #  return (TRUE)
    #}
    ceResult = clarkevans.test(pppPattern, alternative='clustered',nsim=100)
    if(ceResult$p.value<p_th)
      return(TRUE)
  }
  return (FALSE)
}

# to construct polygon owin
tissueList <- function(path, core){
  
  List <- list()
  # core folder 
  #core <- 'A-1'
  coreFolder <- paste(path, core, sep = '')
    
  # count the number of files in the folder
  fileCount <- length(dir(coreFolder))

  if(fileCount !=0){
    for(i in 1:fileCount){
      tissueFile  <- data.frame(t(data.frame(fromJSON(file = paste(coreFolder, '/tissue', i, '.json', sep = ''))$geometry$coordinate[[1]])))
      # reformat
      tissueFile$X1 <- rev(as.numeric(as.character(tissueFile$X1)))
      tissueFile$X2 <- rev(as.numeric(as.character(tissueFile$X2)))
      
      colnames(tissueFile) <- c('x', 'y')
      List[[i]] <- tissueFile
      }
  }
  return(List)
}

# to rescale points
ptsRescale <- function(List, IHC, core){
  
  for(i in 1:length(List)){
    
    # get the pts
    ptsFile <- List[[i]]
    
    # read core dat
    #IHC <- 'CK56'
    metaFile <- read.delim(paste('./IHC_Core_Centroids/', IHC, '.txt', sep = ''), header = T,sep = '\t')
    name <- str_split_fixed(metaFile$Image, '_', n = 6)
    metaFile$Image <- name[,4]
    
    # get the reference coordinates
    ref_x <- metaFile[metaFile$Name == core,]$Centroid.X.um/0.454 - 3030/2 # unit: pixel
    ref_y <- metaFile[metaFile$Name == core,]$Centroid.Y.um/0.454 - 3030/2 # unit: pixel
    
    # rescaling, if the unit is pixel, then do not /0.454
    rescale_x <- ptsFile$x - ref_x
    rescale_y <- ptsFile$y - ref_y  
    rescale_pts <- cbind(rescale_x, rescale_y)
    
    List[[i]] <- rescale_pts
  }
  return(List)
}



# fit subregion to point process model
subRegionFit <- function(points, xwindow, ywindow, alpha){
  n = nrow(points)
  intensity = n/xwindow/ywindow
  p1 = c(n, intensity)
  p2 = c(0, 0, 0, 0,0)
  if(n>0){
    pattern <- ppp(points[,1], points[,2], c(0,xwindow), c(0,ywindow))
    #print(pattern)
    # test CSR
    if(testRejectCSR(pattern,p_th)){
      p2[1] = 1
      p2[2] = 1
      tryCatch( {
        myfit = kppm(pattern, ~1, "Thomas")
        p = parameters(myfit)
        p2[3:5] = c(p$kappa, p$scale, p$mu)
      } , warning = function(w) {
      }, error = function(e) {
      }, finally = {
      })
    }
    else{
      p2[1] = 0
    }
  }
  return(c(p1, p2))
}
# process model fitting for one slide
# get subregions of the slide using a moving window
# fit point pattern in sub regions to spatial process model
# record fitted parameters

# fit spatial point process model and record fitted parameters
spatstatFit <- function(subdir, alpha){
  
  # path to files and load data
  path = paste(subdir, "\\", sep="")
  print(path)
  filePattern = "^\\d+_\\d+.+csv$"
  patchCoords = list.files(path=path, pattern = filePattern
                           , recursive=FALSE, include.dirs = TRUE
                           , full.names = FALSE)
  
  # 7 columns for p_all: x, y, intensity, CSR? kappa, sigma2, mu
  p_all <- matrix(nrow=0,ncol=7)
  header = c("x", "y", "intensity", "notCSR", "kappa", "sigma2", "mu")
  colnames(p_all) = header
  
  for (patch in patchCoords) { #i: x
    coord = as.integer( strsplit(patch,"[_.]")[[1]][1:2])
    i = coord[1]
    j = coord[2]
    # loop through grid
    #print(filename)
    filePath = paste(path, patch, sep="")
    mydata <- read.csv(filePath)
    if(dim(mydata[1])>0){
      mypattern <- ppp(mydata[,1], mydata[,2], c(0,window), c(0,window))
      intensity = nrow(mydata)/window^2
      p1 = c(i, j, intensity)
      p2 = c(0, 0, 0, 0)
      # test CSR
      if(testRejectCSR(mypattern, alpha=alpha)){
        p2[1] = 1
        tryCatch( {
          myfit = kppm(mypattern, ~1, "Thomas")
          p = parameters(myfit)
          p2[2:4] = c(p$kappa, p$scale, p$mu)
        } , warning = function(w) {
        }, error = function(e) {
        }, finally = {
        })
      }
      else{
        p2[4] = 0
      }
      p_all = rbind(p_all, c(p1, p2))
    }
  }
  
  # write results to file
  write.csv(p_all, file = paste(path, "fittedResult_newCSR.csv", sep=""),
            row.names=FALSE)
}

##############################
### batch shape descriptor
##############################


# packages
library(largeVis)
library(alphahull)
library(stats)

# class AshapePol:
# two member variables:
# goodAshape: TRUE/FALSE
# vert: indices of vertices
setClass(Class="AshapePol",
         representation(
           goodAshape="logical",
           vert="numeric"
         )
)

## get the polygon boundary points in right order
library(sp)
library(igraph)

getAlphaShapeInOrder <- function(ashapeX){
  returnVal=new("AshapePol", goodAshape=TRUE, vert=0)
  
  if (nrow(ashapeX$edges)==0) {
    #stop("Graph not connected")
    returnVal@goodAshape = FALSE
  }
  
  ashapeGraph = graph.edgelist(cbind(as.character(ashapeX$edges[, "ind1"]), 
                                     as.character(ashapeX$edges[, "ind2"])), directed = FALSE)
  
  
  
  #plot(ashapeGraph)
  ## check if is closed
  if (!is.connected(ashapeGraph)) {
    #stop("Graph not connected")
    returnVal@goodAshape = FALSE
  }
  if (any(degree(ashapeGraph) != 2)) {
    #stop("Graph not circular")
    returnVal@goodAshape = FALSE
  }
  if (clusters(ashapeGraph)$no > 1) {
    #stop("Graph composed of more than one circle")
    returnVal@goodAshape = FALSE
  }
  
  if(! returnVal@goodAshape){
    return(returnVal)
  }
  
  cutg = ashapeGraph - E(ashapeGraph)[1]
  # find chain end points
  ends = names(which(degree(cutg) == 1))
  path = get.shortest.paths(cutg, ends[1], ends[2])$vpath[[1]]
  # this is an index into the points
  pathX = as.numeric(V(ashapeGraph)[path]$name)
  # join the ends
  pathX = c(pathX, pathX[1])
  
  aShapePath <- Polygon(ashapeX$x[pathX, ])  
  # polygon(X[pathX, 1], X[pathX,2], col = "gray", border = "red")
  
  ## test within polygon
  inAshape = point.in.polygon(ashapeX$x[,1],ashapeX$x[,2],ashapeX$x[pathX,1],ashapeX$x[pathX,2])
  if (any(inAshape==0)){
    #stop("point outside alpha shape")
    returnVal@goodAshape = FALSE
  }else{
    returnVal@vert = pathX
  }
  return(returnVal)
}


pointCharacterize <- function(coords_data, invasive_poly, tumor_poly){
  normal_set <-  matrix(nrow=,ncol=2)
  invasive_set <-  matrix(nrow=0,ncol=2)
  tumor_set <-  matrix(nrow=0,ncol=2)
  
  len <- nrow(coords_data)
  for (idx in seq(1,len,1)){
    x_coords <- coords_data[idx,1]
    y_coords <- coords_data[idx,2]
    Coords <- cbind(x_coords, y_coords)
    # Note: how to determine the boundary condition
    if(point.in.polygon(x_coords, y_coords, invasive_poly[,1], invasive_poly[,2]) != 0){
      invasive_set <- rbind(invasive_set,Coords)
    }
    else if(point.in.polygon(x_coords, y_coords, tumor_poly[,1], tumor_poly[,2]) == 1){
      tumor_set <- rbind(tumor_set,Coords)
    }
    else{
      normal_set <- rbind(normal_set,Coords)
    }
    
  }
  Charac_point <- list('normal' = normal_set, 'tumor' = tumor_set, 'invasive' = invasive_set)
  
  return(Charac_point)
}





bivarAnalysis <- function(IHC_marker1, IHC_marker2, Region){
  
  #IHC_marker1 <- 'Her2Neu'
  #IHC_marker2 <- 'GATA3'
  
  type1 <- get(eval(IHC_marker1))[,1:2]
  colnames(type1) <- c('x', 'y')
  
  type2 <- get(eval(IHC_marker2))[,1:2]
  colnames(type2) <- c('x', 'y')
  
  if(nrow(type1)*nrow(type2) != 0){
    
    # read pts dat
    #Region <- Region_HE

    type1$IHC <- IHC_marker1

    type2$IHC <- IHC_marker2
    
    # create multitype df
    pts_OI <- rbind(type1, type2)
    
    # define the type
    species <- factor(pts_OI$IHC)
    
    # create multitype ppp
    #Region <- Region_CK56
    
    
    # check if empty  
    ppp1 <- ppp(type1$x, type1$y, owin(poly = Region))
    ppp2 <- ppp(type2$x, type2$y, owin(poly = Region))
    
    # prevent NA 

    if(is.empty(ppp1) == 'FALSE' & is.empty(ppp2) == 'FALSE'){
      multitype_ppp <- ppp(pts_OI$x, pts_OI$y, marks = species, owin(poly = Region))
      Gihc <- data.frame(Gcross(multitype_ppp, i = IHC_marker1, j = IHC_marker2))
      # relocat DF

      Gihc <- Gihc[, c(1,2,4)]
      
      Gihc <- Gihc[complete.cases(Gihc),]
      
      # get the 'dense distance'
      diff <- diff(Gihc$theo)/diff(Gihc$r)
      common.distance1 <- Gihc$r[which.max(diff)]
      
      
      # calculate the area (positive - negative )  
      Gihc$rs <- Gihc$rs - Gihc$theo
      
      apprx_index1 <- trapz(Gihc$r, Gihc$rs) 
      
      # 
      Gihc <- data.frame(Gcross(multitype_ppp, i = IHC_marker2, j = IHC_marker1))
      
      Gihc <- Gihc[, c(1,2,4)]
      
      Gihc <- Gihc[complete.cases(Gihc),]
      
      
      diff <- diff(Gihc$theo)/diff(Gihc$r)
      
      common.distance2 <- Gihc$r[which.max(diff)]
      
      # calculate the area (positive - negative )  
      Gihc$rs <- Gihc$rs - Gihc$theo
      Gihc <- Gihc[complete.cases(Gihc),]
      apprx_index2 <- trapz(Gihc$r, Gihc$rs) 
      
      
    } else {
      
      common.distance1 <- 0
      apprx_index1 <- 0
      
      common.distance2 <- 0
      apprx_index2 <- 0
    }
  } else {
    common.distance1 <- 0
    apprx_index1 <- 0
    
    common.distance2 <- 0
    apprx_index2 <- 0
  }
  return(list(0.5*(apprx_index1 + apprx_index2), 0.5*(common.distance1 + common.distance2)))
}
