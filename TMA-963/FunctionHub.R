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


#######################
# Kcross function #####
#######################



bivarAnalysis.Kcross <- function(pts.type1, pts.type2, Region){
  
 # Region <- Region_CK56
  #pts.type1 <- 'Cancer'
  #pts.type2 <- 'lymphocytes'
  
  type1 <- get(eval(pts.type1))[,2:3]
  colnames(type1) <- c('x', 'y')
  
  type2 <- get(eval(pts.type2))[,2:3]
  colnames(type2) <- c('x', 'y')
  
  if(nrow(type1)*nrow(type2) != 0){
    
    # read pts dat
    #Region <- Region_HE

    type1$attr <- pts.type1

    type2$attr <- pts.type2
    
    # create multitype df
    pts_OI <- rbind(type1, type2)
    
    # define the type
    species <- factor(pts_OI$attr)
    
    # create multitype ppp
    #Region <- Region_CK56
    
    
    # check if empty  
    ppp1 <- ppp(type1$x, type1$y, owin(poly = Region))
    ppp2 <- ppp(type2$x, type2$y, owin(poly = Region))
    
    # prevent NA 

    if(is.empty(ppp1) == 'FALSE' & is.empty(ppp2) == 'FALSE'){
      
    
      
      
      multitype_ppp <- ppp(pts_OI$x, pts_OI$y, marks = species, owin(poly = Region))
      K.cross <- data.frame(Kcross(multitype_ppp, i = pts.type1, j = pts.type2, r = seq(0,20,0.1), correction = 'Ripley'))
      
      #plot(Gihc)
      # relocat DF

      K.cross <- K.cross[complete.cases(K.cross),]

      
      
      # calculate the area (positive - negative )  
      K.cross$rs <- K.cross$iso - K.cross$theo
      
      i.to.j.diff.area <- trapz(K.cross$r, K.cross$iso) 
      
      
      # j to i
      
      K.cross <- data.frame(Kcross(multitype_ppp, i = pts.type2, j = pts.type1, r = seq(0,20,0.1), correction = 'Ripley'))
      
      K.cross <- K.cross[complete.cases(K.cross),]
    
      
      
      K.cross$iso <- K.cross$iso - K.cross$theo
      K.cross <- K.cross[complete.cases(K.cross),]
      j.to.i.diff.area <- trapz(K.cross$r, K.cross$iso) 
      
      
    } else {
      
      i.to.j.diff.area <- 0
      j.to.i.diff.area <- 0

    }
  } else {
    i.to.j.diff.area <- 0
    j.to.i.diff.area <- 0
  }
  return(list(i.to.j.diff.area, j.to.i.diff.area))
}



#######################
# Gcross function #####
#######################

bivarAnalysis <- function(IHC_marker1, IHC_marker2, Region){
  
  # Region <- Region_CK56
  #IHC_marker1 <- 'CK20'
  #IHC_marker2 <- 'CK56'
  
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
      
      # number of points within range
      nn.dist.IHC2centered <- nn2(type1[,1:2], query = type2[,1:2], k = nrow(type1), searchtype = 'radius', radius = 20)
      
      nn.dist.IHC1centered <- nn2(type2[,1:2], query = type1[,1:2], k = nrow(type2), searchtype = 'radius', radius = 20)
      
      
      mtrx.IHC1centered <- mean(rowSums(data.frame(nn.dist.IHC1centered$nn.idx) != 0))
      mtrx.IHC2centered <- mean(rowSums(data.frame(nn.dist.IHC2centered$nn.idx) != 0))
      
      
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
      mtrx.IHC1centered <- 0
      
      common.distance2 <- 0
      apprx_index2 <- 0
      mtrx.IHC2centered <- 0
    }
  } else {
    common.distance1 <- 0
    apprx_index1 <- 0
    mtrx.IHC1centered <- 0
    
    common.distance2 <- 0
    apprx_index2 <- 0
    mtrx.IHC2centered <- 0
  }
  return(list(0.5*(apprx_index1 + apprx_index2), 0.5*(common.distance1 + common.distance2), mtrx.IHC1centered, mtrx.IHC2centered))
}




#--------------------------------------------------------------#
#----------------------- SVM  Classifier ----------------------#
#--------------------------------------------------------------#

# ------- Feature selection by t test ---------#



svmClassifier <- function(Data.matrix, fold, kernel){
  
  #Data.matrix <- allFeatures_forModel.spatstat
  #fold <- 5
  set.seed(1)
  inner.cv <- matrix(nrow = 0, ncol = 3)
  
  #kernel <- 'radial'
  folds = createFolds(Data.matrix$Outcome, k = fold)
  # machine learning
  
  
  
  optimization<- matrix(nrow = 0, ncol = 3)
  
  optimized.feature <- list()    
  #glm.fit_list <- list()
  
  iterations <- 0
  
  
  for(gamma in 10^(-10:3)){
    for(cost in 1:100){
      
      Accuracy <- matrix(nrow = 0, ncol = 1)
      
      for (id in 1:length(folds)) {
        #id <- 5
        x = folds[[id]]
        
        # in the next two lines we will separate the Training set into it's 10 pieces
        training_fold = Data.matrix[-x,] # training fold =  training set minus (-) it's sub test fold
        test_fold = Data.matrix[x,] # here we describe the test fold individually
        
        
        # inner cross validation:
        
        #-------------#
        
        
        feature.selection <- matrix(nrow = 0, ncol = 2)
        
        for(feature.id in 1:(ncol(training_fold) -1)){
          
          R <- training_fold[training_fold$Outcome == '0', feature.id]
          NR <- training_fold[training_fold$Outcome == '1', feature.id]
          
          P.value <- wilcox.test(R, NR)$p.value
          
          feature.selection <- data.frame(rbind(feature.selection, cbind(feature.id, P.value)))
        }
        
        # multiple comparison correction
        feature.selection$P.value <- p.adjust(feature.selection$P.value, method = 'BH')
        
        # find all features which P value < 0.05
        feature.selection <- feature.selection[feature.selection$P.value <= 0.05,]
        
        feature.selection <- feature.selection[complete.cases(feature.selection),]
        
        # get id
        select.feature.id <- as.numeric(feature.selection$feature.id)
        
        if(length(select.feature.id) != 0){
          # selected feature matrix
          feature.selected.test <- data.frame(training_fold[,select.feature.id])
          
          # get name
          names.vec <- colnames(training_fold)[select.feature.id]
          
          colnames(feature.selected.test) <- names.vec
          
          # assign outcome
          feature.selected.test$Outcome <- as.factor(training_fold$Outcome)
          
          # now apply (train) the classifer on the training_fold
          optimized.feature[id] <- list(names.vec)
          
          
          classifier = svm(formula = Outcome ~ .,
                           data = feature.selected.test,
                           type = 'C-classification',
                           kernel = kernel,
                           gamma = gamma,
                           cost = cost, probability = TRUE)      
          # next step in the loop, we calculate the predictions and cm and we equate the accuracy
          # note we are training on training_fold and testing its accuracy on the test_fold
          y_pred = predict(classifier, newdata = test_fold[names.vec],  probability = TRUE)
          
          
          
          
          cm = table(pred = y_pred, true=test_fold$Outcome)
          accuracy = sum(diag(cm))/sum(cm)
          
          Accuracy <- rbind(Accuracy, accuracy)
          
          #inner.cv <- rbind(inner.cv, cbind(gamma, cost, accuracy))

          
          #--------- this section for logistic regression -----------#
          #pred_label <- data.frame(as.character(y_pred))
          #names(pred_label) <- 'Predicted'
          
          
          #test_fold_forGLM <- cbind(test_fold[, -ncol(test_fold)], pred_label)
          
          #names(test_fold_forGLM)[ncol(test_fold_forGLM)] <- 'Predicted.Outcome'
          #glm.fit <- glm(Predicted ~., data = test_fold_forGLM, family = binomial)
          
          #glm.fit_list[id] <- list(glm.fit)
          

          
          #optimization <- rbind(optimization, inner.cv[which.max(inner.cv[,3]),])
          
        }
      }
      Accuracy.folds <- mean(Accuracy)
      inner.cv <- rbind(inner.cv, cbind(gamma, cost, Accuracy.folds))
      
      
      iterations <- iterations + 1
      
      print('------------------')
      print(paste('finish iteration ', iterations, sep = ''))
    }
  }
  #return(list(optimization, optimized.feature))
  return(inner.cv)
  
}



#--------------------------------------------------------------#
#----------------------- Random Forest  Classifier ----------------------#
#--------------------------------------------------------------#

# ------- Feature selection by t test ---------#



rfClassifier <- function(Data.matrix, fold, kernel){
  
  #Data.matrix <- allFeatures_forModel.spatstat
  #fold <- 10
  set.seed(100)
  inner.cv <- matrix(nrow = 0, ncol = 3)
  
  folds = createFolds(Data.matrix$Outcome, k = fold)
  # machine learning
  
  
  
  optimization<- matrix(nrow = 0, ncol = 3)
  
  optimized.feature <- list()    
  #glm.fit_list <- list()
  
  iterations <- 0
  
  
  for(maxnode in 2:10){
    for (nodesize in 2:35) {
      
      Accuracy <- matrix(nrow = 0, ncol = 1)
      
      for (id in 1:length(folds)) {
        #id <- 5
        x = folds[[id]]
        
        # in the next two lines we will separate the Training set into it's 10 pieces
        training_fold = Data.matrix[-x,] # training fold =  training set minus (-) it's sub test fold
        test_fold = Data.matrix[x,] # here we describe the test fold individually
        
        
        # inner cross validation:
        
        #-------------#
        
        
        feature.selection <- matrix(nrow = 0, ncol = 2)
        
        for(feature.id in 1:(ncol(training_fold) -1)){
          
          R <- training_fold[training_fold$Outcome == '0', feature.id]
          NR <- training_fold[training_fold$Outcome == '1', feature.id]
          
          P.value <- wilcox.test(R, NR)$p.value
          
          feature.selection <- data.frame(rbind(feature.selection, cbind(feature.id, P.value)))
        }
        
        # multiple comparison correction
        feature.selection$P.value <- p.adjust(feature.selection$P.value, method = 'fdr')
        
        # find all features which P value < 0.05
        feature.selection <- feature.selection[feature.selection$P.value <= 0.05,]
        
        feature.selection <- feature.selection[complete.cases(feature.selection),]
        
        # get id
        select.feature.id <- as.numeric(feature.selection$feature.id)
        
        if(length(select.feature.id) != 0){
          # selected feature matrix
          feature.selected <- data.frame(training_fold[,select.feature.id])
          
          # get name
          names.vec <- colnames(training_fold)[select.feature.id]
          
          colnames(feature.selected) <- names.vec
          
          # assign outcome
          feature.selected$Outcome <- as.factor(training_fold$Outcome)
          
          # now apply (train) the classifer on the training_fold
          optimized.feature[id] <- list(names.vec)
          
        
          classifier = randomForest(Outcome ~ .,
                                       data = feature.selected,
                                       maxnodes = maxnode,
                                       ntree = 100,
                                       nodesize = nodesize,
                                       importance = TRUE) 

                    
          # note we are training on training_fold and testing its accuracy on the test_fold
          y_pred = predict(classifier, newdata = test_fold[names.vec],  probability = TRUE)
          
          
          
          
          cm = table(pred = y_pred, true=test_fold$Outcome)
          accuracy = sum(diag(cm))/sum(cm)
          
          Accuracy <- rbind(Accuracy, accuracy)
          
          #inner.cv <- rbind(inner.cv, cbind(gamma, cost, accuracy))
          
          
        }
      }
      Accuracy.folds <- mean(Accuracy)
      inner.cv <- rbind(inner.cv, cbind(maxnode, nodesize, Accuracy.folds))
      
      
      iterations <- iterations + 1
      
      print('------------------')
      print(paste('finish iteration ', iterations, sep = ''))
    }
  }
  
  #return(list(optimization, optimized.feature))
  return(inner.cv)
  
}



#-------------------------------------------------------------------#
#----------------------- Calculate DoC score  ----------------------#
#-------------------------------------------------------------------#

# ------- Feature selection by t test ---------#



scoreDoC <- function(DF1, DF2,dThresh){
  
  DF1 <- Cancer
  DF2 <- lymphocytes
  

  rho.all <- matrix(nrow = 0, ncol = 1)
  
  for(pts in 1:nrow(DF1)){
    #pts <- 2
    Point <- DF1[pts, 1:2]
    
    # vector definition
    intra.vec <- matrix(nrow = 0, ncol = 1)
    inter.vec <- matrix(nrow = 0, ncol = 1)
    for (d in seq(1, 20, 0.5)) {
      #d <- 15
      # how many points within the range (inter)
      mtrx.inter <- nn2(DF2[,1:2], DF1[pts,1:2], k = nrow(DF2), treetype = 'kd', searchtype = 'radius', radius = d) 
      dens.inter <- rowSums(data.frame(mtrx.inter$nn.idx) != 0)/(pi*d^2)
      
      # how many points within the range (intra)
      mtrx.intra <- nn2(DF1[,1:2], DF1[pts,1:2], k = nrow(DF1), treetype = 'kd', searchtype = 'radius', radius = d) 
      dens.intra <- rowSums(data.frame(mtrx.intra$nn.idx) != 0)/(pi*d^2)
      
      # append density to vector
      intra.vec <- rbind(intra.vec, dens.intra)
      inter.vec <- rbind(inter.vec, dens.inter)
    }  
    rho <- cor.test(as.vector(intra.vec), as.vector(inter.vec), method = 'pearson')
    
    rho.all <- rbind(rho.all, rho$estimate)
    
  }
  
  
}




