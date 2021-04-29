##############################################
#------- Cell classification ----------------#
##############################################

# @Author: Haoyang Mi, JHU

library(ggplot2)

library(igraph)

library(RANN)

library(pracma)

library(ggtern)

library(RANN)

library(mixtools)

library(flexclust)

library(matrixStats)
library(dplyr)

setwd("~/Desktop/Validation-cohort")
source('./FunctionHub.R')



# read qualified cores
CORE.VALID <- read.csv('TMA_HE_Meta.csv')%>%
  filter(Missing == 'False') %>%
  select(Name)

#write.csv(CORE.VALID, '~/Desktop/Validation-cohort/CORE.VALID.csv')


#------- get core area ------#

CORE.Split.area <- read.csv('TMA_allArea.csv')
core.Area <- data.frame(matrix(nrow = 0, ncol = 0))
for(core in unique(CORE.Split.area$TMA.core)){
  
  #core <- 'G-2'
  area.split <- sum(CORE.Split.area[CORE.Split.area$TMA.core == core, 2])
  
  core.Area <- rbind(core.Area, cbind(core, area.split))
}
#write.csv(core.Area, 'TMA_allArea_combined.csv')

Cell.VALID.label <- read.csv('HE_Coordinates_Features.csv')

#---------------------------------------#
# Rescale H&E points for 963 -----------#

TMA.VALID.meta <- read.csv('TMA_HE_Meta.csv')

for(core in CORE.VALID$Name){
  
  
  
  # subset DF, convert to pixel unit
  ptsFile <- Cell.VALID.label[Cell.VALID.label$TMA.core == core, c(2:4)] 
  
  # get ref points
  
  ref_x <- TMA.VALID.meta[TMA.VALID.meta$Name == core,]$Centroid.X.µm - 4340*0.2305/2 # unit: pixel
  ref_y <- TMA.VALID.meta[TMA.VALID.meta$Name == core,]$Centroid.Y.µm - 4340*0.2305/2 # unit: pixel
  
  
  # get the rescaled points
  rescale_x <- ptsFile$Centroid.X.µm - ref_x
  rescale_y <- ptsFile$Centroid.Y.µm - ref_y  
  
  # get type
  Cell.VALID.label[Cell.VALID.label$TMA.core == core, 3:4] <- cbind(rescale_x, rescale_y) 
  
}

#test <- Cell.VALID.label[Cell.VALID.label$TMA.core == 'G-5',]

#plot(test$Centroid.X.µm, test$Centroid.Y.µm)

# prepare tables for python
#colnames(Cell.VALID.label)[3:4] <- c('X.um', 'Y.um')
#write.csv(Cell.VALID.label, '~/Desktop//TMA_DoC/Cell.VALID.label.csv', row.names = FALSE)
#write.csv(Cell.963.label, '~/Desktop/TMA_DoC/Cell.963.label.csv')

# read cell feature files




#Cell.VALID.label <- read.csv('~/Desktop/TMA_DoC/Cell.1042.label.csv')
#Cell.963.label <- read.csv('~/Desktop/TMA_DoC/Cell.963.label.csv')

#-------------------------------------------------------------#
#------------------- Computation starts ----------------------#
#-------------------------------------------------------------#


#TMA.set <- '963'
CORE <- CORE.VALID

CORE.area <- core.Area

Cell.label <- Cell.VALID.label

Region.path <- '~/Desktop/Validation-cohort/TMA_Boundary_Rescaled/'

test <- readRDS('~/Desktop/Validation-cohort/TMA_Boundary_Rescaled/A-2/HE.rds')
test <- lapply(test, '*', 0.2305) ## My old clunky way


points(Cell.label[Cell.label$TMA.core == 'A-2', 3:4])

Mutual.dist <- data.frame(matrix(nrow = 0, ncol = 50))

for(core in CORE$Name){
  
  #core <- 'A-1'
  core.area <- as.numeric(as.character(CORE.area[CORE.area$core == core, ]$area.split))
  
  all.cell <- Cell.label[Cell.label$TMA.core == core,]
  
  # immune cell subset
  lymphocytes <- all.cell[all.cell$Class == 'Immune cells', 3:5]
  
  # cell distance
  Cancer <- all.cell[all.cell$Class == 'Tumor', 3:5]
  
  # stromal
  Stromal <- all.cell[all.cell$Class == 'Stroma', 3:5]
  
  
  # ratios
  # density
  #-------------------------------------------------------------#
  lymphocyte.dens <- nrow(lymphocytes)/core.area
  cancer.dens <- nrow(Cancer)/core.area
  stroma.dens <- nrow(Stromal)/core.area
  
  # how many cancer cells within 20 microns for each lymphocytes?
  #-------------------------------------------------------------#
  cancer.to.lympho.matrix <- nn2(Cancer[,1:2], lymphocytes[,1:2], k = nrow(Cancer), treetype = 'kd', searchtype = 'radius', radius = 20)
  
  cancer.to.lympho.avgCount <- mean(rowSums(data.frame(cancer.to.lympho.matrix$nn.idx) != 0))
  
  
  
  # how many lymphocytes cells within 20 microns for each cancer cell?
  #-------------------------------------------------------------#
  lympho.to.cancer.matrix <- nn2(lymphocytes[,1:2], Cancer[,1:2], k = nrow(lymphocytes), treetype = 'kd', searchtype = 'radius', radius = 20)
  
  lympho.to.cancer.avgCount <- mean(rowSums(data.frame(lympho.to.cancer.matrix$nn.idx) != 0))
  
  
  
  # how many cancer cells within 20 microns for each stromal cells?
  #-------------------------------------------------------------#
  cancer.to.stroma.matrix <- nn2(Cancer[,1:2], Stromal[,1:2], k = nrow(Cancer), treetype = 'kd', searchtype = 'radius', radius = 20)
  
  cancer.to.stroma.avgCount <- mean(rowSums(data.frame(cancer.to.stroma.matrix$nn.idx) != 0))
  
  
  
  # how many stromal cells within 20 microns for each cancer cell?
  #-------------------------------------------------------------#
  stroma.to.cancer.matrix <- nn2(Stromal[,1:2], Cancer[,1:2], k = nrow(Stromal), treetype = 'kd', searchtype = 'radius', radius = 20)
  
  stroma.to.cancer.avgCount <- mean(rowSums(data.frame(stroma.to.cancer.matrix$nn.idx) != 0))
  
  
  
  # how many cancer cells within 20 microns for each lymphocytes?
  #-------------------------------------------------------------#
  stroma.to.lympho.matrix <- nn2(Stromal[,1:2], lymphocytes[,1:2], k = nrow(Stromal), treetype = 'kd', searchtype = 'radius', radius = 20)
  
  stroma.to.lympho.avgCount <- mean(rowSums(data.frame(stroma.to.lympho.matrix$nn.idx) != 0))
  
  
  # how many lymphocytes cells within 20 microns for each cancer cell?
  #-------------------------------------------------------------#
  lympho.to.stroma.matrix <- nn2(lymphocytes[,1:2], Stromal[,1:2], k = nrow(lymphocytes), treetype = 'kd', searchtype = 'radius', radius = 20)
  
  lympho.to.stroma.avgCount <- mean(rowSums(data.frame(lympho.to.stroma.matrix$nn.idx) != 0))
  
  
  
  
  # spatial G-cross function for 
  
  Region_HE <- readRDS(paste(Region.path, core, '/HE.rds', sep = ''))
  #Region_HE <- lapply(Region_HE, '*', 0.454) ## My old clunky way
  
  plot(Region_HE[[1]])
  points(Cancer[,1:2])
  
  list1 <- bivarAnalysis.Kcross('Cancer', 'lymphocytes', Region_HE)
  list2 <- bivarAnalysis.Kcross('Cancer', 'Stromal', Region_HE)
  list3 <- bivarAnalysis.Kcross('lymphocytes', 'Stromal', Region_HE)
  
  #---------------------------------------------------------------#
  #-----------------------DoC score stats-------------------------#
  #---------------------------------------------------------------#
  
  DoC_lympho <- read.csv(paste('~/Desktop/TMA_DoC/TMA-Validation_DoC/', core, '/lympho.csv', sep = ''))
  colnames(DoC_lympho) <- c('x.um', 'y.um', 'Tumor_Lympho', 'Stroma_Lympho')
  
  DoC_Stroma <- read.csv(paste('~/Desktop/TMA_DoC/TMA-Validation_DoC/', core, '/stroma.csv',  sep = ''))
  colnames(DoC_Stroma) <- c('x.um', 'y.um', 'Tumor_Stroma', 'Lympho_Stroma')
  
  DoC_Tumor <- read.csv(paste('~/Desktop/TMA_DoC/TMA-Validation_DoC/', core, '/tumor.csv',  sep = ''))
  colnames(DoC_Tumor) <- c('x.um', 'y.um', 'Lympho_Tumor', 'Stroma_Tumor')
  
  
  DoC_L_max <- t(colMaxs(as.matrix(DoC_lympho[,3:4]), na.rm = T))
  DoC_L_min <- t(colMins(as.matrix(DoC_lympho[,3:4]), na.rm = T))
  DoC_L_mean <- t(colMeans(as.matrix(DoC_lympho[,3:4]), na.rm = T))
  DoC_L_sd <- t(colSds(as.matrix(DoC_lympho[,3:4]), na.rm = T))
  
  DoC_T_max <- t(colMaxs(as.matrix(DoC_Tumor[,3:4]), na.rm = T))
  DoC_T_min <- t(colMins(as.matrix(DoC_Tumor[,3:4]), na.rm = T))
  DoC_T_mean <- t(colMeans(as.matrix(DoC_Tumor[,3:4]), na.rm = T))
  DoC_T_sd <- t(colSds(as.matrix(DoC_Tumor[,3:4]), na.rm = T))
  
  DoC_S_max <- t(colMaxs(as.matrix(DoC_Stroma[,3:4]), na.rm = T))
  DoC_S_min <- t(colMins(as.matrix(DoC_Stroma[,3:4]), na.rm = T))
  DoC_S_mean <- t(colMeans(as.matrix(DoC_Stroma[,3:4]), na.rm = T))
  DoC_S_sd <- t(colSds(as.matrix(DoC_Stroma[,3:4]), na.rm = T))
  
  localized_TL <- as.numeric(as.character(nrow(DoC_lympho[DoC_lympho$Tumor_Lympho >= 0.847, ])/nrow(DoC_lympho)))
  
  localized_SL <- as.numeric(as.character(nrow(DoC_lympho[DoC_lympho$Stroma_Lympho >= 0.847, ])/nrow(DoC_lympho)))
  
  localized_LT <- as.numeric(as.character(nrow(DoC_Tumor[DoC_Tumor$Lympho_Tumor >= 0.847, ])/nrow(DoC_Tumor)))
  
  localized_ST <- as.numeric(as.character(nrow(DoC_Tumor[DoC_Tumor$Stroma_Tumor >= 0.847, ])/nrow(DoC_Tumor)))
  
  localized_TS <- as.numeric(as.character(nrow(DoC_Stroma[DoC_Stroma$Tumor_Stroma >= 0.847, ])/nrow(DoC_Stroma)))
  
  localized_LS <- as.numeric(as.character(nrow(DoC_Stroma[DoC_Stroma$Lympho_Stroma >= 0.847, ])/nrow(DoC_Stroma)))
  
  
  
  DoC_features <- cbind(DoC_L_max, DoC_L_min, DoC_L_mean, DoC_L_sd, DoC_T_max, DoC_T_min, DoC_T_mean, DoC_T_sd, DoC_S_max, DoC_S_min, DoC_S_mean, DoC_S_sd, localized_TL, localized_SL, localized_LT, localized_ST, localized_TS, localized_LS)
  colnames(DoC_features) <- c('DoC_TL_max', 'DoC_SL_max', 'DoC_TL_min', 'DoC_SL_min', 'DoC_TL_mean', 'DoC_SL_mean', 'DoC_TL_Std', 'DoC_SL_Std', 
                              'DoC_LT_max', 'DoC_ST_max', 'DoC_LT_min', 'DoC_ST_min', 'DoC_LT_mean', 'DoC_ST_mean', 'DoC_LT_Std', 'DoC_ST_Std', 
                              'DoC_TS_max', 'DoC_LS_max', 'DoC_TS_min', 'DoC_LS_min', 'DoC_TS_mean', 'DoC_LS_mean', 'DoC_TS_Std', 'DoC_LS_Std', 
                              'localized_TL', 'localized_SL', 'localized_LT', 'localized_ST', 'localized_TS', 'localized_LS')
  #-------------------------------------------------------------------------#
  #----------------------- spatial Shannon entropy -------------------------#
  #-------------------------------------------------------------------------#
  
  
  # lymphocytes variable name: lymphocytes
  
  # cancer variable name: Cancer
  
  # stromal variable name: Stromal
  
  
  lympho_sub <- lymphocytes
  cancer_sub <- Cancer
  stromal_sub <- Stromal
  
  
  type.count <- 3 # number of types
  
  
  
  Total <- nrow(stromal_sub) + nrow(cancer_sub) + nrow(lympho_sub) # total number of cells 
  
  
  #-------------- cancerous  --------------#
  
  cancer_int <- 0
  cancer_ext <- 0
  p <- 0
  
  if(isTRUE(nrow(cancer_sub) != 0)){
    p <- nrow(cancer_sub)/Total
    cancer_int <- mean(as.matrix(dist(cancer_sub[,2:3]))) # remove the id column
    
    if(isTRUE(nrow(stromal_sub)!=0)){
      cancer_ext <- cancer_ext + mean(dist2(cancer_sub[,2:3], stromal_sub[,2:3]))
    }
    
    if(isTRUE(nrow(lympho_sub)!=0)){
      cancer_ext <- cancer_ext + mean(dist2(cancer_sub[,2:3], lympho_sub[,2:3]))
    }
    
  }
  cancer_dat <- cbind(p, cancer_int, cancer_ext/(type.count - 1))
  colnames(cancer_dat) <- c('p', 'int', 'ext')
  
  
  #-------------- lymphocytes --------------#
  lympho_int <- 0
  lympho_ext <- 0
  p <- 0
  if(isTRUE(nrow(lympho_sub) != 0)){
    p <- nrow(lympho_sub)/Total
    lympho_int <- mean(as.matrix(dist(lympho_sub[,2:3])))
    
    if(isTRUE(nrow(cancer_sub)!=0)){
      lympho_ext <- lympho_ext + mean(dist2(lympho_sub[,2:3], cancer_sub[,2:3]))
    }
    
    if(isTRUE(nrow(stromal_sub)!=0)){
      lympho_ext <- lympho_ext + mean(dist2(lympho_sub[,2:3], stromal_sub[,2:3]))
    }
  }
  lympho_dat <- cbind(p, lympho_int, lympho_ext/(type.count - 1))
  colnames(lympho_dat) <- c('p', 'int', 'ext')
  
  
  #-------------- stromal --------------#
  
  stromal_int <- 0
  stromal_ext <- 0
  p <- 0
  if(isTRUE(nrow(stromal_sub) != 0)){
    p <- nrow(stromal_sub)/Total
    stromal_int <- mean(as.matrix(dist(stromal_sub[,2:3])))
    
    if(isTRUE(nrow(cancer_sub)!=0)){
      stromal_ext <- stromal_ext + mean(dist2(stromal_sub[,2:3], cancer_sub[,2:3]))
    }
    
    if(isTRUE(nrow(lympho_sub)!=0)){
      stromal_ext <- stromal_ext + mean(dist2(stromal_sub[,2:3], lympho_sub[,2:3]))
    }
    
  }
  stromal_dat <- cbind(p, stromal_int, stromal_ext/(type.count - 1))
  colnames(stromal_dat) <- c('p', 'int', 'ext')
  
  
  Dat_collect <- rbind(cancer_dat, lympho_dat, stromal_dat)
  
  ShannonH <- 0
  # combine row data
  
  
  for(dat in seq(1, type.count)){
    p <- Dat_collect[dat,1]
    d_int <- Dat_collect[dat, 2]
    d_ext <- Dat_collect[dat, 3]
    d_final <- d_int/d_ext
    
    if(isTRUE(d_int*d_ext == 0)){
      d_final <- 0
    }
    if(isTRUE(p != 0)){
      ShannonH <- -d_final*p*log2(p) + ShannonH
    }
  }
  
  
  
  cancer.ratio <- nrow(cancer_sub)/Total     # for ternary plot
  
  stromal.ratio <- nrow(stromal_sub)/Total
  
  lympho.ratio <- nrow(lympho_sub)/Total
  
  
  #---- combine features ----#
  
  
  # composite
  #Mutual.dist <- rbind(Mutual.dist, cbind(core, TMA.set, lymphocyte.dens, cancer.dens, stroma.dens, cancer.to.lympho.avgCount, lympho.to.cancer.avgCount, cancer.to.stroma.avgCount, stroma.to.cancer.avgCount, stroma.to.lympho.avgCount, lympho.to.stroma.avgCount))
  # composite
  Mutual.dist <- rbind(Mutual.dist, cbind(core, list1[[1]], list1[[2]],  list2[[1]], list2[[2]], list3[[1]], list3[[2]],cancer.ratio, stromal.ratio, lympho.ratio,
                                          lymphocyte.dens, cancer.dens, stroma.dens, cancer.to.lympho.avgCount, lympho.to.cancer.avgCount, cancer.to.stroma.avgCount, stroma.to.cancer.avgCount, stroma.to.lympho.avgCount, lympho.to.stroma.avgCount, ShannonH, DoC_features))
  
  print(core)
  
}

colnames(Mutual.dist)[2:10] <- c('Kcross.Cancer2lympho', 'Kcross.lympho2Cancer', 'Kcross.Cancer2Stromal', 'Kcross.Stromal2Cancer', 'Kcross.Lympho2Stromal', 'Kcross.Stromal2Lympho', 'cancer ratio', 'stromal ratio', 'lymphocytes ratio')


write.csv(Mutual.dist, 'classification_features.csv')



#Mutual.dist$ShannonH <- as.numeric(as.character(Mutual.dist$ShannonH))






#------------------------------------------#
#------------- Ternary plot  --------------#
#------------------------------------------#

# read ternary data
ternary.dat <- Mutual.dist[c('core', 'cancer ratio', 'stromal ratio', 'lymphocytes ratio')]
colnames(ternary.dat)[1] <- 'TMA.core'


# convert factor data to numerical data
ternary.dat$`cancer ratio` <- as.numeric(as.character(ternary.dat$`cancer ratio`))
ternary.dat$`stromal ratio` <- as.numeric(as.character(ternary.dat$`stromal ratio`))
ternary.dat$`lymphocytes ratio` <- as.numeric(as.character(ternary.dat$`lymphocytes ratio`))


# split data to 963 and 1042 and then merge with outcome
ternary.dat.963 <- ternary.dat[ternary.dat$TMA.set == '963',]
ternary.dat.1042 <- ternary.dat[ternary.dat$TMA.set == '1042',]

ternary.dat.963 <- merge(ternary.dat.963, TMA.963.Outcome, by = 'TMA.core')
ternary.dat.1042 <- merge(ternary.dat.1042, TMA.1042.Outcome, by = 'TMA.core')

# recombine
ternary.dat.withOutcome <- rbind(ternary.dat.963, ternary.dat.1042)
ternary.dat.withOutcome[, 3:5] <- ternary.dat.withOutcome[, 3:5]*100



# plot
jpeg('./cellClass_triplot.jpeg', units="in", width=5, height=5, res=300)

ggtern(ternary.dat.withOutcome, aes(x = `cancer ratio`, y = `stromal ratio`, z = `lymphocytes ratio`))+
  # plot points
  geom_point(mapping = aes(fill = as.factor(Outcome), size = 21, shape = as.factor(TMA.set)), size =4)+
  # set shape
  scale_shape_manual(values = c(21, 24)) +
  # background style
  theme_bw()+
  # show arrow
  theme_showarrows() +
  # color scheme for responder/non-responder
  scale_fill_manual(values = c('#a1ffc5', '#ff8a7a'))+
  # customize theme
  theme(tern.axis.line.T = element_line(color='#7a97ff',size=0.5), # ternary axis color
        tern.axis.line.R = element_line(color='#b067ae',size=0.5),
        tern.axis.line.L = element_line(color='#f77c47',size=0.5),
        tern.panel.grid.major.T = element_line(linetype = 'twodash', color= '#7a97ff'), # ternary grid color and line style
        tern.panel.grid.major.R = element_line(linetype = 'twodash', color= '#b067ae'),
        tern.panel.grid.major.L = element_line(linetype = 'twodash', color= '#f77c47'),
        tern.axis.arrow.T       = element_line(color = '#7a97ff', size = 2), # ternary axis arrow color and size
        tern.axis.arrow.R       = element_line(color = '#b067ae', size = 2),
        tern.axis.arrow.L       = element_line(color = '#f77c47', size = 2),
        tern.axis.arrow.text.T       = element_text(color = '#7a97ff', size = 12),# ternary axis arrow text color and size
        tern.axis.arrow.text.R       = element_text(color = '#b067ae', size = 12),
        tern.axis.arrow.text.L       = element_text(color = '#f77c47', size = 12),
        tern.axis.text = element_text(size = 15))+ # ternary axis text size
  # label 
  labs( x = '',
        xarrow  = "Cancerous cell ratio (%)",
        y = '',
        yarrow  = "Stromal cell ratio (%)",
        z = '',
        zarrow  = "Lymphocytes ratio (%)") +
  # remove legends
  guides(fill = 'none', color = 'none', shape = 'none')

dev.off()



#---------------------------------------------------------------------------#
#------------- This section preliminarily test the importance --------------#
#---------------------------------------------------------------------------#

TMA.963.Outcome <- TMA.963.meta[complete.cases(TMA.963.meta), ]
colnames(TMA.963.Outcome)[1] <- 'core'

TMA.1042.Outcome <- TMA.1042.meta[TMA.1042.meta$Outcome != '', ]
colnames(TMA.1042.Outcome)[1] <- 'core'



Mutual.dist.963 <- Mutual.dist[Mutual.dist$TMA.set == '963',]

Mutual.dist.1042 <- Mutual.dist[Mutual.dist$TMA.set == '1042',]


test1 <- merge(TMA.963.Outcome, Mutual.dist.963, by = c('core'))


test2 <- merge(TMA.1042.Outcome, Mutual.dist.1042, by = c('core'))

test3 <- rbind(test1, test2)


#-------------------------------------------------------------------------#
#--------------------------Gaussian mixture model-------------------------#
#-------------------------------------------------------------------------#

GMM <- normalmixEM(test3$ShannonH, k =2)
plot(GMM, which = 2)

#---- get GMM parameters -----#
mu <- GMM[['mu']] #mu
sigma <- GMM[['sigma']] #sigma
lambda <- GMM[['lambda']] # amplitudes



jpeg('./histo_spatialproxi.jpeg', units="in", width=5.5, height=5, res=300)

ggplot(test3, aes(x = ShannonH), color = Outcome) +
  geom_histogram(binwidth = 0.05, fill = NA, color = Outcome) +
  mapply(
    function(mean, sd, lambda, n, binwidth) {
      stat_function(
        fun = function(x) {
          (dnorm(x, mean = mean, sd = sd)) * n * binwidth * lambda
        }
      )
    },
    mean = mu, #mean
    sd = sigma, #standard deviation
    lambda = lambda, #amplitude
    n = length(test3$ShannonH), #sample size
    binwidth = 0.05 #binwidth used for histogram
  )

dev.off()   



write.csv(test3, 'Figure6_plot.csv')


for(colid in 6:26){
  
  
  NR <- as.numeric(as.character(test3[test3$Outcome == '1', ]$ShannonH))
  R <- as.numeric(as.character(test3[test3$Outcome == '0', ]$ShannonH))
  
  print(wilcox.test(NR, R))
  
}

#-------------------------------------------------------#
#------------ This section plot Kcross -----------------#
#-------------------------------------------------------#
Region <- Region_HE
pts.type1 <- 'lymphocytes'
pts.type2 <- 'Cancer'

type1 <- get(eval(pts.type1))[,1:2]
colnames(type1) <- c('x', 'y')

type2 <- get(eval(pts.type2))[,1:2]
colnames(type2) <- c('x', 'y')



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





multitype_ppp <- ppp(pts_OI$x, pts_OI$y, marks = species, owin(poly = Region))
K.cross <- data.frame(Kcross(multitype_ppp, i = pts.type1, j = pts.type2, r = seq(0,20,0.1), correction = 'Ripley'))

#plot(K.cross)
# relocat DF

K.cross <- K.cross[complete.cases(K.cross),]


value.new <- K.cross$iso - K.cross$theo

K.cross.theo <- K.cross[,1:2]
colnames(K.cross.theo) <- c('r', 'value')

K.cross.iso <- K.cross[,c(1,3)]
colnames(K.cross.iso) <- c('r', 'value')

K.cross.theo <- K.cross.theo %>% arrange(desc(row_number()))


kcross.reshape <- rbind(K.cross.theo, rev(K.cross.iso))

jpeg('./Kcross_cancer_lympho.jpeg', units="in", width=5.5, height=5, res=300)

ggplot(data = K.cross) +
  theme_bw()+
  #geom_point(aes(r, theo), fill = '#3894f0', shape = 21, color = 'black', size = 4)+
  #geom_point(aes(r, iso), fill = '#fd7894', shape = 21, color = 'black', size = 4)+
  geom_line(aes(r, iso), color = '#fd7894', size = 2)+
  geom_line(aes(r, theo), color = '#3894f0', size = 2)+
  geom_polygon(data = kcross.reshape, aes(r, value, alpha = 4), fill = 'grey')+
  
  xlab(expression(paste('distance, ', mu, 'm', sep = ''))) +
  ylab(expression(paste('K'[L~','~C, sep = ''], '(r)', sep = ''))) +
  #xlim(0, 20) +
  theme(axis.title = element_text(size = 36),
        axis.title.y = element_text(size = 35),
        axis.text = element_text(size = 25),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1),
        legend.position = 'none') +
  guides(fill = NA)

dev.off()   




#--------------------------------------------#
#------------- histogram plot  --------------#
#--------------------------------------------#
histo.cancer.to.lympho <- cbind(as.numeric(as.character(Mutual.dist$lympho.to.cancer.avgCount)), 'c.to.l')
histo.cancer.to.stromal <- cbind(as.numeric(as.character(Mutual.dist$cancer.to.stroma.avgCount)), 'c.to.s')
histo.lympho.to.stromal <- cbind(as.numeric(as.character(Mutual.dist$cancer.to.stroma.avgCount)), 'l.to.s')


histo.together <- data.frame(rbind(histo.cancer.to.lympho, histo.cancer.to.stromal, histo.lympho.to.stromal))
histo.together$X1 <- as.numeric(as.character(histo.together$X1))
jpeg('./histo_spatialproxi.jpeg', units="in", width=5.5, height=5, res=300)

ggplot(data = histo.together, aes(x = X1, color = X2)) +
  theme_bw()+
  geom_histogram(fill = NA, size = 2)+
  xlab('Cell counts within ROI') +
  ylab('Frequency')+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  
  theme(axis.title = element_text(size = 31),
        axis.title.y = element_text(size = 31),
        axis.text = element_text(size = 30),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1),
        legend.position = 'none') +
  guides(fill = NA)

dev.off()   



#--------------------------------------------#
#------------- DoC plot  --------------#
#--------------------------------------------#
TMA.set <- '1042'
core <- 'L-13'

DoC_lympho <- read.csv(paste('~/Desktop/TMA_DoC/TMA', TMA.set, '_DoC/', core, '/lympho.csv', sep = ''))
colnames(DoC_lympho) <- c('x.um', 'y.um', 'Tumor_Lympho', 'Stroma_Lympho')

DoC_Stroma <- read.csv(paste('~/Desktop/TMA_DoC/TMA', TMA.set, '_DoC/', core, '/stroma.csv',  sep = ''))
colnames(DoC_Stroma) <- c('x.um', 'y.um', 'Tumor_Stroma', 'Lympho_Stroma')

DoC_Tumor <- read.csv(paste('~/Desktop/TMA_DoC/TMA', TMA.set, '_DoC/', core, '/tumor.csv',  sep = ''))
colnames(DoC_Tumor) <- c('x.um', 'y.um', 'Lympho_Tumor', 'Stroma_Tumor')


Tumor_Lympho <- DoC_lympho[,1:3]
Tumor_Lympho <- Tumor_Lympho[complete.cases(Tumor_Lympho),]


Lympho_Tumor <- DoC_Tumor[,1:3]
Lympho_Tumor <- Lympho_Tumor[complete.cases(Lympho_Tumor),]

jpeg('./L-13_DoC.jpeg', units="in", width=7, height=5.2, res=300)
ggplot() +
  theme_bw()+
  #geom_point(data = Tumor_Lympho, aes(x.um, y.um, color = Tumor_Lympho), size = 2) +
  geom_point(data = Lympho_Tumor, aes(x.um, y.um, fill = Lympho_Tumor), size = 4, shape = 21) +
  ylim(1400, 0) +
  xlim(0, 1400) +
  theme(axis.title = element_text(size = 31),
        axis.title.y = element_text(size = 31),
        axis.text = element_text(size = 30),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.position = 'left') +
  scale_fill_gradient2(low = "blue", mid = "white",
                        high = "red", space = "Lab" )+
  guides(fill = guide_colourbar(barwidth = 1, barheight = 20)) +
  labs(fill = "DoC score")
dev.off()
  #geom_point(data = Lympho_Tumor, aes(x.um, y.um), size = 3, fill = Lympho_Tumor) 
  #geom_point(data = DoC_Stroma, aes(x.um, y.um), size = 3, fill = '#4379b9') 
  

#-------------------------------------------------------#
#------------- mean counts with SD  --------------------#
#-------------------------------------------------------#

Cell.1042.label$TMA.core <- as.character(Cell.1042.label$TMA.core)
select.cellTypeDat.1042 <- Cell.1042.label[Cell.1042.label$TMA.core %in% CORE.1042$x,]

n.cancer.all <-matrix(nrow = 0, ncol = 1)
n.lympho.all <-matrix(nrow = 0, ncol = 1)
n.stroma.all <-matrix(nrow = 0, ncol = 1)
for(core in CORE$x){
  
  n.cancer <- nrow(select.cellTypeDat.1042[select.cellTypeDat.1042$Class == 'Tumor' & select.cellTypeDat.1042$TMA.core == core,])
  n.lympho <- nrow(select.cellTypeDat.1042[select.cellTypeDat.1042$Class == 'Immune cells' & select.cellTypeDat.1042$TMA.core == core,])
  n.stroma <- nrow(select.cellTypeDat.1042[select.cellTypeDat.1042$Class == 'Stroma' & select.cellTypeDat.1042$TMA.core == core,])
  
  n.cancer.all <- rbind(n.cancer.all, n.cancer)
  n.stroma.all <- rbind(n.stroma.all, n.stroma)
  n.lympho.all <- rbind(n.lympho.all, n.lympho)
  
  
}
mean(n.lympho.all)
sd(n.lympho.all)

mean(n.cancer.all)
sd(n.cancer.all)

mean(n.stroma.all)
sd(n.stroma.all)


