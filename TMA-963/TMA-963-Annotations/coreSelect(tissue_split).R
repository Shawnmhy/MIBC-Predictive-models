###########################################
## This script is to select valid cores ###
###########################################
library(dplyr)
library(stats)
library(reshape2)
library(gridExtra)
library(cowplot) # 1.0.0

library(scales)
library(stringr)
#pca 
library(ggfortify)
library(ggbiplot)
# t-SNE
library(Rtsne)
library(largeVis)
library(cluster)

set.seed(1)
# set directory
setwd("~/Desktop/TMA_963/TMA-963-Annotations")


#-------------- read all tissue stats -----------------#

# step 1: read csv file
allStats <- read.csv('allTMA_tissueArea.csv')

# step 2: split the image name for matching
allStats$Image <- str_split_fixed(allStats$Image, "_", 6)[,4]


#-------------- t-SNE dat -----------------#
t_SNE_dat <- matrix(nrow = 0, ncol = 2)
colnames(t_SNE_dat) <- c('stain_name', 'V2')

new_stain_dat <- matrix(nrow = 0, ncol = 2)

for(flag in 1:10){ # num of stainings
  # list of stain vectors
  stain_name <- switch (flag, 'HE', 'P16', 'P53', 'P63', 'uv-GATA3', 'Ki67', 'CK20', 'CK5-6', 'Her2Neu', 'CyclinD1')
  
  stain_dat <- allStats[allStats$Image == stain_name,]
  
  new_stain_dat <- matrix(nrow = 0, ncol = 2)
  for(letter in LETTERS[1:12]){
    for(num in 1:10){
      
      # core name     
      coreName <- paste(letter, '-', num, sep = '')
      
      # get the core dat for curremt stain
      stain_core_dat <- stain_dat[stain_dat$TMA.core == coreName,]
      
      if(nrow(stain_core_dat) != 0){
        coreArea <- sum(stain_core_dat$Area.µm.2)*0.454*0.454/10^6
        sum_dat <- cbind(coreName, coreArea)
        colnames(sum_dat) <- c('TMA.core', 'Area')
        
        new_stain_dat <- data.frame(rbind(new_stain_dat, sum_dat))
        
        # remove the core dat from original df
        #stain_dat <- anti_join(stain_dat, stain_core_dat)
      }

      colnames(new_stain_dat) <- c('TMA.core', 'Area')
      new_stain_dat <- unique(new_stain_dat)
    }

    
  }
  if(stain_name == 'HE'){
    Merge_stats <- new_stain_dat
  } else {
    Merge_stats <- merge(Merge_stats, new_stain_dat, by = 'TMA.core')  
  }
  
  
  
  # prepare t-SNE data
  t_SNE_dat <- data.frame(rbind(t_SNE_dat, cbind(stain_name, as.numeric(as.character(new_stain_dat$Area)))))
  
}



Merge_stats[,2:11] <- as.matrix(Merge_stats[,2:11])

# convert unit: um^2 to mm^2

CV <- data.frame(matrix(nrow = 0, ncol = 1))
Mean <- data.frame(matrix(nrow = 0, ncol = 1))
QCoD <- data.frame(matrix(nrow = 0, ncol = 1))

for(i in 1:nrow(Merge_stats)){ # num of stainings
  vec <- as.numeric(as.character(Merge_stats[i, 2:11]))
  CV <- rbind(CV, sd(vec)/mean(vec))
  Mean <- rbind(Mean, mean(vec))
  QCoD <-  rbind(QCoD, (quantile(vec, 0.75) - quantile(vec, 0.25))/(quantile(vec, 0.75) + quantile(vec, 0.25)))
}

Merge_stats <- data.frame(cbind(Merge_stats, Mean, CV, QCoD))
colnames(Merge_stats) <- c('TMA.core', 'H&E', 'P16', 'P53', 'P63', 'uv-GATA3', 'Ki67', 'CK20', 'CK5-6', 'Her2Neu', 'Cyclin', 'Mean', 'CoV', 'QCoD')
Merge_stats$TMA.core <- as.factor(Merge_stats$TMA.core)
# reorder data frame to plot

Merge_stats <- Merge_stats[order(Merge_stats$Mean),]
Merge_stas_mean_order <- Merge_stats[,c(1, 12:14)]
# sort the dataframe for boxplot
Merge_stats_melt <- melt(data = Merge_stats[, c(1:11)], id.vars = 'TMA.core', measure.vars =  c('H&E', 'P16', 'P53', 'P63', 'uv-GATA3', 'Ki67', 'CK20', 'CK5-6', 'Her2Neu', 'Cyclin'))

Merge_stats_melt$value <- as.numeric(Merge_stats_melt$value)
##### seclect cores

#Cores_selected <- Merge_stats[Merge_stats$CoV <- 0.1,]

###############################
## This section is for plot ###
###############################


########################
######## t-SNE #########
########################


tsne_df = setNames(data.frame(t(Merge_stats[,-c(1,12,13,14)])), Merge_stats[,1])

tsne_df <- as.data.frame(lapply(tsne_df, as.numeric)) #<- sapply is here

tsne_model_1 = Rtsne(tsne_df, check_duplicates=FALSE, pca=TRUE, perplexity=1, theta=0.5, dims=2)

d_tsne_1 = as.data.frame(tsne_model_1$Y)  

## keeping original data
d_tsne_1_original <- d_tsne_1

# score vector

Score_df <- matrix(nrow = 0, ncol = 1)

for(i in 1:(nrow(d_tsne_1) - 1)){
  
  # create grid for parameter optimization
  
  # generate clusters
  fit_cluster_kmeans=kmeans(scale(d_tsne_1), i)  
  
  # calculate clusters
  cluster_kmeans <- fit_cluster_kmeans$cluster
  
  # add cluster to pos DF
  dat <- cbind(d_tsne_1, cluster_kmeans)
  if(max(cluster_kmeans) >  1){
    s_score <- summary(silhouette(dat[,3], dist(dat[,1:2])))
    Sil_ind <- mean(s_score$clus.avg.widths)
  } else {
    Sil_ind <- -1
  }
  Score_df <- rbind(Score_df, Sil_ind)
}

optimized_clus <- which.max(Score_df)
fit_cluster_kmeans=kmeans(scale(d_tsne_1), optimized_clus)  
cluster_kmeans <- fit_cluster_kmeans$cluster




d_tsne_1_original$cl_kmeans = factor(fit_cluster_kmeans$cluster)

d_tsne_1_original$marker <- c('H&E', 'P16', 'P53', 'P63', 'uv-GATA3', 'Ki67', 'CK20', 'CK5-6', 'Her2Neu', 'Cyclin')

colnames(d_tsne_1_original) <- c('x', 'y', 'cluster', 'marker')





### create a copy
HE_part_stat <- Merge_stats
HE_part_stat[,2:11] <- as.data.frame(lapply(HE_part_stat[,2:11], as.numeric)) #<- sapply is here

colnames(HE_part_stat) <- c('TMA.core', 'H&E', 'P16', 'P53', 'P63', 'uv-GATA3', 'Ki67', 'CK20', 'CK5-6', 'Her2Neu', 'Cyclin')

## 1st: H&E sections
Sec_3 <- d_tsne_1_original[d_tsne_1_original$cluster == 3,4]

HE_part <- Merge_stats_melt[Merge_stats_melt$variable %in% Sec_3,]

# convert unit: um^2 to mm^2
TMA.core <- HE_part_stat$TMA.core
HE_part_stat <- cbind(TMA.core, HE_part_stat[Sec_3])

CV <- data.frame(matrix(nrow = 0, ncol = 1))
Mean <- data.frame(matrix(nrow = 0, ncol = 1))
QCoD <- data.frame(matrix(nrow = 0, ncol = 1))

for(i in 1:nrow(HE_part_stat)){ # num of stainings
  vec <- as.numeric((HE_part_stat[i, 2:4]))
  CV <- rbind(CV, sd(vec)/mean(vec))
  Mean <- rbind(Mean, mean(vec))
  QCoD <-  rbind(QCoD, (quantile(vec, 0.75) - quantile(vec, 0.25))/(quantile(vec, 0.75) + quantile(vec, 0.25)))
}
stats <- cbind(CV, Mean, QCoD)
colnames(stats) <- c('CV', 'Mean', 'QCoD')

HE_part_stat <- cbind(HE_part_stat, stats)



# create a copy
CK56_part_stat <- Merge_stats
CK56_part_stat[,2:11] <- as.data.frame(lapply(CK56_part_stat[,2:11], as.numeric)) #<- sapply is here

colnames(CK56_part_stat) <- c('TMA.core', 'H&E', 'P16', 'P53', 'P63', 'uv-GATA3', 'Ki67', 'CK20', 'CK5-6', 'Her2Neu', 'Cyclin')

## 2nd: CK56 sections
Sec_2 <- d_tsne_1_original[d_tsne_1_original$cluster == 2,4]

CK56_part <- Merge_stats_melt[Merge_stats_melt$variable %in% Sec_2,]

# convert unit: um^2 to mm^2
TMA.core <- CK56_part_stat$TMA.core
CK56_part_stat <- cbind(TMA.core, CK56_part_stat[Sec_2])

CV <- data.frame(matrix(nrow = 0, ncol = 1))
Mean <- data.frame(matrix(nrow = 0, ncol = 1))
QCoD <- data.frame(matrix(nrow = 0, ncol = 1))

for(i in 1:nrow(CK56_part_stat)){ # num of stainings
  vec <- as.numeric(CK56_part_stat[i, 2:6])
  CV <- rbind(CV, sd(vec)/mean(vec))
  Mean <- rbind(Mean, mean(vec))
  QCoD <-  rbind(QCoD, (quantile(vec, 0.75) - quantile(vec, 0.25))/(quantile(vec, 0.75) + quantile(vec, 0.25)))
}

stats <- cbind(CV, Mean, QCoD)
colnames(stats) <- c('CV', 'Mean', 'QCoD')

CK56_part_stat <- cbind(CK56_part_stat, stats)





##############
## 3rd: Ki67 #
##############
# create a copy
Ki67_part_stat <- Merge_stats
Ki67_part_stat[,2:11] <- as.data.frame(lapply(Ki67_part_stat[,2:11], as.numeric)) #<- sapply is here

colnames(Ki67_part_stat) <- c('TMA.core', 'H&E', 'P16', 'P53', 'P63', 'uv-GATA3', 'Ki67', 'CK20', 'CK5-6', 'Her2Neu', 'Cyclin')


Sec_1 <- d_tsne_1_original[d_tsne_1_original$cluster == 1,4]

Ki67_part <- Merge_stats_melt[Merge_stats_melt$variable %in% Sec_1,]

# convert unit: um^2 to mm^2
TMA.core <- Ki67_part_stat$TMA.core
Ki67_part_stat <- cbind(TMA.core, Ki67_part_stat[Sec_1])

Range <- data.frame(matrix(nrow = 0, ncol = 1))

for(i in 1:nrow(Ki67_part_stat)){ # num of stainings
  vec <- as.numeric(Ki67_part_stat[i, 2:3])
  Range <- rbind(Range, max(vec) - min(vec))
}


Ki67_part_stat <- cbind(Ki67_part_stat, Range)
colnames(Ki67_part_stat) <- c('TMA.core', 'Ki67', 'Cyclin', 'Range')



#------------- Find qualified positions ----------------#
HE_qualifiedPos <- HE_part_stat[HE_part_stat$QCoD < 0.1,]
HE_part_stat_list <- HE_qualifiedPos$TMA.core

CK56_qualifiedPos <- CK56_part_stat[CK56_part_stat$QCoD < 0.1,]
CK56_part_stat_list <- CK56_qualifiedPos$TMA.core

Ki67_qualifiedPos <- Ki67_part_stat[Ki67_part_stat$QCoD < 0.1,]
Ki67_part_stat_list <- Ki67_part_stat$TMA.core


Pos_list <- Reduce(intersect,list(HE_part_stat_list,CK56_part_stat_list))

write.csv(Pos_list, 'qualified_TMAcores.csv')


# -------------- Correlation analysis -------------------#


cor.test(HE_part_stat[,2], HE_part_stat[,6],  method = "spearman")

ggplot(HE_part_stat, aes(HE_part_stat[,2], HE_part_stat[,6])) +
  geom_point(shape = 21, colour = 'black', fill = 'white', size = 8) + 
  geom_smooth(method = "lm", se = FALSE)



#--------------------This section is for plotting--------------------#

#------- pca plot -------$

ggplot(d_pca_original, aes(x=x, y=y, color = as.factor(marker))) +  
  geom_point(shape = 21, size=4) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  xlab("") + ylab("") +
  ggtitle("PCA") +
  theme_bw(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  #stat_ellipse(aes(x, y, group = cluster, colour = factor(cluster)),size = 0.5) +
  geom_text(aes(x = x, y = y + 0.1, label = marker))
#scale_colour_brewer(palette = "Set2")

#------- t-SNE plot ------#
ggplot(d_tsne_1_original, aes(x=scale(x), y=scale(y), fill = as.factor(marker))) +  
  geom_point(size=15, shape = 21, stroke = 1) +
  guides(colour=guide_legend(override.aes=list(size=6), title = 'Marker')) +
  xlab("") + ylab("") +
  ggtitle("t-SNE") +
  theme_bw(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
#  stat_ellipse(aes(x, y, group = cluster, colour = factor(cluster)),size = 0.5) +
  geom_text(aes(x=scale(x), y=scale(y) - 0.1, label = marker, size = 10)) +
  labs(fill='Marker') 
#scale_colour_brewer(palette = "Set2")


#------- boxplot all ------#
write.csv(Merge_stats, 'Merge_stats.csv')

ratio <- max(Merge_stats_melt$value)/max(Merge_stats$QCoD)

jpeg( './coreSelection.jpeg',units="in", width=16, height=10, res=300)

ggplot() + 
  theme_bw(base_rect_size = 1) +
  geom_boxplot(data = Merge_stats_melt, aes(x = reorder(TMA.core, value, FUN = function(x){(quantile(x, 0.75) - quantile(x, 0.25))/(quantile(x, 0.75) + quantile(x, 0.25))}), y = value), outlier.shape = NA) +
  #geom_jitter(data = Merge_stats_melt, aes(x = reorder(TMA.core, value, FUN = function(x){sd(x)/mean(x)}), y = value),color = 'black', fill = 'white', shape=16, position=position_jitter(0.2)) +
  #geom_point(data = Merge_stats_melt,aes(x = reorder(TMA.core, value, FUN = function(x){sd(x)/mean(x)}), y = value),color = 'black', size = 1) +
  geom_dotplot(data = Merge_stats_melt,aes(x = reorder(TMA.core, value, FUN = function(x){(quantile(x, 0.75) - quantile(x, 0.25))/(quantile(x, 0.75) + quantile(x, 0.25))}), y = value),
               binaxis = 'y', stackdir = 'center', dotsize = 0.2, color = 'black', fill = 'white') +
  geom_line(data = Merge_stats, aes(x = reorder(TMA.core, QCoD), y = QCoD*ratio, group = 1), size =0.5, linetype = 1) +
  geom_point(data = Merge_stats, aes(x = reorder(TMA.core, QCoD), y = QCoD*ratio, group = 1), size = 1) +
  #geom_hline(yintercept = 0.1*ratio, linetype = 4) +
  #geom_vline(xintercept = 118, linetype = 4) +
  geom_point(data = Merge_stats, aes(x = reorder(TMA.core, QCoD), y = Mean, group = 1), color = 'red') +
  #geom_point(aes(x = 118, y = 0.1*ratio), color = 'red', size = 1) +
  scale_y_continuous(
    name = expression(paste('Core area (')~ mm^{2}~ ')'),
    sec.axis = sec_axis(~./ratio, name = 'Quartile Coeefficient of Dispersion')
  ) +
  xlab('Core labels') +
  theme(axis.text = element_text(angle = 90), axis.text.y = element_text(size = 10), axis.title = element_text(size = 15))
dev.off()


#------- boxplot HE plot ------#
ratio <- max(HE_part$value)/max(HE_part_stat$QCoD)

jpeg( './coreSelection.jpeg',units="in", width=16, height=10, res=300)

ggplot() + 
  theme_bw(base_rect_size = 1) +
  geom_boxplot(data = HE_part, aes(x = reorder(TMA.core, value, FUN = function(x){(quantile(x, 0.75) - quantile(x, 0.25))/(quantile(x, 0.75) + quantile(x, 0.25))}), y = value), outlier.shape = NA) +
  #geom_jitter(data = Merge_stats_melt, aes(x = reorder(TMA.core, value, FUN = function(x){sd(x)/mean(x)}), y = value),color = 'black', fill = 'white', shape=16, position=position_jitter(0.2)) +
  #geom_point(data = Merge_stats_melt,aes(x = reorder(TMA.core, value, FUN = function(x){sd(x)/mean(x)}), y = value),color = 'black', size = 1) +
  geom_dotplot(data = HE_part,aes(x = reorder(TMA.core, value, FUN = function(x){(quantile(x, 0.75) - quantile(x, 0.25))/(quantile(x, 0.75) + quantile(x, 0.25))}), y = value),
               binaxis = 'y', stackdir = 'center', dotsize = 0.2, color = 'black', fill = 'white') +
  geom_line(data = HE_part_stat, aes(x = reorder(TMA.core, QCoD), y = QCoD*ratio, group = 1), size =0.5, linetype = 1) +
  geom_point(data = HE_part_stat, aes(x = reorder(TMA.core, QCoD), y = QCoD*ratio, group = 1), size = 1) +
  geom_hline(yintercept = 0.1*ratio, linetype = 4) +
  #geom_vline(xintercept = 118, linetype = 4) +
  geom_point(data = HE_part_stat, aes(x = reorder(TMA.core, QCoD), y = Mean, group = 1), color = 'red') +
  #geom_point(aes(x = 118, y = 0.1*ratio), color = 'red', size = 1) +
  scale_y_continuous(
    name = expression(paste('Core area (')~ mm^{2}~ ')'),
    sec.axis = sec_axis(~./ratio, name = 'Quartile Coeefficient of Dispersion')
  ) +
  xlab('Core labels') +
  theme(axis.text = element_text(angle = 90), axis.text.y = element_text(size = 10), axis.title = element_text(size = 15))
dev.off()


#------- boxplot-CK5-6 plot ------#
ratio <- max(CK56_part$value) / max(CK56_part_stat$QCoD)

jpeg( './coreSelection.jpeg',units="in", width=16, height=10, res=300)

ggplot() + 
  theme_bw(base_rect_size = 1) +
  geom_boxplot(data = CK56_part, aes(x = reorder(TMA.core, value, FUN = function(x){(quantile(x, 0.75) - quantile(x, 0.25))/(quantile(x, 0.75) + quantile(x, 0.25))}), y = value), outlier.shape = NA) +
  #geom_jitter(data = Merge_stats_melt, aes(x = reorder(TMA.core, value, FUN = function(x){sd(x)/mean(x)}), y = value),color = 'black', fill = 'white', shape=16, position=position_jitter(0.2)) +
  #geom_point(data = Merge_stats_melt,aes(x = reorder(TMA.core, value, FUN = function(x){sd(x)/mean(x)}), y = value),color = 'black', size = 1) +
  geom_dotplot(data = CK56_part,aes(x = reorder(TMA.core, value, FUN = function(x){(quantile(x, 0.75) - quantile(x, 0.25))/(quantile(x, 0.75) + quantile(x, 0.25))}), y = value, fill = as.factor(CK56_part$variable)),
               binaxis = 'y', stackdir = 'center', dotsize = 0.1, color = 'black', fill = 'white') +
  geom_line(data = CK56_part_stat, aes(x = reorder(TMA.core, QCoD), y = QCoD*ratio, group = 1), size =0.5, linetype = 1) +
  geom_point(data = CK56_part_stat, aes(x = reorder(TMA.core, QCoD), y = QCoD*ratio, group = 1), size = 1) +
  geom_hline(yintercept = 0.1*ratio, linetype = 4) +
  #geom_vline(xintercept = 117, linetype = 4) +
  geom_point(data = CK56_part_stat, aes(x = reorder(TMA.core, QCoD), y = Mean, group = 1), color = 'red') +
  #geom_point(aes(x = 117, y = 0.1*ratio), color = 'red', size = 1) +
  scale_y_continuous(
    name = expression(paste('Core area (')~ mm^{2}~ ')'),
    sec.axis = sec_axis(~./ratio, name = 'Quartile Coeefficient of Dispersion')
  ) +
  xlab('Core labels') +
  theme(axis.text = element_text(angle = 90), axis.text.y = element_text(size = 10), axis.title = element_text(size = 15))
dev.off()

#------- boxplot- Ki67 plot ------#
source("http://www.sthda.com/upload/rquery_cormat.r")

require("corrplot")

corrplot(corr =cor(CK56_part_stat[, 2:6]),order="AOE",type="upper",tl.pos="tp")

corrplot(corr = cor(CK56_part_stat[, 2:6]),add=TRUE, type="lower", method="number",order="AOE", col="black",diag=FALSE,tl.pos="n", cl.pos="n")


Merge_stats[,2:11] <- as.data.frame(lapply(Merge_stats[,2:11], as.numeric)) #<- sapply is here

#cor.test(Merge_stats[,2], Merge_stats[,11],  method = "spearman")

ggplot(Ki67_part_stat, aes(Ki67_part_stat[,2], Ki67_part_stat[,3])) +
  geom_point(shape = 21, colour = 'black', fill = '#6CB0D6', size = 8) + 
  geom_smooth(method = "lm", se = FALSE, color = 'black') +
  xlab("Ki67") + ylab("Cyclin") +
  theme(axis.title = element_text(size = 30)) +
  #theme_bw(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) 


