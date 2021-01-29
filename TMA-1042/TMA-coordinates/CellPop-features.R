# This file is to classify cells
library(Rtsne)
library(ggplot2)
library(cluster)

cellFeatures1 <- read.csv('./cellFeatures/K-2/cellFeatures.csv', row.names = 1)
cellFeatures2 <- read.csv('./cellFeatures/C-4/cellFeatures.csv', row.names = 1)
cellFeatures3 <- read.csv('./cellFeatures/L-2/cellFeatures.csv', row.names = 1)
cellFeatures4 <- read.csv('./cellFeatures/H-5/cellFeatures.csv', row.names = 1)
cellFeatures5 <- read.csv('./cellFeatures/I-1/cellFeatures.csv', row.names = 1)
cellFeatures6 <- read.csv('./cellFeatures/E-9/cellFeatures.csv', row.names = 1)
cellFeatures7 <- read.csv('./cellFeatures/A-7/cellFeatures.csv', row.names = 1)

cellFeatures <- rbind(cellFeatures1, cellFeatures2, cellFeatures3, cellFeatures4, cellFeatures5, cellFeatures6, cellFeatures8)

tsne_df <- as.data.frame(lapply(cellFeatures[,c(4,5)], as.numeric)) #<- sapply is here

tsne_model_1 = Rtsne(tsne_df, check_duplicates=FALSE)

d_tsne_1 = as.data.frame(tsne_model_1$Y)  

## keeping original data
d_tsne_1_original <- d_tsne_1
colnames(d_tsne_1_original) <- c('x', 'y')
# score vector

Score_df <- matrix(nrow = 0, ncol = 1)

for(i in 1:5){
  
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


pca_plot <- prcomp(cellFeatures[,c(4,5, 6, 8, 9)], scale = TRUE)

plot(pca_plot$x)
#------- t-SNE plot ------#
ggplot(d_tsne_1_original, aes(x=scale(x), y=scale(y))) +  
  geom_point() +
  guides(colour=guide_legend(override.aes=list(size=6), title = 'Marker')) +
  xlab("") + ylab("") +
  ggtitle("t-SNE") +
  theme_bw(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  #  stat_ellipse(aes(x, y, group = cluster, colour = factor(cluster)),size = 0.5) +
#  geom_text(aes(x=scale(x), y=scale(y) - 0.1, label = marker, size = 10)) +
  labs(fill='Marker') 
#scale_colour_brewer(palette = "Set2")