###########################################################
# This script is used to rescale the points for each core #
###########################################################

library(jpeg)
library(grid)

setwd("~/Desktop/TMA - coregistration")

Ki67_img <- readJPEG('Core_images/D-12/Ki67.jpg')






for_plot <-  ggplot() +
  annotation_custom(rasterGrob(Ki67_img, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc"),interpolate = FALSE))
plot(for_plot)
                    .#0,  31.872, -18.024, 0) 
  #xlab('x,mm')+
  #ylab('y,mm')