## Hunter Colegrove
## FA dynamics
##
## Reads in position_data_all from compileVisFiles.R
## Outputs png frames of simulations for either animation or bulk viewing
## across replicates
## use compile_images.py for bulk viewing


require(foreach)
require(tidyverse)
require(glue)
library(gridExtra)
require(grid)
library(imager)
library(gtools)
library(png)


cluster = FALSE

project_path = "2024-01-12-SingleInjection_50years_DispersionSpacing"

## read in data from compileVisFiles.R
if(cluster){
  alldat <- read_csv(glue("/net/feder/vol1/project/HomeostaticEpithelium/dat/{project_path}/position_data_All.csv"))
}else{
  alldat <- read_csv(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/position_data_All.csv"))
}


#######################################
#########  Save frames for figure 3
#######################################

plotInitialCorrection <- function(rep_filter, sigma_filter, dose_filter, year_filter, corrBlock_filter){
  
  df_subset <- alldat %>% filter(y==0 & sigma == sigma_filter & dose == dose_filter & corrBlock == corrBlock_filter & replicate == rep_filter & year < year_filter)
  df_subset %>% filter(mut == 666) %>% print()
  plot <- ggplot(df_subset, aes(x = x, y = z, fill = factor(injSite))) +
    geom_tile(width = 1, height = 1) +
    scale_fill_manual(values = c("0" = "#44AA99", "-1" = "#DDDDDD", "1" = "#882255", "2" = "#88CCEE")) +
    labs(x = NULL, y = NULL, fill = "Corrected") +
    #labs(x = "X", y = "Y", fill = "Corrected") +
    geom_hline(yintercept = seq(-0.5,99.5), color = "black", size = 0.03) +
    geom_vline(xintercept = seq(-0.5,99.5), color = "black", size = 0.03) +
    theme(legend.position="none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = grid::unit(c(0,0,-1,-1), "mm"),
          plot.background=element_rect(fill="black", colour=NA)) +
    #labs(title = "Year: {as.integer(previous_state)}") +
    coord_cartesian(xlim = c(36, 58), ylim = c(36, 58)) # Setting plot limits for zoom in center
    #coord_cartesian(xlim = c(4, 95), ylim = c(4, 95)) # Setting plot limits for 100x100
  #coord_cartesian(xlim = c(2.6, 66.4), ylim = c(2.6, 66.4))  ## plot limits for 70x70
  show(plot)
  
  frame_out_dir <- "/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/"
  #frame <- glue("postCorrection_sigma_{sigma_filter}_dose_{dose_filter}.png")
  frame <- glue("postCorrection_sigma_{sigma_filter}_dose_{dose_filter}_thin.png")
  frame_out <- glue(frame_out_dir,frame)
  #ggsave(frame_out, plot=p, width=0.65, height=0.65, dpi=300)
  ggsave(glue({frame_out}), plot, width = 0.446, height = 0.446, units = "in", dpi = 300)
}

## dose 3 sigma 2
## dose 10 sigma 2
## dose 30 sigma 2

## dose 10 sigma 2
## dose 10 sigma 10
## dose 10 sigma 20

rep_filter <- 1
dose_filter <- 10
sigma_filter <- 20
year_filter <- 1 #(less than year 1 for injection time)
corrBlock_filter <- 0

alldat %>% filter(y==0 & sigma == sigma_filter & dose == dose_filter & corrBlock == corrBlock_filter & replicate == rep_filter & year < year_filter) %>% filter(x > 45 & x < 55 & z > 45 & z < 55) %>% print(n=Inf)
print(unique(alldat$replicate))
plotInitialCorrection(rep_filter, sigma_filter, dose_filter, year_filter, corrBlock_filter)

#############################################


# 
# if(cluster){
#   output_folder <- glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/animation_frames/")
# }else{
#   output_folder <- glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/animation_frames/")
# }
# 
# plotCorrectionEvents(tempdat, output_folder)


#########
#### Animation use ffmpeg on command line
#########

# ffmpeg -r 1 -i "/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-03-27-correction_tp53_competition_50years/results/block_0.01_rep_1/frame_%d_block_0.01_rep_1.png" -c:v libx264 -vf fps=1 -pix_fmt yuv420p "/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-03-27-correction_tp53_competition_50years/results/block_0.01_rep_1.mp4"
