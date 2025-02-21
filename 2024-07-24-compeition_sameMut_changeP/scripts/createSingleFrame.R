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


cluster = TRUE

project_path = "2024-07-24-competition_sameMut_changeP"

## read in data from compileVisFiles.R
if(cluster){
  alldat <- read_csv(glue("/net/feder/vol1/project/HomeostaticEpithelium/dat/2024-07-24-tp53_competition_sameMut_changeP/position_data_All.csv"))
}else{
  #alldat <- read_csv(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/position_data_All.csv"))
  alldat <- read_csv(glue("/Volumes/feder-vol1/project/HomeostaticEpithelium/dat/2024-07-24-tp53_competition_sameMut_changeP/position_data_All.csv"))
}

print(alldat)

#########################################
#########  Pull a frame for visualization
#########################################

## adjust injSite values to easily identify corrected and tp53
## start visualizing at year 1 exclude initial time point
tempdat <- alldat
tempdat$injSite[tempdat$tp53Clone != -1] <- 1
tempdat$injSite[tempdat$mut == 777] <- 2
tempdat <- tempdat %>%
  filter(!(year < 1)) %>%
  filter(!(corrBlock == 0 & corrTime == 1))
#b
# replicates_noCorr <- c(84, 91, 51, 36, 77, 78)
# replicates_0.01 <- c(48, 81, 59, 50, 92, 96)
# replicates_0.1 <- c(51, 72, 66, 62,  3, 50)
# replicates_0.01_noCorrSampled <- c(45, 31, 13, 35)
# corrBlocks_list <- c(0,0.01,0.1,0.011)

# set 1
# replicates_noCorr <- c(97, 13, 60, 14,  9, 19)
# replicates_0.01 <- c(14, 51, 50, 47, 60, 18)
# replicates_0.1 <- c(44, 46, 28, 56, 76, 62)
# replicates_0.01_noCorrSampled <- c(91, 60, 54,  4, 44)
# corrBlocks_list <- c(0,0.01,0.1,0.011)

# set 2
replicates_noCorr <- c(65, 46, 59, 74, 17, 84)
replicates_0.01 <- c(63, 40, 60, 88, 11, 81)
replicates_0.1 <- c(55, 36, 65, 32, 80, 20)
replicates_0.01_noCorrSampled <- c(45, 59, 35, 16, 27, 15)
replicates_0.1_noCorrSampled <- c(81)
corrBlocks_list <- c(0,0.01,0.1,0.011, 0.11)

plotCorrectionEvents <- function(df, output_folder){

  for (k in corrBlocks_list){
    if(k == 0){
      k_original = 0
      reps_list <- replicates_noCorr
    }
    if(k ==0.01){
      k_original=0.01
      reps_list <- replicates_0.01
    }
    if(k==0.1){
      k_original = 0.1
      reps_list <- replicates_0.1
    }
    if(k==0.011){ ## if 0.01 samples do not have timepoints, resample from no correction
      reps_list <- replicates_0.01_noCorrSampled
      k_original = 0.01
    }
    if(k==0.11){ ## if 0.01 samples do not have timepoints, resample from no correction
      reps_list <- replicates_0.1_noCorrSampled
      k_original = 0.1
    }
    for (j in reps_list){
      if(k==0.011 | k==0.11){
        k = 0
      }
      df_subset <- subset(df, year==46 & replicate==j & corrBlock==k & FAp53block == 0.04)
      ## dont make frame if there is no timepoint
      if ( dim(df_subset)[1] == 0){
        print(paste("Timepoint doesn't exist: corrBlock", k, " rep", j))
        next
      }
      
      p <- ggplot(df_subset, aes(x = x, y = z, fill = factor(injSite))) +
        geom_tile(width = 1, height = 1) +
        scale_fill_manual(values = c("0" = "#44AA99", "-1" = "#DDDDDD", "1" = "#882255", "2" = "#88CCEE")) +
        labs(x = NULL, y = NULL, fill = "Corrected") +
        #labs(x = "X", y = "Y", fill = "Corrected") +
        #geom_hline(yintercept = seq(-0.5,99.5), color = "black", size = 0.025) +
        #geom_vline(xintercept = seq(-0.5,99.5), color = "black", size = 0.025) +
        theme(legend.position="none",
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              plot.margin = grid::unit(c(0,0,-1,-1), "mm"),
              plot.background=element_rect(fill="black", colour=NA)) +
        #labs(title = "Year: {as.integer(previous_state)}") +
        #coord_cartesian(xlim = c(36, 58), ylim = c(36, 58)) # Setting plot limits for zoom in center
        #coord_cartesian(xlim = c(4, 95), ylim = c(4, 95)) # Setting plot limits for 100x100
        coord_cartesian(xlim = c(2.6, 66.4), ylim = c(2.6, 66.4))  ## plot limits for 70x70
      #frame_out <- paste(output_folder,"figureFrame_", "_block_", k_original, "_rep_", j, ".png", sep="")
      frame_out <- paste(output_folder,"figureFrame_set2_", "_block_", k_original, "_rep_", j, ".png", sep="")
      ggsave(frame_out, plot=p, width=0.25, height=0.25, dpi=300)
      #show(frame_out)
    }
  }
}

if(cluster){
  output_folder <- glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/figure5_frames/")
}else{
  output_folder <- glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/figure5_frames/")
}

plotCorrectionEvents(tempdat, output_folder)


#########
#### Animation use ffmpeg on command line
#########

# ffmpeg -r 1 -i "/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-03-27-correction_tp53_competition_50years/results/block_0.01_rep_1/frame_%d_block_0.01_rep_1.png" -c:v libx264 -vf fps=1 -pix_fmt yuv420p "/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-03-27-correction_tp53_competition_50years/results/block_0.01_rep_1.mp4"