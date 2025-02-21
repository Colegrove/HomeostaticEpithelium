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

project_path = "2024-06-13-tp53_competition_changeMut_sameP"

## read in data from compileVisFiles.R
if(cluster){
  alldat <- read_csv(glue("/net/feder/vol1/project/HomeostaticEpithelium/dat/{project_path}/position_data_All.csv"))
}else{
  alldat <- read_csv(glue("/Volumes/feder-vol1/project/HomeostaticEpithelium/dat/{project_path}/position_data_All.csv"))
}


#######################################
#########  Save frames for animation
#######################################

## adjust injSite values to easily identify corrected and tp53
## start visualizing at year 1 exclude initial time point
tempdat <- alldat
tempdat$injSite[tempdat$tp53Clone != -1] <- 1
tempdat$injSite[tempdat$mut == 777] <- 2
tempdat <- tempdat %>%
  filter(!(year < 1)) %>%
  filter(!(corrBlock == 0 & corrTime == 1))

plotCorrectionEvents <- function(df, output_folder){
  print(unique(df$mutRate))
  if (!dir.exists(output_folder)){
    dir.create(output_folder)
  }
  for (k in unique(df$corrBlock)){
    for (m in unique(df$mutRate)){
      for (j in unique(df$replicate)){
        cat(glue("Saving corrected blocking probability: {k}, mutrate: {m}, replicate: {j}"),"\n")
        count = 0
        for (i in unique(df$year)){
          count = count + 1
          df_subset <- subset(df, year==i & replicate==j & corrBlock==k & mutRate==m)
          
          ## dont make frame if there is no timepoint
          if ( dim(df_subset)[1] == 0){
            next
          }
          
          p <- ggplot(df_subset, aes(x = x, y = z, fill = factor(injSite))) +
            geom_tile(width = 1, height = 1) +
            scale_fill_manual(values = c("0" = "#44AA99", "-1" = "#DDDDDD", "1" = "#882255", "2" = "#88CCEE")) +
            labs(x = NULL, y = NULL, fill = "Corrected") +
            #labs(x = "X", y = "Y", fill = "Corrected") +
            geom_hline(yintercept = seq(-0.5,99.5), color = "black", size = 0.5) +
            geom_vline(xintercept = seq(-0.5,99.5), color = "black", size = 0.5) +
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
          #geom_segment(data = segment_lines, aes(x = x, y = y, xend = xend, yend = yend))
          frame_folder <- paste(paste(output_folder, "block_", k, "_mutRate_", m, "_rep_", j, "/", sep = ""))
          #cat(paste(frame_folder, "\n"))
          frame_out <- paste(frame_folder, paste("frame_", count, "_block_", k, "_mutRate_", m, "_rep_", j, ".png", sep=""), sep="")
          #cat(paste(frame_out, "\n"))
          if (!dir.exists(frame_folder)){
            dir.create(frame_folder)
          }
          ggsave(frame_out, plot=p, width=6, height=6, dpi=300)
        }
      }
    }
  }
}

if(cluster){
  output_folder <- glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/animation_frames/")
}else{
  output_folder <- glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/animation_frames/")
}

plotCorrectionEvents(tempdat, output_folder)


#########
#### Animation use ffmpeg on command line
#########

# ffmpeg -r 1 -i "/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-03-27-correction_tp53_competition_50years/results/block_0.01_rep_1/frame_%d_block_0.01_rep_1.png" -c:v libx264 -vf fps=1 -pix_fmt yuv420p "/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-03-27-correction_tp53_competition_50years/results/block_0.01_rep_1.mp4"
