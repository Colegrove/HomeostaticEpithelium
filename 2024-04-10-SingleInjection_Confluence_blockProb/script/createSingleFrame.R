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

project_path = "2024-04-10-SingleInjection_Confluence_blockProb"

## read in data from compileVisFiles.R
if(cluster){
  alldat <- read_csv(glue("/net/feder/vol1/project/HomeostaticEpithelium/dat/{project_path}/position_data_All.csv"))
}else{
  #alldat <- read_csv(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/position_data_All.csv"))
  alldat <- read_csv(glue("/Volumes/feder-vol1/project/HomeostaticEpithelium/dat/{project_path}/position_data_All.csv"))
}

# print(unique(alldat$mutRate))
# alldat <- alldat %>%
#   filter(mutRate == 0.5) %>%
#   filter(year == 46) %>%
#   group_by(corrBlock) %>%
#   summarise(unique_replicates = list(unique(replicate)))
# 
# print(alldat$unique_replicates[1])
# print(alldat$unique_replicates[2])
# print(alldat$unique_replicates[3])
# q()
# mutRate_unique <- unique(alldat$mutRate)


#########################################
#########  Pull a frame for visualization
#########################################

## adjust injSite values to easily identify corrected and tp53
## start visualizing at year 1 exclude initial time point
tempdat <- alldat
#tempdat$injSite[tempdat$tp53Clone != -1] <- 1
#tempdat$injSite[tempdat$mut == 777] <- 2
# tempdat <- tempdat %>%
#   filter((year < 1))

unique(tempdat$corrBlock)
unique(tempdat$year)
unique(tempdat$replicate)
corr1 <- tempdat %>% filter(corrBlock == 1)
unique(corr1$year)
plotCorrectionEvents <- function(df, output_folder, t, rep, corr, confluence){
      df_subset <- subset(df, year==t & replicate==rep & corrBlock==corr)
      # dont make frame if there is no timepoint
      if ( dim(df_subset)[1] == 0){
        print(paste("Timepoint doesn't exist: corrBlock", corr, " rep", rep, " year", t))

      }

      p <- ggplot(df_subset, aes(x = x, y = z, fill = factor(injSite))) +
        geom_tile(width = 1, height = 1) +
        scale_fill_manual(values = c("0" = "#44AA99", "-1" = "#DDDDDD", "1" = "#882255", "2" = "#88CCEE")) +
        labs(x = NULL, y = NULL, fill = "Corrected") +
        #labs(x = "X", y = "Y", fill = "Corrected") +
        geom_hline(yintercept = seq(-0.5,99.5), color = "black", size = 0.025) +
        geom_vline(xintercept = seq(-0.5,99.5), color = "black", size = 0.025) +
        theme(legend.position="none",
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              plot.margin = grid::unit(c(0,0,-1,-1), "mm"),
              plot.background=element_rect(fill="black", colour=NA)) +
        #labs(title = "Year: {as.integer(previous_state)}") +
        #coord_cartesian(xlim = c(36, 58), ylim = c(36, 58)) # Setting plot limits for zoom in center
        coord_cartesian(xlim = c(4, 95), ylim = c(4, 95)) # Setting plot limits for 100x100
        #coord_cartesian(xlim = c(2.6, 66.4), ylim = c(2.6, 66.4))  ## plot limits for 70x70
      frame_out <- paste(output_folder,"figureFrame_", "_block_", corr, "_rep_", rep, "_year_", t, "_", confluence , ".png", sep="")
      ggsave(frame_out, plot=p, width=0.75, height=0.75, dpi=300)
  
}

if(cluster){
  output_folder <- glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/basal_frames/")
}else{
  output_folder <- glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/basal_frames/")
}


rep = 2
corrBlock_val = 0.2

all_timepoints <- tempdat %>% filter(replicate == rep & corrBlock == corrBlock_val)
all_timepoints
unique(all_timepoints$year)
max_year <- max(unique(all_timepoints$year))
min_year <- min(unique(all_timepoints$year))
last_timepoint <- all_timepoints %>% filter(year == max_year)
first_timepoint <- all_timepoints %>% filter(year == min_year)

nCorr_last <- last_timepoint %>% filter(injSite >= 0) %>% summarise(n())
if(nCorr_last >= 8000){
  confluence = "CONFLUENCE"
}else{confluence = "LOSS"}
confluence_first = "START"

plotCorrectionEvents(last_timepoint, output_folder, max_year, rep, corrBlock, confluence)

plotCorrectionEvents(first_timepoint, output_folder, min_year, rep, corrBlock, confluence_first)

#########
#### Animation use ffmpeg on command line
#########

# ffmpeg -r 1 -i "/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-03-27-correction_tp53_competition_50years/results/block_0.01_rep_1/frame_%d_block_0.01_rep_1.png" -c:v libx264 -vf fps=1 -pix_fmt yuv420p "/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-03-27-correction_tp53_competition_50years/results/block_0.01_rep_1.mp4"