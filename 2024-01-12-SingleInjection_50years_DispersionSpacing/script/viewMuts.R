
require(foreach)
require(tidyverse)
require(glue)
install.packages("ggplot2")
install.packages("gganimate")
require(gganimate)


## Project information and parameters used in experiment (enter from snakefile)
project_path = "2024-01-12-SingleInjection_50years_DispersionSpacing"
sigmas = c(2, 10, 20)
doses = c(3, 10, 30)
block_values = c(0, 0.001, 0.01, 0.1, 0.2, 0.5, 1)
print(p53_block_values)
corrTime = c(1)
nreplicates = 1
replicate = 0:nreplicates
#mutRates = c(4,40)
#replicate = c(4)

### adjust p53 block values to fit filename outputs
# scientific_threshold = 0.0001
# format_p53 <- function(x){
#   if (x >= scientific_threshold){
#     return(format(x, scientific = FALSE))
#   } else if(x == 0e+00){
#     return(format(x, scientific=FALSE))
#     } else {
#     return(format(x, scientific = TRUE))
#   }
# }
# p53_block_values <- sapply(p53_block_values, format_p53)

print(p53_block_values)


## Timesteps to be selected
timesteps = c(2, 1*364*4.5)
#timesteps = c(35217, 41769, 61425, 74529, 81081)
#timesteps = c(7371, 9009, 13923, 16380, 18018, 19656, 21294, 25389, 27027)


## Path to directory of raw data visualization files
path = glue("/Volumes/feder-vol1/project/HomeostaticEpithelium/dat/{project_path}/results/")
print(path)

## combines all data into a large tibble
alldat <- foreach(block = block_values, .combine = "rbind") %:%
            foreach(rep = replicate, .combine = "rbind") %:%
              foreach(timestep = timesteps, .combine = "rbind") %:%
                foreach(sigma = sigmas, .combine = "rbind") %:%
                  foreach(dose = doses, .combine = "rbind") %do%  {
                    
  f <- glue("VisFile_block_{block}_growth_1.0_sigma_{sigma}_dose_{dose}_correctionTime_{corrTime[1]}_{rep}.{timestep}.txt")
  cat("Opening file:", f, "\n")
  
  #paste("Opening file: fakedat/", f)
  full_path <- paste0(path, f)
  print(full_path)
  column_names = c('x','z','y', 'h','s','v','alpha','mut','injSite','tp53Clone')
  if (file.exists(full_path)){
    data <- tibble(read.csv(full_path, sep = '\t', header = FALSE, col.names = column_names)) %>% 
      filter(y == 0) %>%
      mutate(replicate = rep) %>%
      mutate(corrBlock = block) %>%
      #mutate(tp53Block = p53_val) %>%
      mutate(year = timestep/4.5/364) %>%
      mutate(sigma = sigma) %>%
      mutate(dose = dose)
      #mutate(mutRate = mutRate)
      #filter(tp53Clone != -1)
  }
}

print(alldat)

## filter for tp53 clones greater than a detection threshold
# result <- alldat %>%
#   group_by(tp53Block, replicate, tp53Clone, year, mutRate) %>%
#   filter(tp53Clone != -1) %>%
#   summarize(count = n()) %>%
#   print(n=Inf) %>%
#   filter(count >= 400) %>%
#   print() %>%
#   group_by(tp53Block, replicate, year, mutRate) %>%
#   summarize(unique_tp53Clones = n()) %>%
#   print() %>%
#   group_by(tp53Block, year, mutRate) %>%
#   summarize(total_unique_tp53Clones = sum(unique_tp53Clones)) %>%
#   mutate(tp53clone_frequency = total_unique_tp53Clones / (nreplicates+1)) %>%
#   print(n = Inf)
# result <- result %>%
#   select(-total_unique_tp53Clones)
# result <- result %>%
#   rename(blockProb = tp53Block) %>%
#   rename(age = year) %>%
#   rename(AverageMutationCount = tp53clone_frequency)
# write.csv(result, glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/countClones_tp53_R.csv"), row.names = FALSE)


## VAF
# resultVAF <- alldat %>%
#   group_by(tp53Block, replicate, tp53Clone, year, mutRate) %>%
#   filter(tp53Clone != -1) %>%
#   summarize(count = n()) %>%
#   print() %>%
#   filter(count >= 400) %>%
#   print() %>%
#   #mutate(VAF = count/(141*141/2)) %>%
#   mutate(VAF = count/(4900)) %>%
#   print(n=Inf)
# 
# resultVAF <- resultVAF %>%
#   rename(blockProb = tp53Block) %>%
#   rename(age = year)
# write.csv(resultVAF, glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/VAF_tp53_R.csv"), row.names = FALSE)




#### Grid of tissue sections
# alldat %>% ggplot(aes(x = x, y = z, fill = as.factor(tp53Clone))) +
#   geom_tile(color = "white") +
#   facet_grid(tp53Block ~ replicate) +
#   labs(x = "X", y = "Y") +
#   theme_minimal()





##########
#### Single tissue TestOuput
##########

#path = glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/TestOutput/VisFiles/")
# f <- glue("VisFile_block_{block_values[1]}_growth_1.0_sigma_{sigma[1]}_dose_{dose[1]}_p53_1e-06_correctionTime_{corrTime[1]}_3.{timestep[1]}.txt")
# 
# cat("Opening file:", f, "\n")
# full_path <- paste0(path, f)
# print(full_path)
# 
# column_names = c('x','z','y', 'h','s','v','alpha','mut','injSite','tp53Clone')
# if (file.exists(full_path)){
#   data <- tibble(read.csv(full_path, sep = '\t',header = FALSE, col.names = column_names)) 
#   print(data)
# }
# data <- data %>% filter(y==0)
# print(data)

##########
#### Visualize grid tissue basal layer using alldat tibble
##########

alldat 
write.csv(alldat, glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/position_data_All.csv"), row.names = FALSE)


alldat <- read_csv(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/position_data_All.csv"))
alldat
print(unique(alldat$year))
alldat

tempdat <- alldat %>%
  filter(replicate %in% c(0)) %>% ## filter for vis replicates
  filter(year < 1) %>%
  #filter(year == 2/365/4.5) %>%
  filter(corrBlock == 1) %>%
  filter(sigma == 2) %>%
  filter(dose == 3)
  #filter(tp53Block == 0.005) %>%
  #filter(mutRate == 4) %>%
  #filter(tp53Clone > -1) %>%
  #mutate(replicate = as.factor(replicate))
print(tempdat)

# Calculate the count of each unique value
# value_counts <- tempdat %>%
#   filter(tp53Clone != -1) %>%
#   group_by(tp53Clone) %>%
#   summarize(count = n())
# print(value_counts)

tempdat %>% filter(injSite == 0)

plotCorrectionEvents <- function(df){
  
  df <- df %>%
    mutate(injSite = ifelse((x == 49 & z == 49) |
                              (x == 51 & z == 50), -1, injSite))
  df <- df %>%
    mutate(injSite = ifelse((x == 50 & z == 46) |
                              (x == 49 & z == 45), 0, injSite))

  # df <- df %>% 
  #   mutate(injSite = -1)
  
  p <- ggplot(df, aes(x = x, y = z, fill = factor(injSite))) +
    geom_tile(width = 1, height = 1) +
    scale_fill_manual(values = c("0" = "#44AA99", "-1" = "#DDDDDD")) +
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
    coord_cartesian(xlim = c(36, 58), ylim = c(36, 58)) # Setting plot limits
    #coord_cartesian(xlim = c(4, 95), ylim = c(4, 95)) # Setting plot limits
    #geom_segment(data = segment_lines, aes(x = x, y = y, xend = xend, yend = yend))
  
  # anim <- p + transition_time(year) +
  #   ease_aes('linear') +
  #   shadow_mark(past = TRUE, future = FALSE)
  
  return(p)
}

plotCorrectionEvents(tempdat)
plot <- plotCorrectionEvents(tempdat)
animate(plot)

ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/preCorrection_small.png"), plot, width = 6, height = 6, units = "in", dpi = 300)
ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/postCorrection_big.png"), plot, width = 6, height = 6, units = "in", dpi = 300)

ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/postCorrection_sigma_2_dose_3.png"), plot, width = 6, height = 6, units = "in", dpi = 300)
ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/postCorrection_sigma_20_dose_30.png"), plot, width = 6, height = 6, units = "in", dpi = 300)
ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/postCorrection_sigma_20_dose_3.png"), plot, width = 6, height = 6, units = "in", dpi = 300)
ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/postCorrection_sigma_20_dose_10.png"), plot, width = 6, height = 6, units = "in", dpi = 300)
ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/postCorrection_sigma_10_dose_3.png"), plot, width = 6, height = 6, units = "in", dpi = 300)
ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/postCorrection_sigma_10_dose_30.png"), plot, width = 6, height = 6, units = "in", dpi = 300)
ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/postCorrection_sigma_10_dose_10.png"), plot, width = 6, height = 6, units = "in", dpi = 300)
ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/postCorrection_sigma_2_dose_30.png"), plot, width = 6, height = 6, units = "in", dpi = 300)
ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/postCorrection_sigma_2_dose_10.png"), plot, width = 6, height = 6, units = "in", dpi = 300)


# Create a color palette based on the count of unique values
# if(length(value_counts > 0)){
#   color_palette <- scales::hue_pal()(length(value_counts$count))
# }else{color_palette = character(0)}
# print(color_palette)
# # Map colors to each unique value based on its count
# color_key <- setNames(color_palette, value_counts$tp53Clone)
# color_key[["-1"]] <- "grey"
# print(color_key)
# 
# ggplot(tempdat, aes(x = x, y = z, fill = as.factor(tp53Clone))) +
#   geom_tile(color = "black", width = 1, height = 1, show.legend=FALSE) +
#   labs(x = "X", y = "Y") +
#   scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
#   scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
#   scale_fill_manual(values = color_key) +
#   theme_void() + 
#   #facet_wrap(~ replicate, scales = "free", ncol = 2) + 
#   theme(aspect.ratio = 1)



data %>% filter(y == 0) %>%
  filter(mut == 68) %>%
  print(n = Inf)



## plot viable cells
viableOnly <- data %>%
  mutate(mut = ifelse(mut>= 0,0,mut)) %>%
  print(n = Inf)

ggplot(viableOnly, aes(x = x, y = z, fill = as.factor(mut))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("black", "grey"), labels = c("Dead", "Alive")) +
  labs(x = "X", y = "Y") +
  theme_minimal()
  