## Hunter Colegrove
## FA dynamics
##
## Reads in position_data_all from compileVisFiles.R
## Plots number of tp53 mutations per replicate


require(foreach)
require(tidyverse)
require(glue)
library(gridExtra)
require(grid)
library(imager)
library(gtools)
library(png)
library(ggtext)


cluster = TRUE
#threshold_vector = c(10,50,100,400)
threshold_vector = c(108)
project_path = "2024-07-24-tp53_competition_sameMut_changeP"

## read in data from compileVisFiles.R
if(cluster){
  alldat <- read_csv(glue("/net/feder/vol1/project/HomeostaticEpithelium/dat/{project_path}/position_data_All.csv"))
}else{
  #alldat <- read_csv(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/position_data_All.csv"))
  alldat <- read_csv(glue("/Volumes/feder-vol1/project/HomeostaticEpithelium/dat/{project_path}/position_data_All.csv"))
}


######## 
## Average number of total tp53 mutations per replicate
## Includes tp53 + FA or tp53 + correction
## change to mut == 68 for only tp53 + FA
########


## filter for replicates that were not lost
## find n replicates for each tp53Block/year
completed_reps <- alldat %>%
  group_by(corrBlock, FAp53block, corrTime, replicate) %>%
  mutate(max_year = max(year)) %>%
  filter(max_year == 46) %>%
  group_by(corrBlock, year, FAp53block, corrTime) %>%
  summarise(completed_replicates = n_distinct(replicate)) %>%
  filter(!(year<1)) %>%
  filter(!(corrBlock == 0 & corrTime == 1)) %>%
  ungroup()
completed_reps %>% print(n=Inf)

## counts the number of cells that make up each tp53 clone within each replicate
## change to mut == 68 for only tp53 + FA
tp53clones <- alldat %>%
  filter(!(corrBlock == 0 & corrTime == 1)) %>% ## plots no correction case with corrBlock == 0
  filter(!(year < 1)) %>%
  dplyr::select(tp53Clone, replicate, corrBlock, year, mut, FAp53block) %>%
  filter(mut %in% c(68, 777)) %>% # tp53 + FA or tp53 + correction
  #filter(mut ==68) %>% # tp53 + FA only (no correction)
  group_by(tp53Clone, replicate, corrBlock, year, mut, FAp53block) %>%
  summarise(cells_count = n())
tp53clones

## filter data based on detection thresholds and combine all dataframes
df_list <- list()
for (val in threshold_vector) {
  duplicated_df <- tp53clones %>%
    filter(cells_count >= val) %>%
    mutate(threshold = val)
  df_list <- append(df_list, list(duplicated_df))
}
tp53_thresholds <- do.call(rbind, df_list)

## counts the unique number of tp53 mutations within each replicate at each time
unique_tp53Clone <- tp53_thresholds %>%
  group_by(corrBlock, year, replicate, threshold, FAp53block) %>%
  summarise(unique_count = n_distinct(tp53Clone)) %>%
  ungroup() %>%
  group_by(corrBlock, year, threshold, FAp53block) %>%
  summarise(total_unique_count = sum(unique_count)) %>%
  ungroup()
unique_tp53Clone %>% print(n=Inf)

## remove columns from completed dataframe
# completed_reps <- completed_reps %>%
#   dplyr::select(-mutRate, -corrTime)
# completed_reps %>% print(n=Inf)

## combine tp53 mut counts with the number of replicates completed in experiment
unique_tp53Clone_reps <- completed_reps %>%
  left_join(unique_tp53Clone, by = c("corrBlock", "year", "FAp53block"))
unique_tp53Clone_reps
## average mutations per replicate
ave_unique_count <- unique_tp53Clone_reps %>%
  mutate(ave_unique_count = total_unique_count/completed_replicates) %>%
  drop_na()
print(ave_unique_count)

## label mutation modifier values to read as an "X increase"
# custom_labels <- c("1" = "4X", "2" = "2X", "2.666" = "1.5X")
# blue_shades <- c("1" = "#03045e", "2" = "#0077b6", "2.666" = "#48cae4")
# 
# ## extra shades
# custom_labels <- c("0.375" = "10.6X", "0.5" = "8X", "1" = "4X", "2" = "2X", "2.666" = "1.5X")
# blue_shades <- c("0.375" = "#012A4A", "0.5" = "#01497c", "1" = "#2a6f97", "2" = "#468faf", "2.666" = "#89c2d9")

custom_labels <- c("0.01" = "0.01", "0.015" = "0.015", "0.02" = "0.02")
blue_shades <- c("0.02" = "#882255", "0.015" = "#CC6677", "0.01" = "#EE99AA")

## construct a label df of the final timepoint for line labels
dat_label <- ave_unique_count %>%
  group_by(corrBlock, FAp53block) %>%
  filter(year == max(year)) %>%
  ungroup()
dat_label

## facet labels
corrBlock_labels <- c("0" = "No correction", "0.01" = "0.01", "0.1" = "0.1")

## Plot using ggplot
plot <- ggplot(ave_unique_count, aes(x = year, y = ave_unique_count, color = factor(FAp53block))) +
  geom_line(size = 2) +
  labs(x = "Year", y = "TP53 mutation frequency", color = "FA p53 persistence") +
  theme_minimal() + 
  #theme(aspect.ratio = 0.75) + 
  scale_x_continuous(breaks = seq(1, 50, 5)) +
  facet_wrap(~corrBlock, scales = "fixed", nrow = 1, labeller = as_labeller(corrBlock_labels)) +
  theme(text = element_text(size = 8), panel.spacing = unit(2,'lines')) +
  scale_color_manual(values = blue_shades, labels = custom_labels) +
  #theme(legend.position = c(0.8, 0.8))
  theme(legend.position = "none") +
  geom_richtext(data = dat_label, aes(label = custom_labels[as.character(FAp53block)]), 
                color = "black", label.size = NA, fill = NA,
                label.margin = unit(4, "pt"),
                label.padding = unit(3, "pt"),
                hjust = 0, show.legend = FALSE) +
  coord_cartesian(clip = "off") +
  ylim(0,0.4) +
  theme(plot.margin = margin(5.5, 35, 5.5, 5.5))
  
show(plot)
## save plot to file
if(cluster){
  #ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/tp53_mutations_per_simulation_thresholds.png"), plot, width = 15, height = 6, units = "in", dpi = 300)
  ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-07-24-competition_sameMut_changeP/results/tp53_mutations_per_simulation_thresholds.png"), plot, width = 4.5, height = 2, units = "in", dpi = 300)
  
}else{
  ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/tp53_mutations_per_simulation_thresholds.png"), plot, width = 4.5, height = 2, units = "in", dpi = 300)
}

q()
############
##### Fill missing values
############

### Determine the missing replicates to be added.
summary_data <- alldat %>%
  filter(!(year < 1)) %>%
  filter(!(corrBlock == 0 & corrTime == 1)) %>%
  group_by(corrBlock, year, replicate, mutRate) %>%
  summarize(observations = n())
summary_data %>% print()
print(unique(summary_data$mutRate))


max_replicates = n_distinct(unique(summary_data$replicate[summary_data$corrBlock ==0]))

#summary_data$replicate[summary_data$corrBlock == 0]

missing_replicates <- lapply(years_total, function(year_value) {
  summary_data %>%
    filter(year == year_value) %>%
    group_by(corrBlock, mutRate) %>%
    print() %>%
    summarize(available_replicates = n_distinct(replicate)) %>%
    mutate(missing_replicates = max_replicates - available_replicates) %>%
    mutate(missing_replicates = ifelse(missing_replicates < 0, 0, missing_replicates)) %>%
    mutate(year = year_value)
})

missing_replicates <- do.call(rbind, missing_replicates)
print(missing_replicates, n=Inf)
missing_replicates <- missing_replicates %>% filter(missing_replicates >0)

tp53clones <- alldat %>%
  filter(!(corrBlock == 0 & corrTime == 1)) %>% ## plots no correction case with corrBlock == 0
  filter(!(year < 1)) %>%
  dplyr::select(tp53Clone, replicate, corrBlock, year, mut, mutRate) %>%
  filter(mut %in% c(68, 777)) %>% # tp53 + FA or tp53 + correction
  #filter(mut ==68) %>% # tp53 + FA only (no correction)
  group_by(tp53Clone, replicate, corrBlock, year, mut, mutRate) %>%
  summarise(cells_count = n())
tp53clones
tp53clones %>% print(n=Inf)
tp53clones

missing_replicates
for (i in 1:nrow(missing_replicates)){
  n_missing <- missing_replicates$missing_replicates[i]
  corrBlock_replace <- missing_replicates$corrBlock[i]
  year <- missing_replicates$year[i]
  mutRate <- missing_replicates$mutRate[i]
  for (j in 1:n_missing){
    n_reps = 100
    random_rep = sample(0:n_reps-1, size = 1)
    sampled_rep <- tp53clones %>%
      filter(corrBlock==0, year == year, replicate == random_rep, mutRate == mutRate)
    #print(sampled_rep)
    
    sampled_rep <- mutate(sampled_rep, corrBlock = corrBlock_replace)
    #print(corrBlock_replace)
    #print(sampled_rep)
    tp53clones <- bind_rows(tp53clones, sampled_rep)
    
  }
}
tp53clones

unique_reps <- tp53clones %>%
  group_by(corrBlock, year, mutRate) %>%
  summarize(unique_reps = n_distinct(replicate))
unique_reps


## filter for replicates that were not lost
## find n replicates for each tp53Block/year
completed_reps <- alldat %>%
  group_by(corrBlock, mutRate, corrTime, replicate) %>%
  mutate(max_year = max(year)) %>%
  filter(max_year == 46) %>%
  group_by(corrBlock, year, mutRate, corrTime) %>%
  summarise(completed_replicates = n_distinct(replicate)) %>%
  filter(!(year<1)) %>%
  filter(!(corrBlock == 0 & corrTime == 1)) %>%
  ungroup()
completed_reps %>% print(n=Inf)
### APPLY ONLY IF RANDOM SAMPLING
completed_reps <- mutate(completed_reps, completed_replicates = 100)


## counts the number of cells that make up each tp53 clone within each replicate
## change to mut == 68 for only tp53 + FA
# tp53clones <- alldat %>%
#   filter(!(corrBlock == 0 & corrTime == 1)) %>% ## plots no correction case with corrBlock == 0
#   filter(!(year < 1)) %>%
#   select(tp53Clone, replicate, corrBlock, year, mut) %>%
#   filter(mut %in% c(68, 777)) %>% # tp53 + FA or tp53 + correction
#   #filter(mut ==68) %>% # tp53 + FA only (no correction)
#   group_by(tp53Clone, replicate, corrBlock, year, mut) %>%
#   summarise(cells_count = n())
# tp53clones %>% print(n=500)


## filter data based on detection thresholds and combine all dataframes
df_list <- list()
for (val in threshold_vector) {
  duplicated_df <- tp53clones %>%
    filter(cells_count >= val) %>%
    mutate(threshold = val)
  df_list <- append(df_list, list(duplicated_df))
}
tp53_thresholds <- do.call(rbind, df_list)

## counts the unique number of tp53 mutations within each replicate at each time
unique_tp53Clone <- tp53_thresholds %>%
  group_by(corrBlock, year, replicate, threshold, mutRate) %>%
  summarise(unique_count = n_distinct(tp53Clone)) %>%
  ungroup() %>%
  group_by(corrBlock, year, threshold, mutRate) %>%
  summarise(total_unique_count = sum(unique_count)) %>%
  ungroup()
unique_tp53Clone %>% print(n=Inf)

## remove columns from completed dataframe
# completed_reps <- completed_reps %>%
#   dplyr::select(-mutRate, -corrTime)
# completed_reps %>% print(n=Inf)

## combine tp53 mut counts with the number of replicates completed in experiment
unique_tp53Clone_reps <- completed_reps %>%
  left_join(unique_tp53Clone, by = c("corrBlock", "year", "mutRate"))

## average mutations per replicate
ave_unique_count <- unique_tp53Clone_reps %>%
  mutate(ave_unique_count = total_unique_count/completed_replicates) %>%
  drop_na()
print(ave_unique_count)
ave_unique_count <- ave_unique_count %>%
  filter(threshold == 108)
ave_unique_count

#### Add 0's to data to show up on plot even if no mutations
#### must manually check data first if 0s can be added.
# start_years <- c(1, 6, 11)
# unique_corrBlock <- unique(ave_unique_count$corrBlock)
# unique_threshold <- unique(ave_unique_count$threshold)
# 
# # Create a dataframe with the additional rows
# additional_rows <- expand.grid(year = start_years, 
#                                ave_unique_count = 0, 
#                                corrBlock = unique_corrBlock, 
#                                threshold = unique_threshold)
# 
# # Combine with the original dataframe
# ave_unique_count_extended <- bind_rows(ave_unique_count, additional_rows) %>%
#   arrange(year, corrBlock, threshold)
# 
# print(ave_unique_count_extended, n=Inf)


## Plot using ggplot
num_colors <- length(unique(ave_unique_count$corrBlock))
my_palette <- rev(viridis(num_colors))
my_palette[1] <- "#BADE28FF"
plot <- ggplot(ave_unique_count, aes(x = year, y = ave_unique_count, color = factor(mutRate))) +
  geom_line(size = 3) +
  ylim(0, 0.45) +
  #labs(x = "Year", y = "TP53 mutation frequency", color = "FA mutant mutation rate") +
  labs(x = "Year", y = expression(paste(italic("TP53"), " mutation frequency")), color = "p") +
  theme_minimal() + 
  scale_color_manual(values = my_palette, 
                     breaks = unique(ave_unique_count$corrBlock), 
                     labels = unique(ave_unique_count$corrBlock)) +
  #theme(aspect.ratio = 0.75) + 
  scale_x_continuous(breaks = seq(1, 50, 5)) +
  #facet_wrap(~corrBlock, scales = "fixed", nrow = 1, labeller = as_labeller(corrBlock_labels)) +
  facet_wrap(~corrBlock, scales = "fixed", nrow = 1) +
  scale_color_manual(values = blue_shades, labels = custom_labels) +
  theme(legend.position = c(0.8, 0.8)) +
  theme(
    axis.line = element_line(color="black"), 
    axis.text = element_text(color="black", size=28),
    axis.title.x = element_text(color="black", size=28),
    axis.title.y = element_text(color="black", size=28),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(color="black", size=28), 
    legend.title = element_text(color="black", face="italic", size=28),
    #strip.background = element_blank(),
    #strip.text.x = element_blank(),
  )


show(plot)

## save plot to file
if(cluster){
  ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/tp53_mutations_per_simulation_thresholds.png"), plot, width = 20, height = 10, units = "in", dpi = 300)
}else{
  ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/tp53_mutations_per_simulation_thresholds_poster.png"), plot, width = 10, height = 6, units = "in", dpi = 300)
}


