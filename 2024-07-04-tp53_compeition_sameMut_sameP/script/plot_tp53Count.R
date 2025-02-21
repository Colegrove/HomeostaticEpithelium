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
require(viridis)


cluster = TRUE
#threshold_vector = c(10,50,100,400)
threshold_vector = c(108)
project_path = "2024-07-04-tp53_competition_sameMut_sameP"

## read in data from compileVisFiles.R
if(cluster){
  alldat <- read_csv(glue("/net/feder/vol1/project/HomeostaticEpithelium/dat/{project_path}/position_data_All.csv"))
}else{
  #alldat <- read_csv(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/position_data_All.csv"))
  alldat <- read_csv(glue("/Volumes/feder-vol1/project/HomeostaticEpithelium/dat/{project_path}/position_data_All.csv"))
}
alldat
#alldat1 <- read_csv(glue("/Volumes/feder-vol1/project/HomeostaticEpithelium/dat/2024-06-13-tp53_competition_sameMut_sameP/position_data_All.csv"))
#alldat2 <- read_csv(glue("/Volumes/feder-vol1/project/HomeostaticEpithelium/dat/2024-07-03-tp53_competition_sameMut_sameP/position_data_All.csv"))
alldat1 <- read_csv(glue("/net/feder/vol1/project/HomeostaticEpithelium/dat/2024-06-13-tp53_competition_sameMut_sameP/position_data_All.csv"))
alldat1
alldat2 <- read_csv(glue("/net/feder/vol1/project/HomeostaticEpithelium/dat/2024-07-03-tp53_competition_sameMut_sameP/position_data_All.csv"))
alldat2
alldat1 <- alldat1 %>% mutate(replicate = replicate + 100)
alldat2 <- alldat2 %>% mutate(replicate = replicate + 200)
alldat <- bind_rows(alldat, alldat1, alldat2)


######## 
## Average number of total tp53 mutations per replicate
## Incluces tp53 + FA or tp53 + correction
## change to mut == 68 for only tp53 + FA
########

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
## counts the number of cells that make up each tp53 clone within each replicate
## change to mut == 68 for only tp53 + FA
tp53clones <- alldat %>%
  filter(!(corrBlock == 0 & corrTime == 1)) %>% ## plots no correction case with corrBlock == 0
  filter(!(year < 1)) %>%
  dplyr::select(tp53Clone, replicate, corrBlock, year, mut) %>%
  filter(mut %in% c(68, 777)) %>% # tp53 + FA or tp53 + correction
  #filter(mut ==68) %>% # tp53 + FA only (no correction)
  group_by(tp53Clone, replicate, corrBlock, year, mut) %>%
  summarise(cells_count = n())

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
  group_by(corrBlock, year, replicate, threshold) %>%
  summarise(unique_count = n_distinct(tp53Clone)) %>%
  ungroup() %>%
  group_by(corrBlock, year, threshold) %>%
  summarise(total_unique_count = sum(unique_count)) %>%
  ungroup()
unique_tp53Clone %>% print(n=Inf)

## remove columns from completed dataframe
completed_reps <- completed_reps %>%
  dplyr::select(-mutRate, -corrTime)
completed_reps %>% print(n=Inf)

## combine tp53 mut counts with the number of replicates completed in experiment
unique_tp53Clone_reps <- completed_reps %>%
  left_join(unique_tp53Clone, by = c("corrBlock", "year"))
unique_tp53Clone_reps
## average mutations per replicate
ave_unique_count <- unique_tp53Clone_reps %>%
  mutate(ave_unique_count = total_unique_count/completed_replicates) %>%
  drop_na()
print(ave_unique_count)



## label mutation modifier values to read as an "X increase"
custom_labels <- c("0" = "0", "0.01" = "0.01", "0.1" = "0.1")
green_shades <- c("0.1" = "#44aa99", "0.01" = "#67c2b3", "0" = "#bbe4dd")


plot <- ggplot(ave_unique_count, aes(x = year, y = ave_unique_count, color = factor(corrBlock))) +
  geom_line(size = 1.5) +
  labs(x = "Year", y = "TP53 mutation frequency", color = "FA mutant mutation rate") +
  theme_minimal(base_size = 8) + 
  #theme(aspect.ratio = 0.75) + 
  scale_x_continuous(breaks = seq(1, 50, 5)) +
  theme(text = element_text(size = 8), panel.spacing = unit(2,'lines')) +
  scale_color_manual(values = green_shades, labels = custom_labels) +
  theme(legend.position = c(0.8, 0.8)) +
  theme(plot.margin = margin(5.5, 35, 5.5, 5.5))

## Plot using ggplot
# plot <- ggplot(ave_unique_count, aes(x = year, y = ave_unique_count, color = factor(corrBlock))) +
#   geom_line() +
#   labs(x = "Year", y = "Average number of tp53 mutations per replicate", color = "Blocking Probability b") +
#   theme_minimal() + 
#   #theme(aspect.ratio = 0.75) + 
#   scale_x_continuous(breaks = seq(1, 50, 5)) +
#   facet_wrap(~threshold, scales = "fixed", nrow = 1) +
#   theme(text = element_text(size = 8))
#   #theme(strip.text = element_text(size = 20))




show(plot)
## save plot to file
if(cluster){
  ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/tp53_mutations_per_simulation_thresholds.png"), plot, width = 1.5, height = 2, units = "in", dpi = 300)
}else{
  ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/tp53_mutations_per_simulation_thresholds.png"), plot, width = 1.5, height = 2, units = "in", dpi = 300)
}


############
##### Fill missing values
############

reps_num = 100
reps_num = 300 # use if combining all experiments

### Determine the missing replicates to be added.
summary_data <- alldat %>%
  filter(!(year < 1)) %>%
  filter(!(corrBlock == 0 & corrTime == 1)) %>%
  group_by(corrBlock, year, replicate) %>%
  summarize(observations = n())
#summary_data %>% print(n=1000)

years_total = unique(summary_data$year)

max_replicates = n_distinct(unique(summary_data$replicate[summary_data$corrBlock ==0]))

missing_replicates <- lapply(years_total, function(year_value) {
  summary_data %>%
    filter(year == year_value) %>%
    group_by(corrBlock) %>%
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
  dplyr::select(tp53Clone, replicate, corrBlock, year, mut) %>%
  filter(mut %in% c(68, 777)) %>% # tp53 + FA or tp53 + correction
  #filter(mut ==68) %>% # tp53 + FA only (no correction)
  group_by(tp53Clone, replicate, corrBlock, year, mut) %>%
  summarise(cells_count = n())
tp53clones %>% print(n=Inf)
tp53clones

missing_replicates
for (i in 1:nrow(missing_replicates)){
  n_missing <- missing_replicates$missing_replicates[i]
  corrBlock_replace <- missing_replicates$corrBlock[i]
  year <- missing_replicates$year[i]
  for (j in 1:n_missing){
    n_reps = reps_num
    random_rep = sample(0:n_reps-1, size = 1)
    sampled_rep <- tp53clones %>%
      filter(corrBlock==0, year == year, replicate == random_rep)
    sampled_rep <- mutate(sampled_rep, corrBlock = corrBlock_replace)

    tp53clones <- bind_rows(tp53clones, sampled_rep)
    
  }
}
tp53clones

unique_reps <- tp53clones %>%
  group_by(corrBlock, year) %>%
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
completed_reps <- mutate(completed_reps, completed_replicates = reps_num)


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
  group_by(corrBlock, year, replicate, threshold) %>%
  summarise(unique_count = n_distinct(tp53Clone)) %>%
  ungroup() %>%
  group_by(corrBlock, year, threshold) %>%
  summarise(total_unique_count = sum(unique_count)) %>%
  ungroup()
unique_tp53Clone %>% print(n=Inf)

## remove columns from completed dataframe
completed_reps <- completed_reps %>%
  dplyr::select(-mutRate, -corrTime)
completed_reps %>% print(n=Inf)

## combine tp53 mut counts with the number of replicates completed in experiment
unique_tp53Clone_reps <- completed_reps %>%
  left_join(unique_tp53Clone, by = c("corrBlock", "year"))

## average mutations per replicate
ave_unique_count <- unique_tp53Clone_reps %>%
  mutate(ave_unique_count = total_unique_count/completed_replicates) %>%
  drop_na()
print(ave_unique_count)
ave_unique_count <- ave_unique_count %>%
  filter(threshold == 108)

#### Add 0's to data to show up on plot even if no mutations
#### must manually check data first if 0s can be added.
start_years <- c(1, 6, 11)
start_years <- c(1)
unique_corrBlock <- unique(ave_unique_count$corrBlock)
unique_threshold <- unique(ave_unique_count$threshold)

 # Create a dataframe with the additional rows
additional_rows <- expand.grid(year = start_years,
                                ave_unique_count = 0,
                                corrBlock = unique_corrBlock,
                                threshold = unique_threshold)

#Combine with the original dataframe
ave_unique_count_extended <- bind_rows(ave_unique_count, additional_rows) %>%
  arrange(year, corrBlock, threshold)
print(ave_unique_count_extended, n=Inf)

## Plot using ggplot
num_colors <- length(unique(ave_unique_count$corrBlock))
my_palette <- rev(viridis(num_colors))
my_palette[1] <- "#BADE28FF"


custom_labels <- c("0" = "0", "0.01" = "0.01", "0.1" = "0.1")
green_shades <- c("0.1" = "#44aa99", "0.01" = "#67c2b3", "0" = "#bbe4dd")

plot <- ggplot(ave_unique_count_extended, aes(x = year, y = ave_unique_count, color = factor(corrBlock))) +
  geom_line(size = 1.5) +
  ylim(0, 0.45) +
  #labs(x = "Year", y = "TP53 mutation frequency", color = "p") +
  labs(x = "Year", y = expression(paste(italic("TP53"), " mutation frequency")), color = "p") +
  theme_minimal() + 
  scale_color_manual(values = green_shades, 
                     breaks = unique(ave_unique_count$corrBlock), 
                     labels = custom_labels) +
  #theme(aspect.ratio = 0.75) + 
  scale_x_continuous(breaks = seq(1, 50, 10)) +
  facet_wrap(~threshold, scales = "fixed", nrow = 1) +
  #theme(legend.position = "top") +
  theme(legend.position = c(0.35, 0.75)) +
  theme(
    axis.line = element_line(color="black"), 
    axis.text = element_text(color="black", size=8),
    axis.title.x = element_text(color="black", size=8),
    axis.title.y = element_text(color="black", size=8),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(color="black", size=8), 
    legend.title = element_text(color="black", face="italic", size=8),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
  )


show(plot)

## save plot to file
if(cluster){
  ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/tp53_mutations_per_simulation_thresholds_color.png"), plot, width = 1.5, height = 2, units = "in", dpi = 300)
}else{
  ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/tp53_mutations_per_simulation_thresholds_poster.png"), plot, width = 10, height = 6, units = "in", dpi = 300)
}


