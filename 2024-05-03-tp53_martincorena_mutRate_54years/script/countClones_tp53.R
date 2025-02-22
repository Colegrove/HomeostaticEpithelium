## Hunter Colegrove
## FA dynamics
## 
## Outputs cloneSizes_tp53.csv to track number of mutations

## Update parameters from snakefile prior to running

require(foreach)
require(tidyverse)
require(glue)

cluster = TRUE ## change to FALSE if running locally
#clone_detection_thresholds = c(1, 5, 10, 20, 50, 100, 200, 300, 400, 500)
#clone_detection_thresholds = c(10,88,100,400)
clone_detection_thresholds = c(108)


## Project information and parameters used in experiment (enter from snakefile)
project_path = "2024-05-03-tp53_martincorena_mutRate_54years"
sigma = c(2)
dose = c(30)
block_values = c(0) # corrected block value
p53_block_values = c(0, 0.007, 0.01, 0.0125, 0.015)
corrTime = c(55*4.5*364)
nreplicates = 235
replicate = 0:nreplicates
mutRates = c(1,4)



### adjust p53 block values to fit filename outputs if necessary
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


## Timesteps to be selected
timesteps = c(35217, 41769, 61425, 74529, 81081, 87633) # martincorena data time ranges

## Path to directory of raw data visualization files
if(cluster){
  path = glue("/net/feder/vol1/project/HomeostaticEpithelium/dat/{project_path}/results/")
} else {
  path = glue("/Volumes/feder-vol1/project/HomeostaticEpithelium/dat/{project_path}/results/")
}


## combines all data into a large tibble
alldat <- foreach(p53_val = p53_block_values, .combine = "rbind") %:%
            foreach(rep = replicate, .combine = "rbind") %:%
              foreach(timestep = timesteps, .combine = "rbind") %:%
                foreach(mutRate = mutRates, .combine = "rbind") %do% {
  f <- glue("VisFile_block_{block_values[1]}_growth_1.0_sigma_{sigma[1]}_dose_{dose[1]}_p53_{p53_val}_mutrate_{mutRate}_correctionTime_{corrTime[1]}_{rep}.{timestep}.txt")
  cat("Opening file:", f, "\n")
  #paste("Opening file: fakedat/", f)
  full_path <- paste0(path, f)
  print(full_path)
  column_names = c('x','z','y', 'h','s','v','alpha','mut','injSite','tp53Clone')
  if (file.exists(full_path)){
    data <- tibble(read.csv(full_path, sep = '\t', header = FALSE, col.names = column_names)) %>% 
      filter(y == 0) %>%
      mutate(replicate = rep) %>%
      mutate(tp53Block = p53_val) %>%
      mutate(year = timestep/4.5/364) %>%
      mutate(mutRate = mutRate)
      #filter(tp53Clone != -1)
  }
}

## output all positional data from all runs
if(cluster){
  write.csv(alldat, glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/position_data_All.csv"), row.names = FALSE)
} else {
  write.csv(alldat, glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/position_data_All.csv"), row.names = FALSE)
}


##########
#### Number of tp53 mutations
##########

## find n replicates for each tp53Block/year
## this number may be variable if data analyzed mid run. 
completed_reps <- alldat %>%
  group_by(tp53Block, year, mutRate) %>%
  summarise(completed_replicates = n_distinct(replicate))

## filter for tp53 clones greater than a detection threshold
filter_threshold <- function(data, threshold){
  data %>%
    group_by(tp53Block, replicate, tp53Clone, year, mutRate) %>%
    filter(tp53Clone != -1) %>% # exclude non-tp53 mutations
    summarize(count = n()) %>%
    filter(count >= threshold) %>% # apply detection threshold
    group_by(tp53Block, replicate, year, mutRate) %>%
    summarize(unique_tp53Clones = n()) %>%
    group_by(tp53Block, year, mutRate) %>%
    summarize(total_unique_tp53Clones = sum(unique_tp53Clones)) %>%
    mutate(threshold = threshold)
}

## look across a range of detection thresholds and combine results
result <- map_dfr(clone_detection_thresholds, ~filter_threshold(alldat, .x))

## combine tp53 mut counts with the number of replicates completed in experiment
result <- result %>%
  left_join(completed_reps %>% select(tp53Block, year, mutRate, completed_replicates), 
            by = c("tp53Block", "year", "mutRate"))
## calculate the frequency of mutations per tissue section (replicate)
result <- result %>%
  mutate(tp53clone_frequency = total_unique_tp53Clones / completed_replicates)

## output clone count information
result <- result %>%
  select(-total_unique_tp53Clones)
result <- result %>%
  rename(blockProb = tp53Block) %>%
  rename(age = year) %>%
  rename(AverageMutationCount = tp53clone_frequency)
result

if(cluster){
  write.csv(result, glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/countClones_tp53_R.csv"), row.names = FALSE)
} else {
  write.csv(result, glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/countClones_tp53_R.csv"), row.names = FALSE)
}



##########
#### VAF calculations
##########

## VAF

resultVAF <- alldat %>%
  group_by(tp53Block, replicate, tp53Clone, year, mutRate) %>%
  filter(tp53Clone != -1) %>%
  summarize(count = n()) %>%
  print() %>%
  filter(count >= 88) %>%
  mutate(VAF = count/(70*70*2)) %>%
  print(n=Inf)

resultVAF <- resultVAF %>%
  rename(blockProb = tp53Block) %>%
  rename(age = year)


if(cluster){
  write.csv(resultVAF, glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/VAF_tp53_R.csv"), row.names = FALSE)
} else {
  write.csv(resultVAF, glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/VAF_tp53_R.csv"), row.names = FALSE)
}



##########
#### Single tissue TestOuput
##########

# #path = glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/TestOutput/VisFiles/")
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
#### Grid tissue TestOuput using alldat tibble
##########

# alldat
# write.csv(alldat, glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/position_data_All.csv"), row.names = FALSE)
# 
# tempdat <- alldat %>%
#   filter(replicate %in% c(20:30)) %>% ## filter for vis replicates
#   filter(year == 45.5) %>%
#   filter(tp53Block == 0.005) %>%
#   filter(mutRate == 4) %>%
#   filter(tp53Clone > -1) %>%
#   mutate(replicate = as.factor(replicate))
# print(tempdat, n=Inf)
# 
# # Calculate the count of each unique value
# value_counts <- tempdat %>%
#   filter(tp53Clone != -1) %>%
#   group_by(tp53Clone) %>%
#   summarize(count = n())
# print(value_counts)
# 
# # Create a color palette based on the count of unique values
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
  