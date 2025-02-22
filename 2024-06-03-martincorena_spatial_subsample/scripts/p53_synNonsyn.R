## Hunter Colegrove
## FA dynamics
## 
## Plot mutation frequency per over time and compare to Martincorena esophagus

## First, run subsectEsophagus.R to generate random subsections of martincorena data
## Run, countClones_tp53.R to generate tp53 clone sizes across a range of detection thresholds


cluster = FALSE
vis_mutRate = 4

require(tidyverse)
require(viridis)
require(glue)

## project path of subsampled reconstructed tissue
subsample_path = "2024-06-03-martincorena_spatial_subsample"

## project path of simulations
project_path = "2024-05-03-tp53_martincorena_mutRate_54years"

if(cluster){
  ## old subsect way
  #subsectfilepath = glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/subsect_clones_tp53.csv")
  
  ## new subsample way
  subsamplefilepath = glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{subsample_path}/dat/subsamples/mutations_subsamples.csv")
  
  ## tp53 unique clones across detection thresholds
  simfilepath = glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/countClones_tp53_R.csv")
  
  ## VAF path
  #vaffilepath = glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/VAF_tp53_R.csv")
}else{
  
  ## old subsect way
  #subsectfilepath = glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/subsect_clones_tp53.csv")
  
  ## new subsample way
  subsamplefilepath = glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{subsample_path}/dat/subsamples/mutations_subsamples.csv")
  
  ## tp53 unique clones across detection thresholds
  simfilepath = glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/countClones_tp53_R.csv")
  
  ## VAF path
  #vaffilepath = glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/VAF_tp53_R.csv")
}

#nsubsections = 9 ## for 70x70 tissue sections

## read in files
tp53sim = read_delim(simfilepath, delim = ",")
tp53subsample = read_delim(subsamplefilepath, delim = ",")
#tp53vafSim = read_delim(vaffilepath, delim = ",")


## number of samples collected for each age range pulled from Martincorena data
#nSamples <- c(84, 93, 96, 94, 95, 98, 95, 95, 94)
nSamples <- c(84, 93, 96, 94, 95, 98)

## calculate the average number of mutations per tissue sample
tp53_persample <- tp53subsample %>%
  group_by(age) %>%
  summarise(MutationCount = n()) %>%
  mutate(nSamples = nSamples) %>%
  mutate(AverageMutationCount = MutationCount / (nSamples)) %>%
  mutate(blockProb = "IM_data")

tp53_persample

## line plot of martincorena data. Mutation frequency per sample by age range.
ggplot(tp53_persample, aes(x = age, y = AverageMutationCount, group = 1)) +
  geom_line() +
  labs(x = "Age Range", y = "Mutation frequency per subsample") +
  theme_minimal()

#######################
##### simulation data

## adjust age ranges to match martincorena data
tp53simAge <- tp53sim %>%
  #rename(age = time) %>% ## toggle on if using python
  mutate(age = case_when(
    age == 21.5 ~ "20-23",
    age == 25.5 ~ "24-27",
    age == 37.5 ~ "36-39",
    age == 45.5 ~ "44-47",
    age == 49.5 ~ "48-51",
    age == 53.5 ~ "52-55",
    ## below age conversion uses schenk default time
    # age == 4.5 ~ "20-23",
    # age == 5.5 ~ "24-27",
    # age == 8.5 ~ "36-39",
    # age == 10 ~ "44-47",
    # age == 11 ~ "48-51",
    # age == 12 ~ "52-55",
    # age == 13 ~ "56-59",
    # age == 15.5 ~ "68-71",
    # age == 16.5 ~ "72-75",
    TRUE ~ as.character(age)  # Keep the original value if none of the conditions match
  )) %>%
  mutate(blockProb = as.character(blockProb))


## filter the mutation rate if multiple run
tp53simAge <- tp53simAge %>% 
  filter(mutRate == vis_mutRate)

## duplicate martincorena results to have same number as thresholds in simulations
threshold_vector = unique(tp53simAge$threshold)
df_list <- list()
for (val in threshold_vector) {
  duplicated_df <- tp53_persample
  duplicated_df$threshold <- val
  df_list <- append(df_list, list(duplicated_df))
}
tp53_persample <- do.call(rbind, df_list)

## combine the data between martincorena and simulated tissue
combined_tibble <- bind_rows(tp53simAge, tp53_persample)
combined_tibble <- combined_tibble %>%
  mutate(blockProb = as.numeric(blockProb))

## lineplot of average mutation counts per tissue section - martincorena and simulated tissue
plot <- ggplot(combined_tibble, aes(x = age, y = AverageMutationCount, group = as.factor(blockProb), color = as.factor(blockProb))) +
  geom_point(size = 5) +
  geom_line() +
  labs(x = "Age Range", y = "Mutations per Sample") +
  theme_minimal() +
  #facet_wrap(~threshold, scales = "free", nrow = 2) +
  theme(strip.text = element_text(size = 20))

## save plot to file
if(cluster){
  ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{subsample_path}/dat/martincorena_comparison_thresholds_mutRate_{vis_mutRate}.png"), plot, width = 12, height = 6, units = "in", dpi = 300)
}else{
  ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{subsample_path}/dat/martincorena_comparison_thresholds.png"), plot, width = 12, height = 6, units = "in", dpi = 300)
}


############### square difference