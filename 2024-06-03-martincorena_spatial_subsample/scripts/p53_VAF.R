## Hunter Colegrove
## FA dynamics
## 
## Plot mutation frequency per over time and compare to Martincorena esophagus

## First, run subsectEsophagus.R to generate random subsections of martincorena data
## Run, countClones_tp53.R to generate tp53 clone sizes across a range of detection thresholds


cluster = TRUE
vis_mutRate = 4
lower_detection_limit = 108 # cells

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
  vaffilepath = glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/VAF_tp53_R.csv")
}else{
  ## old subsect way
  #subsectfilepath = glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/subsect_clones_tp53.csv")
  ## new subsample way
  subsamplefilepath = glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{subsample_path}/dat/subsamples/mutations_subsamples.csv")
  ## tp53 unique clones across detection thresholds
  simfilepath = glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/countClones_tp53_R.csv")
  ## VAF path
  vaffilepath = glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/VAF_tp53_R.csv")
}


## read in files
#tp53sim = read_delim(simfilepath, delim = ",")
tp53subsample = read_delim(subsamplefilepath, delim = ",")
tp53vafSim = read_delim(vaffilepath, delim = ",")

########
### VAF of martincorena SNVs
########

#### plot martincorena vaf values only from subsected samples
# ggplot(tp53subsect, aes(x = age, y = vaf)) +
#   geom_point() + #position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 2) +
#   labs(x = "Age",
#        y = "VAF") +
#   theme_minimal()
# print(tp53subsect)

######### VAF combined simulations and martincorena 

## re-assign age ranges based on martincorena data
tp53vafSim <- tp53vafSim %>%
  mutate(age = case_when(
    age == 21.5 ~ "20-23",
    age == 25.5 ~ "24-27",
    age == 37.5 ~ "36-39",
    age == 45.5 ~ "44-47",
    age == 49.5 ~ "48-51",
    age == 53.5 ~ "52-55",
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
  ))


## only grab needed columns and assign IM label to martincorena data
tp53vaf_subsample = tp53subsample %>%
  select(vaf, age, size) %>%
  mutate(blockProb = "IM_data") %>%
  rename(count = size)

## order simulation data from weakest to strongest block and filter LLOD
tp53vafSim <- tp53vafSim %>%
  arrange(desc(blockProb)) %>%
  mutate(blockProb = as.factor(blockProb)) %>%
  rename(vaf = VAF) %>%
  filter(count >= lower_detection_limit)

## combine simulation data and IM data
tp53_allVAF <- bind_rows(tp53vafSim, tp53vaf_subsample)

### plot VAF of simulations and martincorena
factor_levels <- c("0.005","0.007", "0.01")

# Define custom color palette using a colorblind-friendly palette
my_palette <- c("IM_data" = "black",
                setNames(c("#DCE319FF", "#1F968DFF", "#88CCEE"), factor_levels))

tp53_allVAF <- tp53_allVAF %>%
  filter(!(blockProb == 0.0125 | blockProb == 0.015))

## plot
plot <- ggplot(tp53_allVAF, aes(x = age, y = vaf, color = blockProb)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.01, dodge.width = .9), size = 2) +
  labs(x = "Age",
       y = "VAF",
       color = "TP53 Block Probability") +
  scale_color_manual(values = my_palette) +
                     #breaks = c(NA, factor_levels),  # Set the breaks in numerical order
                     #labels = c("IM_data", factor_levels)) +  # Set the labels in numerical order
  theme_minimal() +
  theme(aspect.ratio = 0.5)


## save plot to file
if(cluster){
  ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{subsample_path}/dat/martincorena_comparison_thresholds_VAF.png"), plot, width = 10, height = 5, units = "in", dpi = 300)
}else{
  ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{subsample_path}/dat/martincorena_comparison_thresholds_VAF.png"), plot, width = 10, height = 5, units = "in", dpi = 300)
}
