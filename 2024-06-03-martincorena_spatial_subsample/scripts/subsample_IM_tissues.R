## Hunter Colegrove
## FA dynamics
## 24 May 2024

## Subsample spatially reconstructed IM tissues 
## First, generate reconstructed tissues using assignSpatialMutations.R



library(tidyverse)
library(glue)
library(ggplot2)
library(RColorBrewer)

project_path = "2024-06-03-martincorena_spatial_subsample"

#####################################

## Set the seed
set.seed(11)

#filenames <- list.files("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-04-22-tp53_martincorena_54years/dat", pattern="*.csv", full.names=TRUE)
filenames <- list.files("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{2024-06-03-martincorena_spatial_subsample}/dat/reconstructed_tissue", pattern="*.csv", full.names=TRUE)


## Grid size and subsample size
grid_size <- 173
sample_size <- 70
cell_density <- 15000
mutation_threshold <- 108 # (0.0018 *2*2*15000) 0.0018 smallest detected vaf in any gene

final_df <- tibble(sampleID = character(), mut = numeric(), vaf = numeric(), age = character(), size = numeric(), area = numeric())

for(file in filenames){
  df <- read_csv(file)
  split <- unlist(strsplit(file, split = "_"))
  pid <- split[7]
  age <- split[8]
  age <- unlist(strsplit(age, split = ".csv"))
  
  ## Randomly select the starting point for the subsection
  start_x <- sample(1:(grid_size - sample_size + 1), 1)
  start_y <- sample(1:(grid_size - sample_size + 1), 1)
  
  # Define the ending points for the subsection
  end_x <- start_x + sample_size - 1
  end_y <- start_y + sample_size - 1
  
  # Filter out the rest of spaces not in the subsection
  sampled_df <- df %>%
    filter(x >= start_x & x <= end_x & y >= start_y & y <= end_y)

  sample_vaf_df <- sampled_df %>%
    filter(Mut1 != 0) %>%
    group_by(Mut1) %>%
    summarize(count = n()) %>%
    mutate(vaf = count / (2*(sample_size*sample_size/cell_density)*cell_density),
           sampleID = pid,
           mut = Mut1, 
           age = age,
           size = count, 
           area = count/cell_density) %>%
    select(sampleID, mut, vaf, age, size, area)
  
  # Append the results to the final dataframe
  final_df <- bind_rows(final_df, sample_vaf_df)
  
  table(sampled_df$Mut1)
  
  ##############
  #### Visualize
  ##############
  

  df <- df %>%
    mutate(Mut1 = as.factor(df$Mut1)) %>%
    mutate(Mut2 = as.factor(df$Mut2))

  sampled_df <- sampled_df %>%
    mutate(Mut1 = as.factor(sampled_df$Mut1)) %>%
    mutate(Mut2 = as.factor(sampled_df$Mut2))

  p1 <- ggplot(df, aes(x = x, y = y, fill = Mut1)) +
    geom_tile() +
    scale_fill_brewer(palette = "Set3") +
    labs(title = "Mutation Distribution on 173x173 Grid",
         x = "Column",
         y = "Row",
         fill = "Mutation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    coord_fixed() +
    geom_rect(aes(xmin = start_x - 0.5, xmax = end_x + 0.5, ymin = start_y - 0.5, ymax = end_y + 0.5),
              color = "red", fill = NA, linewidth = 1)

  p2 <- ggplot(sampled_df, aes(x = x, y = y, fill = Mut1)) +
    geom_tile() +
    scale_fill_brewer(palette = "Set3") +
    labs(title = "Mutation Distribution on 173x173 Grid",
         x = "Column",
         y = "Row",
         fill = "Mutation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    coord_fixed()

  #ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/subsamples/{pid}_{age}_sample_cutout_plot.png"), p1, width = 10, height = 5, units = "in", dpi = 300)
  #ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/subsamples/{pid}_{age}_subsample_plot.png"), p2, width = 10, height = 5, units = "in", dpi = 300)

  ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/subsamples/{pid}_{age}_sample_cutout_plot.png"), p1, width = 10, height = 5, units = "in", dpi = 300)
  ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/subsamples/{pid}_{age}_subsample_plot.png"), p2, width = 10, height = 5, units = "in", dpi = 300)

}

# adjust for the lower limit of detection
final_df <- final_df %>% filter(size >= mutation_threshold)
write.csv(final_df, glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/subsamples/mutations_subsamples.csv"), row.names = FALSE)

###################################