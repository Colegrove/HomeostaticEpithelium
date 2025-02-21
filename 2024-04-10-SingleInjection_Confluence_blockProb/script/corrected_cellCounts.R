## Hunter Colegrove
## FA dynamics
## 29 Mar 2024

## This file compiles the cell counts at each time point.
## The output file corrCellCounts.csv can be used to run speedArea.R

require(glue)
require(foreach)
require(tidyverse)


## Project information and parameters used in experiment (enter from snakefile)
project_path = "2024-04-10-SingleInjection_Confluence_blockProb"
sigmas = c(2)
doses = c(10)
block_values = c(0, 0.001, 0.01, 0.1, 0.2, 0.5, 1.0)
corrTimes = c(1)
nreplicates = 100
replicate = 0:nreplicates

## Path to directory of raw data visualization files
path = glue("/net/feder/vol1/project/HomeostaticEpithelium/dat/{project_path}/results/")

## Timesteps to be selected 
## generic timesteps at each 6 month interval
#timesteps <- c(2, seq(1, 100) * 182 * 4.5)
## add unique timesteps associated with confluence/loss
all_files <- list.files(path)
vis_files <- all_files[str_starts(all_files, "VisFile")]

timepoints <- sapply(vis_files, function(f) {
  parts <- strsplit(f, "_")[[1]]
  last_part <- parts[length(parts)]
  last_split <- strsplit(last_part, "[.]")[[1]]
  timepoint <- last_split[2]
  return(timepoint)
})

timesteps <- sort(unique(as.integer(timepoints)))

## combines all data into a large tibble
alldat <- foreach(block = block_values, .combine = "rbind") %:%
  foreach(rep = replicate, .combine = "rbind") %:%
  foreach(timestep = timesteps, .combine = "rbind") %do%  {
    f <- glue("VisFile_block_{block}_growth_1.0_sigma_{sigmas[1]}_dose_{doses[1]}_correctionTime_{corrTimes[1]}_{rep}.{timestep}.txt")
    cat("Opening file:", f, "\n")

    full_path <- paste0(path, f)
    print(full_path)
    column_names = c('x','z','y', 'h','s','v','alpha','mut','injSite','tp53Clone')
    if (file.exists(full_path)){

      data <- tibble(read.csv(full_path, sep = '\t', header = FALSE, col.names = column_names)) %>% 
        filter(y == 0) %>%
        filter(injSite == 0) %>%
        mutate(replicate = rep) %>%
        mutate(corrBlock = block) %>%
        mutate(year = timestep/4.5/364) %>%
        
        ##### use for counting how many corrected cells script
        mutate(cellCount = n()) %>% 
        group_by(corrBlock, replicate, year) %>% 
        summarize(cellCount = unique(cellCount))
    }
  }

write.csv(alldat, glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/corrCellCounts.csv"), row.names = FALSE)

