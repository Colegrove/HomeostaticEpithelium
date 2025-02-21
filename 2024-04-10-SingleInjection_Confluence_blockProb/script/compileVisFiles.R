## Hunter Colegrove
## FA dynamics
##
## Reads in raw vis files, filters only basal layer and outputs to combined file


require(foreach)
require(tidyverse)
require(glue)
library(gridExtra)
require(grid)
library(imager)
library(gtools)
library(png)


cluster = TRUE

## Project information and parameters used in experiment (enter from snakefile)
project_path = "2024-04-10-SingleInjection_Confluence_blockProb"
sigmas = c(2)
doses = c(10)
block_values = c(0, 0.001, 0.01, 0.1, 0.2, 0.5, 1.0)
corrTimes = c(1)

nreplicates = 100
nreplicates = 2
replicate = 0:nreplicates


## Timesteps to be selected
timesteps <- c(2, seq(1, 50, 5) * 364 * 4.5)


## Path to directory of raw data visualization files
if(cluster){
  path = glue("/net/feder/vol1/project/HomeostaticEpithelium/dat/{project_path}/results/")
}else{
  path = glue("/Volumes/feder-vol1/project/HomeostaticEpithelium/dat/{project_path}/results/")
}

# Find unique timestpes
all_files <- list.files(path, pattern = "VisFile_.*\\.txt", full.names = TRUE)
# Extract unique timestep values from filenames
timesteps <- unique(as.numeric(sub(".*\\.(\\d+)\\.txt$", "\\1", all_files)))
print(timesteps)


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
      filter(y == 0) %>% # basal only
      mutate(replicate = rep) %>%
      mutate(corrBlock = block) %>%
      mutate(year = timestep/4.5/364)
  }
}


## output all basal layer positional data
if(cluster){
  #write.csv(alldat, glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/position_data_All.csv"), row.names = FALSE)
  write.csv(alldat, glue("/net/feder/vol1/project/HomeostaticEpithelium/dat/{project_path}/position_data_All.csv"), row.names = FALSE)
}else{
  write.csv(alldat, glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/position_data_All.csv"), row.names = FALSE)
}