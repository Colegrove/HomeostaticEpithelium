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
project_path = "2024-07-04-tp53_competition_sameMut_sameP"
sigmas = c(2)
doses = c(30)
block_values = c(0, 0.01, 0.1)
p53block = 0.01
corrTimes = c(1, (50+1)*4.5*364)
nreplicates = 100
#nreplicates = 3
replicate = 0:nreplicates
mutRates = c(4)

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


## Timesteps to be selected
#timesteps = c(1*364*4.5, 11*364*4.5, 20*364*4.5, 30*364*4.5, 40*364*4.5, 50*364*4.5)
#timesteps <- seq(1, 100) * 182 * 4.5
timesteps <- c(2, seq(1, 50, 5) * 364 * 4.5)

## Path to directory of raw data visualization files
if(cluster){
  path = glue("/net/feder/vol1/project/HomeostaticEpithelium/dat/{project_path}/results/")
}else{
  path = glue("/Volumes/feder-vol1/project/HomeostaticEpithelium/dat/{project_path}/results/")
}


## combines all data into a large tibble
alldat <- foreach(block = block_values, .combine = "rbind") %:%
            foreach(rep = replicate, .combine = "rbind") %:%
              foreach(corrTime = corrTimes, .combine = "rbind") %:%
                foreach(timestep = timesteps, .combine = "rbind") %do%  {

  f <- glue("VisFile_block_{block}_growth_1.0_sigma_{sigmas[1]}_dose_{doses[1]}_p53_{p53block}_mutrate_{mutRates[1]}_correctionTime_{corrTime}_{rep}.{timestep}.txt")
  cat("Opening file:", f, "\n")
  
  full_path <- paste0(path, f)
  print(full_path)
  column_names = c('x','z','y', 'h','s','v','alpha','mut','injSite','tp53Clone')
  if (file.exists(full_path)){
    data <- tibble(read.csv(full_path, sep = '\t', header = FALSE, col.names = column_names)) %>% 
      filter(y == 0) %>% # basal only
      mutate(replicate = rep) %>%
      mutate(corrBlock = block) %>%
      #mutate(tp53Block = p53_val) %>%
      mutate(year = timestep/4.5/364) %>%
      # mutate(sigma = sigma) %>%
      # mutate(dose = dose)
      mutate(mutRate = mutRates[1]) %>%
      mutate(corrTime = corrTime)
  }
}


## output all basal layer positional data
if(cluster){
  #write.csv(alldat, glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/position_data_All.csv"), row.names = FALSE)
  write.csv(alldat, glue("/net/feder/vol1/project/HomeostaticEpithelium/dat/{project_path}/position_data_All.csv"), row.names = FALSE)
}else{
  write.csv(alldat, glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/position_data_All.csv"), row.names = FALSE)
}