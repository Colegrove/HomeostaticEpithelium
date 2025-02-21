## FA dynamics
## clone size distribution plot
## 18 Jan 2024

## This script uses cloneSizes.csv data generated from clonesSizes.py
## This script plots the clone size distribution at a given timepoint across
## the range of persistence coefficient values.

require(tidyverse)
require(khroma)
require(cowplot)
require(RColorBrewer)
require("jsonlite")
require(glue)

TIME_POINT = 1 ## year to plot

project_path = "2024-04-10-SingleInjection_Confluence_blockProb"
## read in data csv as tibble
#datdir <- glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/")
datdir <- glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/")
file <- "cloneSizes.csv"
filepath <- paste(datdir,file, sep="")
outcomes <- tibble(read.csv(filepath))

## convert the list of clone sizes from chr to a list of ints
outcomes$cloneSizes <- lapply(outcomes$cloneSizes, function(x) fromJSON(x))

## compute the mean clone size and add as a column
outcomes <- outcomes %>%
    rowwise() %>%
    mutate(aveCloneSize = mean(as.numeric(cloneSizes)))
outcomes
## filter df for only the time point for plotting
outcomes <- outcomes %>%
    filter(year == TIME_POINT)
outcomes
## find the maximum clone size and add 1 to ensure an x-intercept
max_value <- max(sapply(outcomes$cloneSizes, max)) + 1

## calculate how many clones are greater or equal to x up to the max value.
outcomes <- outcomes %>%
    rowwise() %>%
    mutate(cloneCounts = list(sapply(1:max_value, function(x) sum(as.numeric(cloneSizes) >=x))))


## explode cloneSizes list into long form
outcomes <- outcomes %>%
    unnest(cloneCounts)
## add column with associated x-axis value
outcomes <- outcomes %>% 
    group_by(blockProb) %>%
    mutate(xVal = row_number() - 1) %>%
    ungroup()

## Create a color palette associated with category types

muted <- color("muted")(9)[c(1:3, 5:9)] # skip green, keep teal
    names(muted) <- c("0", "0.001", "0.01", "0.1", "0.2", "0.5", "1")

# outcomes <- outcomes %>%
#     filter(blockProb == 0 | blockProb == 0.001 | blockProb == 0.01)

x_tick <- 200
plot <- outcomes %>%
  ggplot(aes(x = xVal, y = cloneCounts, color=as.factor(blockProb))) +
  geom_line(linewidth=1) +
  scale_color_manual("Persistence\ncoefficient", values = muted) +
  theme_classic() +
  scale_y_continuous(breaks = seq(0,100,10)) +
  scale_x_continuous(breaks = seq(0,max_value+x_tick*2,x_tick), limits = c(0, max_value + x_tick)) +
  theme(legend.position = c(.9,.67)) +
  labs(x = "Clone size (x)",
       y = "Number of clones >= x",
       title = paste("year:", TIME_POINT))

## save individual plot
ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/cloneSize_year_{TIME_POINT}.png"), plot = plot, height = 3.5, width = 6)