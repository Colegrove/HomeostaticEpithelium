## FA dynamics
## clone size distribution plot
## 18 Jan 2024


require(tidyverse)
require(khroma)
require(cowplot)
require(RColorBrewer)
require("jsonlite")


TIME_POINT = 10 ## year to plot


## read in data csv as tibble
datdir <- "/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-05-SingleInjection_50years_DispersionSpacing/dat/"
datdir <- "/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-05-SingleInjection_50years_DispersionSpacing/dat/"
file <- "cloneSizes.csv"
filepath <- paste(datdir,file, sep="")
outcomes <- tibble(read.csv(filepath))

## convert the list of clone sizes from chr to a list of ints
outcomes$cloneSizes <- lapply(outcomes$cloneSizes, function(x) fromJSON(x))
## compute the mean clone size and add as a column
outcomes <- outcomes %>%
    rowwise() %>%
    mutate(aveCloneSize = mean(as.numeric(cloneSizes)))
## filter df for only the time point for plotting
outcomes <- outcomes %>%
    filter(year == TIME_POINT)
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

outcomes <- outcomes %>%
    filter(blockProb == 0 | blockProb == 0.001 | blockProb == 0.01)

x_tick <- 200
plot <- outcomes %>%
  ggplot(aes(x = xVal, y = cloneCounts, color=as.factor(blockProb))) +
  geom_line(linewidth=1) +
  scale_color_manual("Blocking\nProbability", values = muted) +
  theme_classic() +
  scale_y_continuous(breaks = seq(0,100,10)) +
  scale_x_continuous(breaks = seq(0,max_value+x_tick,x_tick)) +
  theme(legend.position = c(.91,.67)) +
  ylim(0,100) +
  labs(x = "Clone size (x)",
       y = "Number of clones >= x",
       title = paste(TIME_POINT,"year"))


ggsave("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/results/cloneSize_neutral.png", plot = plot, height = 3.5, width = 6)





# ## Remove ongoing datapoints from tibble, but add a fake datapoint to plot the legend label on panel A
# outToPlot <- outcomes %>% filter(Confluence != "Ongoing") 
# outToPlot <- bind_rows(outToPlot, tibble(X = "", Block_prob = 1, tConfluence = -100, Confluence = "Ongoing"))
# ## Find the median Confluence or loss for a given blocking probability - used for plotting medians
# dataSummaries <- outToPlot %>% group_by(Confluence, Block_prob) %>% summarize(tConfluence = median(tConfluence))


# ## Panels
# conf_loss_panel <- outToPlot %>% 
#     mutate(Confluence = factor(Confluence,levels = c("Confluence", "Ongoing", "Loss"))) %>%
#     ggplot(aes(x = factor(Block_prob), 
#                y = tConfluence, 
#                group = Confluence, 
#                color = Confluence)) + 
#     geom_jitter(height = 0, width = 0.25, alpha = 0.2) +
# #    stat_summary(aes(color = Confluence), fun.data = 'med_cl_boot') +
#     labs(x = "Blocking probability (b)",
#          y = "Time (Years)") +
#     theme_classic() + 
#     geom_path(aes(color = Confluence), data = dataSummaries) + 
#     geom_point(aes(color = Confluence), data = dataSummaries, size = 3) +
#     theme(legend.position = c(.85,.85)) +
#      theme(legend.title=element_blank()) + scale_color_manual(values = outcome_cols) + coord_cartesian(ylim = c(0, 50))

# conf_loss_panel 





# outcome_counts <-  outcomes %>% group_by(Block_prob, Confluence) %>% 
#     summarize(count = n()) 

# outcome_bar_panel <- outcome_counts %>% ggplot(aes(x = factor(Block_prob), y = count, fill = factor(Confluence, levels = c("Confluence", "Ongoing", "Loss")))) + 
#     geom_bar(stat = "identity") + labs(x = "Blocking probability (b)",
#          y = "Count")  + theme_classic() + theme(legend.position = "none")  + scale_fill_manual(values = outcome_cols)


# f1_panA_B <- plot_grid(conf_loss_panel, outcome_bar_panel, ncol = 1, labels = c("A", "B"), rel_heights = c(1, 0.65))


# ggsave("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-05-SingleInjection_50years_blockingProb/results/figs/fig2_panelsA_B.png", plot = f1_panA_B, height = 5, width = 4)
