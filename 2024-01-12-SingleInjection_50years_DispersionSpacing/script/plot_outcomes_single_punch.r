## FA dynamics
## Figure 2 plot
## 18 Jan 2024


require(tidyverse)
require(cowplot)
require(RColorBrewer)

## read in data csv as tibble
outcomes <- tibble(read.csv("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/dat/timeToConfluence.csv"))

### Facet wrap data into 3 separate plots
# outcomes %>%
#     ggplot(aes(x = factor(Block_prob), 
#                y = tConfluence, group = Confluence)) + 
#     geom_jitter(height = 0, width = 0.25, alpha = 0.5) +
#     stat_summary(fun.data = 'mean_cl_boot', color = "red") +
#     facet_wrap(~Confluence, ncol = 1, scales = "free_y") +
#     labs(x = "Blocking probability (b)",
#          y = "Time (Years)") +
#     theme_classic() + scale_y_log10()

#interquartile range
    ## geom_pointrange(
    ##     stat = "summary",
    ##     fun.min = function(z) {quantile(z,0.25)},
    ##     fun.max = function(z) {quantile(z,0.75)},
##     fun = median, color = "blue") +

## Create a color palette associated with category types
outcome_cols <- brewer.pal(5, "Set1")[c(1, 2, 5)]
    names(outcome_cols) <- c("Confluence", "Loss", "Ongoing")


## Remove ongoing datapoints from tibble, but add a fake datapoint to plot the legend label on panel A
outToPlot <- outcomes %>% filter(Confluence != "Ongoing") 
outToPlot <- bind_rows(outToPlot, tibble(X = "", Block_prob = 1, tConfluence = -100, Confluence = "Ongoing"))
## Find the median Confluence or loss for a given blocking probability - used for plotting medians
dataSummaries <- outToPlot %>% group_by(Confluence, Block_prob) %>% summarize(tConfluence = median(tConfluence))


## Panels
conf_loss_panel <- outToPlot %>% 
    mutate(Confluence = factor(Confluence,levels = c("Confluence", "Ongoing", "Loss"))) %>%
    ggplot(aes(x = factor(Block_prob), 
               y = tConfluence, 
               group = Confluence, 
               color = Confluence)) + 
    geom_jitter(height = 0, width = 0.25, alpha = 0.2) +
#    stat_summary(aes(color = Confluence), fun.data = 'med_cl_boot') +
    labs(x = "Blocking probability (b)",
         y = "Time (Years)") +
    theme_classic() + 
    geom_path(aes(color = Confluence), data = dataSummaries) + 
    geom_point(aes(color = Confluence), data = dataSummaries, size = 3) +
    theme(legend.position = c(.85,.85)) +
     theme(legend.title=element_blank()) + scale_color_manual(values = outcome_cols) + coord_cartesian(ylim = c(0, 50))

conf_loss_panel 





outcome_counts <-  outcomes %>% group_by(Block_prob, Confluence) %>% 
    summarize(count = n()) 

outcome_bar_panel <- outcome_counts %>% ggplot(aes(x = factor(Block_prob), y = count, fill = factor(Confluence, levels = c("Confluence", "Ongoing", "Loss")))) + 
    geom_bar(stat = "identity") + labs(x = "Blocking probability (b)",
         y = "Count")  + theme_classic() + theme(legend.position = "none")  + scale_fill_manual(values = outcome_cols)


f1_panA_B <- plot_grid(conf_loss_panel, outcome_bar_panel, ncol = 1, labels = c("A", "B"), rel_heights = c(1, 0.65))


ggsave("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/results/figs/fig2_panelsA_B.png", plot = f1_panA_B, height = 5, width = 4)
