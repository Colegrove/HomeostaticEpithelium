## FA dynamics
## Figure 2 plot
## 18 Jan 2024


require(tidyverse)
require(cowplot)
require(RColorBrewer)
require(viridis)

## read in data csv as tibble
outcomes <- tibble(read.csv("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/dat/timeToConfluence.csv"))


proportionSummaryLoss <- outcomes %>% 
    select(-tConfluence) %>% 
    group_by(Block_prob, Sigma, Dose) %>% 
    summarize(pLoss = sum(Confluence == "Loss") /n())

dataSummaries <- outcomes %>% group_by(Confluence, Block_prob, Sigma, Dose) %>% summarize(tConfluence = median(tConfluence))

dataSummariesConfluence <- dataSummaries %>% 
    filter(Confluence == "Confluence")

dataSummariesLoss <- dataSummaries %>% 
    filter(Confluence == "Loss")

heatmapConfluence <- dataSummariesConfluence %>%
  ggplot(aes(x = as.factor(Sigma), y = as.factor(Block_prob), fill = tConfluence)) +
  geom_tile() +
  facet_wrap(~ Dose, scales = "free", ncol = 3, drop=FALSE) + 
  theme_classic() +
  labs(title = "dose", x = "Sigma", y="Block_prob") +
  geom_text(aes(label = round(tConfluence, 1))) +
  #scale_fill_viridis_c()
  scale_fill_viridis_c(limits = c(0, 50), oob = scales::squish)

ggsave("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/results/tConfluence_heatmap.png", plot = heatmapConfluence, height = 4, width = 10)

heatmapLoss <- proportionSummaryLoss %>%
  ggplot(aes(x = as.factor(Sigma), y = as.factor(Block_prob), fill = pLoss)) +
  geom_tile() +
  facet_wrap(~ Dose, scales = "free", ncol = 3, drop=FALSE) +
  theme_classic() +
  labs(title = "dose", x = "Sigma", y="Block_prob", fill="pLoss") + 
  geom_text(aes(label = round(pLoss, 2))) +
  scale_fill_viridis_c()
  #scale_fill_viridis_c(limits = c(0, 5), oob = scales::squish)

ggsave("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/results/pLoss_heatmap.png", plot = heatmapLoss, height = 4, width = 10)



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


# ggsave("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/results/figs/tConfluence_heatmap.png", plot = f1_panA_B, height = 5, width = 4)
