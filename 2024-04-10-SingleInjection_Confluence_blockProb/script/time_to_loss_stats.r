## Hunter Colegrove
## FA dynamics
## Figure 2 plot
## 18 Jan 2024

## This script uses the timeToConfluence.csv data generated from timeToConfluence.py
## This script finds the proportion of lost clones that were lost within 1 year

require(tidyverse)
require(cowplot)
require(RColorBrewer)
require(glue)
require(jsonlite)
require(scales)


## project
project_path = "2024-04-10-SingleInjection_Confluence_blockProb"
## read in data csv as tibble
#outcomes <- tibble(read.csv(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/timeToConfluence.csv")))
outcomes <- tibble(read.csv(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/timeToConfluence.csv")))


loss_proportion_by_block <- outcomes %>%
  group_by(Block_prob) %>%
  summarise(
    total_loss = sum(Confluence == "Loss"),  # Total Loss events for this Block_prob
    loss_within_tConfluence = sum(Confluence == "Loss" & tConfluence <= 1),  # Loss events within tConfluence <= 1
    proportion_loss = loss_within_tConfluence / total_loss  # Proportion of Loss within tConfluence <= 1
  ) %>%
  ungroup() 


loss_proportion_by_block

total_loss_within_tConfluence <- sum(loss_proportion_by_block$loss_within_tConfluence)
total_loss <- sum(loss_proportion_by_block$total_loss)

total_proportion_loss <- total_loss_within_tConfluence / total_loss
total_proportion_loss




##################################################################

##################################################################

##################################################################

# ## Create a color palette associated with category types
# outcome_cols <- brewer.pal(5, "Set1")[c(1, 2, 5)]
# names(outcome_cols) <- c("Confluence", "Loss", "Ongoing")
# 
# ## Remove ongoing datapoints from tibble, but add a decoy datapoint to plot the legend label on panel A
# outToPlot <- outcomes %>% filter(Confluence != "Ongoing") 
# outToPlot <- bind_rows(outToPlot, tibble(X = "", Block_prob = 1, tConfluence = -100, Confluence = "Ongoing"))


# ## Find the median Confluence or loss for a given blocking probability - used for plotting medians
# dataSummaries <- outToPlot %>% group_by(Confluence, Block_prob) %>% summarize(tConfluence = stats::median(tConfluence))
# 
# 
# 
# ## Panels
# conf_loss_panel <- outToPlot %>% 
#     mutate(Confluence = factor(Confluence,levels = c("Confluence", "Ongoing", "Loss"))) %>%
#     ggplot(aes(x = factor(Block_prob), 
#                y = tConfluence, 
#                group = Confluence, 
#                color = Confluence)) + 
#     geom_jitter(height = 0, width = 0.25, alpha = 0.2) +
# #    stat_summary(aes(color = Confluence), fun.data = 'med_cl_boot') +
#     labs(x = "Persistence coefficient (p)",
#          y = "Years") +
#     theme_classic() + 
#     geom_path(aes(color = Confluence), data = dataSummaries) + 
#     geom_point(aes(color = Confluence), data = dataSummaries, size = 2) +
#     theme(legend.position = c(.85,.85)) +
#     theme(legend.title=element_blank()) + 
#     scale_color_manual(values = outcome_cols) + 
#     coord_cartesian(ylim = c(0, 50)) + 
#     theme(text=element_text(size=8), 
#         legend.background = element_blank(),
#         legend.text=element_text(size=8),
#         legend.key.size = unit(0.5, "lines"),
#         legend.key.height = unit(0.5, "lines"),
#         legend.margin = margin(0, 0, 0, 0),
#         legend.key.spacing.y = unit(0.0001, "cm"),
#         legend.key.spacing.x = unit(0.0001, "cm"),
#         axis.title = element_text(size = 8),
#         axis.text = element_text(size = 8, color = 'black'))
# 
# outcome_counts <-  outcomes %>% group_by(Block_prob, Confluence) %>% 
#     summarize(count = n()) 
# 
# outcome_bar_panel <- outcome_counts %>% ggplot(aes(x = factor(Block_prob), y = count, fill = factor(Confluence, levels = c("Confluence", "Ongoing", "Loss")))) + 
#     geom_bar(stat = "identity") + labs(x = "Persistence coefficient (p)",
#          y = "Count")  + theme_classic() + theme(legend.position = "none")  + scale_fill_manual(values = outcome_cols) + 
#   theme(text=element_text(size=8),
#         axis.title = element_text(size = 8),
#         axis.text = element_text(size = 8, color = 'black'))
# 
# f2_panA_B <- plot_grid(outcome_bar_panel, conf_loss_panel, ncol = 1, rel_heights = c(0.6, 1))
# show(f2_panA_B)
# 
# ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/fig2_panelsA.png"), plot = outcome_bar_panel, height = 2, width = 3.75)
# ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/fig2_panelsB.png"), plot = conf_loss_panel, height = 3, width = 3.75)
# ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/fig2_panelsA_B.png"), plot = f2_panA_B, height = 5, width = 3.75)
# 
# 
# ######################################
# ###### panel C clone size distribution
# ######################################
# datdir <- glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/")
# file <- "cloneSizes.csv"
# filepath <- paste(datdir,file, sep="")
# outcomes <- tibble(read.csv(filepath))
# TIME_POINT = 1
# outcomes$cloneSizes <- lapply(outcomes$cloneSizes, function(x) fromJSON(x))
# outcomes <- outcomes %>%
#   rowwise() %>%
#   mutate(aveCloneSize = mean(as.numeric(cloneSizes)))
# outcomes <- outcomes %>%
#   filter(year == TIME_POINT)
# max_value <- max(sapply(outcomes$cloneSizes, max)) + 1
# outcomes <- outcomes %>%
#   rowwise() %>%
#   mutate(cloneCounts = list(sapply(1:max_value, function(x) sum(as.numeric(cloneSizes) >=x))))
# outcomes <- outcomes %>%
#   unnest(cloneCounts)
# outcomes <- outcomes %>% 
#   group_by(blockProb) %>%
#   mutate(xVal = row_number() - 1) %>%
#   ungroup()
# muted <- c('#bbe4dd', '#91d3c8', '#67c2b3', '#44aa99', '#338073', '#22564d', '#122c28')
# names(muted) <- c("0", "0.001", "0.01", "0.1", "0.2", "0.5", "1")
# 
# #x_tick <- 400
# cloneSizeDist <- outcomes %>%
#   ggplot(aes(x = xVal, y = cloneCounts, color=as.factor(blockProb))) +
#   geom_line(linewidth=1) +
#   scale_color_manual("p", values = muted) +
#   theme_classic() +
# 
#   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^floor(x)),
#                 labels = trans_format("log10", math_format(10^.x))) +
#   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                 labels = trans_format("log10", math_format(10^.x))) +
#   
#   theme(legend.position = c(.22,.375), legend.title.align=0.5, legend.title = element_text(face = "italic")) +
#   guides(color = guide_legend(nrow=4)) +
#   
#   #scale_y_continuous(breaks = seq(0,100,20)) +
#   #scale_x_continuous(breaks = seq(0,max_value+x_tick*2,x_tick), limits = c(0, max_value + x_tick)) +
#   #theme(legend.position = c(.9,.55), legend.title.align=0.5, legend.title = element_text(face = "italic")) +
#   
#   
#   labs(x = "Clone size (x)",
#        y = "Clones >= x",
#        #title = paste("year:", TIME_POINT)
#        ) + 
#   theme(text=element_text(size=8), 
#                 legend.background = element_blank(),
#                 legend.key.size = unit(0.5, "lines"),
#                 legend.key.height = unit(0.5, "lines"),
#                 legend.margin = margin(0, 0, 0, 0),
#                 legend.key.spacing.y = unit(0.0001, "cm"),
#                 legend.key.spacing.x = unit(0.0001, "cm"),
#                 axis.title = element_text(size = 8),
#                 axis.text = element_text(size = 8, color = 'black'),
#                 legend.text=element_text(size=8))
# 
# show(cloneSizeDist)
# f1_panA_B_C <- plot_grid(conf_loss_panel, cloneSizeDist,outcome_bar_panel, ncol = 2, labels = c("A", "B", "C"), rel_heights = c(1, 0.65))
# f1_panA_B_C
# 
# ######################################
# ###### panel D area of clone spread
# ######################################
# areaSpread %>% filter(corrBlock == 0.1)
# areaSpread <- tibble(read.csv(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/averageAreaSpread.csv")))
# muted <- c('#bbe4dd', '#91d3c8', '#67c2b3', '#44aa99', '#338073', '#22564d', '#122c28')
# names(muted) <- c("0", "0.001", "0.01", "0.1", "0.2", "0.5", "1")
# print(areaSpread, n=Inf)
# ave_countOnly <- ggplot(areaSpread, aes(x = year, y = avg_tissueArea, color = as.factor(corrBlock))) +
#   geom_line(size = 1.6) +
#   scale_color_manual("p", values = muted) +
#   labs(x = "Year", y = expression(paste('Area (mm' ^2*')' )), color = "CorrBlock") +
#   theme_classic() +
#   theme(legend.position = c(0.9,0.4), legend.title.align=0.5, legend.title = element_text(face = "italic"))  + 
#   theme(legend.background=element_rect(fill = alpha("white", 0))) +
#   ylim(0,0.6) + 
#   theme(text=element_text(size=8), 
#         legend.background = element_blank(),
#         legend.key.size = unit(0.5, "lines"),
#         legend.key.height = unit(0.5, "lines"),
#         legend.margin = margin(0, 0, 0, 0),
#         legend.key.spacing.y = unit(0.0001, "cm"),
#         legend.key.spacing.x = unit(0.0001, "cm"),
#         axis.title = element_text(size = 8),
#         axis.text = element_text(size = 8, color = 'black'),
#         legend.text=element_text(size=8))
# show(ave_countOnly)
# 
# f1_panA_B_C_D <- plot_grid(outcome_bar_panel, cloneSizeDist, conf_loss_panel,  ave_countOnly, ncol = 2, labels = c("B", "D", "C", "E"), label_size = 12, rel_widths = c(1,1))
# show(f1_panA_B_C_D)
# #f1_panA_B_C_D <- plot_grid(outcome_bar_panel,  cloneSizeDist, conf_loss_panel,ave_countOnly, ncol = 2, rel_heights = c(1,1))
# 
# f2_panC_D <- plot_grid(cloneSizeDist, ave_countOnly, ncol = 2, rel_heights = c(1,1))
# 
# ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/fig2_panelsC_D.png"), plot = f2_panC_D, height = 2.5, width = 7.5)
# 
# show(f1_panA_B_C_D)
# #ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/fig2_panelsA_B_C_D.png"), plot = f1_panA_B_C_D, height = 5, width = 4)
# ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/fig2_panelsA_B_C_D.png"), plot = f1_panA_B_C_D, height = 2.75, width = 5)
