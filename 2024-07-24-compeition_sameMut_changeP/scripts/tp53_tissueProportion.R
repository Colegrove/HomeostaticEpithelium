## Hunter Colegrove
## FA dynamics
##
## Reads in position_data_all from compileVisFiles.R
## Plots number of tp53 mutations per replicate


require(foreach)
require(tidyverse)
require(glue)
library(gridExtra)
require(grid)
library(imager)
library(gtools)
library(png)
library(ggtext)


cluster = TRUE
#cluster = FALSE
#threshold_vector = c(10,50,100,400)
threshold_vector = c(108)
project_path = "2024-07-24-tp53_competition_sameMut_changeP"
set.seed(223)
## read in data from compileVisFiles.R
if(cluster){
  alldat <- read_csv(glue("/net/feder/vol1/project/HomeostaticEpithelium/dat/{project_path}/position_data_All.csv"))
}else{
  #alldat <- read_csv(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/position_data_All.csv"))
  alldat <- read_csv(glue("/Volumes/feder-vol1/project/HomeostaticEpithelium/dat/{project_path}/position_data_All.csv"))
}

######## 
## Average number of total tp53 mutations per replicate
## Includes tp53 + FA or tp53 + correction
## change to mut == 68 for only tp53 + FA
########

## filter for replicates that were not lost
## find n replicates for each tp53Block/year

completed_reps <- alldat %>%
  group_by(corrBlock, FAp53block, corrTime, replicate) %>%
  mutate(max_year = max(year)) %>%
  filter(max_year == 46) %>%
  group_by(corrBlock, year, FAp53block, corrTime) %>%
  summarise(completed_replicates = n_distinct(replicate)) %>%
  filter(!(year<1)) %>%
  filter(!(corrBlock == 0 & corrTime == 1)) %>%
  filter(!(FAp53block == 0.01)) %>%
  ungroup()

cat('completed_reps')
completed_reps %>% print(n=Inf)

total_reps = completed_reps %>% filter(corrBlock == 0 & year == 46 & FAp53block == 0.02) %>%
  pull(completed_replicates)

tissue_size = 4900 # 70x70 basal layer = 4900
total_size = tissue_size*total_reps

## counts the number of tp53 cells within each replicate/year
## change to mut == 68 for only tp53 + FA
tp53clones <- alldat %>%
  filter(!(corrBlock == 0 & corrTime == 1)) %>%
  filter(!(year < 1)) %>%
  filter(!(FAp53block == 0.01)) %>%
  filter(mut %in% c(68, 777)) %>% # tp53 + FA or tp53 + correction
  group_by(replicate, corrBlock, year, FAp53block) %>%
  summarize(count = n(), .groups = 'drop') %>%
  filter(count >= threshold_vector) %>% ## lower limit of detection
  mutate(proportion = count/tissue_size)

cat("tp53clones")
print(tp53clones)
write.csv(tp53clones, glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-07-24-competition_sameMut_changeP/dat/tp53clones_df.csv"))

### Fill in 0's for simulations that acquired no tp53 mutations 
all_years <- seq(1, 46, by = 5)
tp53clones_noCorrection <- tp53clones %>% filter(corrBlock == 0) %>%
  group_by(year, FAp53block) %>%
  mutate(replicate = row_number()-1) %>%
  ungroup

expanded_tp53clones <- expand_grid(replicate = seq(0,total_reps-1), 
                                   corrBlock = 0, 
                                   year = all_years, 
                                   FAp53block = unique(tp53clones_noCorrection$FAp53block))

# Merge with the original data and fill missing values with 0
tp53clones_full <- expanded_tp53clones %>%
  left_join(tp53clones_noCorrection, by = c("replicate", "corrBlock", "year", "FAp53block")) %>%
  replace_na(list(count = 0, proportion = 0))

cat("tp53clones_full")
print(tp53clones_full)

### randomly sample lost simulations from the no correction case
missing_reps <- completed_reps %>% filter(corrBlock != 0) %>%
  mutate(missing_replicates = total_reps - completed_replicates) %>%
  mutate(missing_replicates = if_else(missing_replicates < 1, 0, missing_replicates))


results <- list()
for (i in 1:nrow(missing_reps)){
  filter_corrBlock <- missing_reps$corrBlock[i]
  filter_year <- missing_reps$year[i]
  missing_replicates <- missing_reps$missing_replicates[i]
  filter_FAp53block <- missing_reps$FAp53block[i]
  
  sampled_count <- tp53clones_full %>%
    filter(year == filter_year) %>%
    filter(FAp53block == filter_FAp53block) %>%
    filter(corrBlock == 0) %>%
    pull(count) %>%
    sample(missing_replicates, replace=TRUE)
  
  #print(sampled_count)
  # Append results to the list
  results[[i]] <- data.frame(
    corrBlock = rep(filter_corrBlock, missing_replicates),
    year = rep(filter_year, missing_replicates),
    FAp53block = rep(filter_FAp53block, missing_replicates),
    count = sampled_count
  )
}

sampled_dataframe <- bind_rows(results)
cat("sampled_dataframe")
print(sampled_dataframe)
tp53clones_sampled <- bind_rows(tp53clones, sampled_dataframe)

cat("tp53clones_sampled")
print(tp53clones_sampled)

tp53clones_sampled_summarized <- tp53clones_sampled %>%
  group_by(corrBlock, year, FAp53block) %>%
  summarize(total_count = sum(count), .groups = 'drop') %>%
  mutate(total_proportion = total_count/total_size)
cat("tp53clones_sampled_summarized")
print(tp53clones_sampled_summarized)

## if no mutations at a given timepoint, set proportion to 0

combinations <- expand.grid(
  corrBlock = unique(tp53clones_sampled_summarized$corrBlock),
  FAp53block = unique(tp53clones_sampled_summarized$FAp53block),
  year = all_years
)


#### for full non-summarized data
combinations_full <- expand.grid(
  corrBlock = unique(tp53clones_sampled_summarized$corrBlock),
  FAp53block = unique(tp53clones_sampled_summarized$FAp53block),
  year = all_years, 
  replicate = seq(0,total_reps-1)
)

full_reps_summary <- combinations_full %>%
  left_join(tp53clones_sampled, by = c("corrBlock", "year", "FAp53block", "replicate"))

mean_proportions <- combinations %>%
  left_join(tp53clones_sampled_summarized, by = c("corrBlock", "year", "FAp53block")) %>%
  mutate(
    total_proportion = ifelse(is.na(total_proportion), 0, total_proportion),
    total_count = ifelse(is.na(total_count), 0, total_count),
  )
print(mean_proportions)

## plot the mean proportion
#muted <- c('#bbe4dd', '#91d3c8', '#67c2b3', '#44aa99', '#338073', '#22564d', '#122c28')
#names(muted) <- c("0", "0.001", "0.01", "0.1", "0.2", "0.5", "1")

#custom_labels <- c("0.375" = "10.6X", "0.5" = "8X", "1" = "4X", "2" = "2X", "2.666" = "1.5X")
#blue_shades <- c("0.375" = "#012A4A", "0.5" = "#01497c", "1" = "#2a6f97", "2" = "#468faf", "2.666" = "#89c2d9")


custom_labels <- c( "0.04" = "4X", "0.02" = "2X", "0.015" = "1.5X")
red_shades <- c("0.04" = "#882255", "0.02" = "#CC6677", "0.015" = "#EE99AA")

custom_breaks <- seq(0,46,by=10)


## for line labels
dat_label <- mean_proportions %>%
  group_by(corrBlock, FAp53block) %>%
  filter(year == max(year)) %>%
  ungroup()
#dat_label

df_labs <- data.frame(
  corrBlock = c(0, 0.01, 0.1),  # The unique values of corrBlock
  corrBlock_label = c("No~correction", "italic(p)[italic(FANC)^'+']^italic(TP53)^'+' == 0.01", "italic(p)[italic(FANC)^'+']^italic(TP53)^'+' == 0.1")  # Math expressions
)
mean_proportions <- merge(mean_proportions, df_labs, by = "corrBlock", all.x = TRUE)
mean_proportions$corrBlock_label <- factor(mean_proportions$corrBlock_label,
                                           levels = c("No~correction", "italic(p)[italic(FANC)^'+']^italic(TP53)^'+' == 0.01", "italic(p)[italic(FANC)^'+']^italic(TP53)^'+' == 0.1"))
mean_proportions$FAp53block <- factor(mean_proportions$FAp53block, 
                                      levels = c("0.04", "0.02", "0.015"))  # Order by 4X, 2X, 1.5X

write.csv(mean_proportions, glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-07-24-competition_sameMut_changeP/dat/mean_proportion_plot_data.csv"))

write.csv(tp53clones_sampled, glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-07-24-competition_sameMut_changeP/dat/raw_proportion_plot_data.csv"))



#mean_proportions <- read_csv(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-07-24-competition_sameMut_changeP/dat/mean_proportion_plot_data.csv"))

mean_proportions$corrBlock_label <- factor(mean_proportions$corrBlock_label,
                                          levels = c("No~correction", "italic(p)[italic(FANC)^'+']^italic(TP53)^'+' == 0.01", "italic(p)[italic(FANC)^'+']^italic(TP53)^'+' == 0.1"))
mean_proportions$FAp53block <- factor(mean_proportions$FAp53block,
                                     levels = c("0.04", "0.02", "0.015"))  # Order by 4X, 2X, 1.5X
mean_proportions

plot <- ggplot(mean_proportions, aes(x = year, y = total_proportion, color = as.factor(FAp53block), group = FAp53block)) +
  geom_line(size = 1) +
  facet_wrap(~corrBlock_label, scales = "fixed", nrow = 1, labeller = label_parsed) +
  labs(x = "Year", y = "Tissue proportion", color = bquote(italic(p) [italic(FANC)^"-"] ^italic(TP53)^"-" ~ "increase")) +
  scale_color_manual(values = red_shades, labels = custom_labels) +
  theme_minimal() +
  #scale_x_continuous(limits = c(1, 46), breaks = custom_breaks) +
  #scale_x_continuous(breaks = seq(1, 50, 10)) +
  scale_x_continuous(breaks = c(1, 10, 20, 30, 40, 50)) +
  scale_y_continuous(limits = c(0, 0.4), breaks = seq(0,0.4, by = 0.1)) +
  #theme(legend.position = c(0.2, 0.55),
  theme(legend.position = c(0.1, 0.6))+
        #legend.title = element_text(face = "italic")) +
  coord_cartesian(clip = "off") +
  theme(
    #legend.text=element_text(size=8),
    legend.key.size = unit(0.5, "lines"),
    legend.key.height = unit(0.5, "lines"),
    legend.margin = margin(0, 0, 0, 0),
    legend.key.spacing.y = unit(0.0001, "cm"),
    legend.key.spacing.x = unit(0.0001, "cm"),
    
    strip.text = element_text(size = 8, color='black'),
    axis.text = element_text(size = 8, color='black'),
    axis.title = element_text(size = 8, color='black'), 
    legend.text = element_text(size = 8, color='black'), 
    legend.title = element_text(size = 8, color='black'),
    axis.line = element_line(color="black"),
    axis.ticks = element_line(color = "black"),
    plot.margin = margin(0.5, 25.5, 0.5, 5.5), 
    panel.spacing = unit(1.5, "lines")
  )
show(plot)

plot_noLabels <- plot + 
  theme(strip.text = element_blank(), 
        axis.title.y = element_text(hjust = 0.75),
        plot.margin = margin(10, 25.5, 0.5, 5.5), 
        legend.title = element_text(margin = margin(b = 0))
  )

show(plot_noLabels)
if(cluster){
  ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-07-24-competition_sameMut_changeP/results/tp53_tissueProportion.png"), plot, width = 6, height = 1.35, units = "in", dpi = 300)
}else{
  ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-07-24-competition_sameMut_changeP/results/tp53_tissueProportion.png"), plot, width = 5.8, height = 1.35, units = "in", dpi = 300)
  ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-07-24-competition_sameMut_changeP/results/tp53_tissueProportion_noFacetLabels.png"), plot_noLabels, width = 5.8, height = 1, units = "in", dpi = 300)
  
  }

q()

#########################################################
############
#########################################################

# ## filter for replicates that were not lost
# ## find n replicates for each tp53Block/year
# alldat %>% print()
# completed_reps <- alldat %>%
#   group_by(corrBlock, FAp53block, corrTime, replicate) %>%
#   mutate(max_year = max(year)) %>%
#   filter(max_year == 46) %>%
#   group_by(corrBlock, year, FAp53block, corrTime) %>%
#   summarise(completed_replicates = n_distinct(replicate)) %>%
#   filter(!(year<1)) %>%
#   filter(!(corrBlock == 0 & corrTime == 1)) %>%
#   ungroup()
# completed_reps %>% print(n=Inf)
# 
# total_reps = completed_reps %>% filter(corrBlock == 0 & year == 46 & FAp53block == 0.1) %>%
#   pull(completed_replicates)
# 
# 
# 
# tissue_size = 4900 # 70x70 basal layer = 4900
# total_size = tissue_size*total_reps
# cat("total_reps")
# cat(total_reps)
# ## counts the number of tp53 cells within each replicate/year
# ## change to mut == 68 for only tp53 + FA
# tp53clones <- alldat %>%
#   filter(!(corrBlock == 0 & corrTime == 1)) %>% ## plots no correction case with corrBlock == 0
#   filter(!(year < 1)) %>%
#   filter(mut %in% c(68, 777)) %>% # tp53 + FA or tp53 + correction
#   group_by(replicate, corrBlock, year, FAp53block) %>%
#   summarize(count = n(), .groups = 'drop') %>%
#   filter(count >= threshold_vector) %>% ## lower limit of detection
#   mutate(proportion = count/tissue_size)
# tp53clones
# ### Fill in 0's for simulations that acquired no tp53 mutations 
# all_years <- seq(1, 46, by = 5)
# tp53clones_noCorrection <- tp53clones %>% filter(corrBlock == 0)
# 
# cat("total_reps-1")
# cat(total_reps-1)
# expanded_tp53clones <- expand_grid(replicate = seq(0,total_reps-1), 
#                                    corrBlock = 0, 
#                                    year = all_years, 
#                                    FAp53block = unique(tp53clones_noCorrection$FAp53block))
# 
# # Merge with the original data and fill missing values with 0
# tp53clones_full <- expanded_tp53clones %>%
#   left_join(tp53clones_noCorrection, by = c("replicate", "corrBlock", "year", "FAp53block")) %>%
#   replace_na(list(count = 0, proportion = 0))
# 
# ### randomly sample lost simulations from the no correction case
# missing_reps <- completed_reps %>% filter(corrBlock != 0) %>%
#   mutate(missing_replicates = total_reps - completed_replicates)
# 
# 
# results <- list()
# for (i in 1:nrow(missing_reps)){
#   filter_corrBlock <- missing_reps$corrBlock[i]
#   filter_year <- missing_reps$year[i]
#   missing_replicates <- missing_reps$missing_replicates[i]
#   filter_FAp53block <- missing_reps$FAp53block[i]
#   
#   sampled_count <- tp53clones_full %>%
#     filter(year == filter_year) %>%
#     filter(FAp53block == filter_FAp53block) %>%
#     pull(count) %>%
#     sample(missing_replicates, replace=TRUE)
#   
#   # Append results to the list
#   results[[i]] <- data.frame(
#     corrBlock = rep(filter_corrBlock, missing_replicates),
#     year = rep(filter_year, missing_replicates),
#     FAp53block = rep(filter_FAp53block, missing_replicates),
#     count = sampled_count
#   )
# }
# 
# sampled_dataframe <- bind_rows(results)
# tp53clones_sampled <- bind_rows(tp53clones, sampled_dataframe)
# 
# 
# tp53clones_sampled_summarized <- tp53clones_sampled %>%
#   group_by(corrBlock, year, FAp53block) %>%
#   summarize(total_count = sum(count), .groups = 'drop') %>%
#   mutate(total_proportion = total_count/total_size)
# tp53clones_sampled_summarized
# 
# ## if no mutations at a given timepoint, set proportion to 0
# 
# combinations <- expand.grid(
#   corrBlock = unique(tp53clones_sampled_summarized$corrBlock),
#   FAp53block = unique(tp53clones_sampled_summarized$FAp53block),
#   year = all_years
# )
# 
# mean_proportions <- combinations %>%
#   left_join(tp53clones_sampled_summarized, by = c("corrBlock", "year", "FAp53block")) %>%
#   mutate(
#     total_proportion = ifelse(is.na(total_proportion), 0, total_proportion),
#     total_count = ifelse(is.na(total_count), 0, total_count),
#   )
# mean_proportions
# 
# ## plot the mean proportion
# muted <- c('#bbe4dd', '#91d3c8', '#67c2b3', '#44aa99', '#338073', '#22564d', '#122c28')
# names(muted) <- c("0", "0.001", "0.01", "0.1", "0.2", "0.5", "1")
# 
# custom_labels <- c("0.375" = "10.6X", "0.5" = "8X", "1" = "4X", "2" = "2X", "2.666" = "1.5X")
# blue_shades <- c("0.375" = "#012A4A", "0.5" = "#01497c", "1" = "#2a6f97", "2" = "#468faf", "2.666" = "#89c2d9")
# 
# custom_breaks <- seq(0,46,by=10)
# corrBlock_labels <- c(
#   "0" = "No correction",
#   "0.01" = "0.01",
#   "0.1" = "0.1"
# )
# 
# ## for line labels
# dat_label <- mean_proportions %>%
#   group_by(corrBlock, FAp53block) %>%
#   filter(year == max(year)) %>%
#   ungroup()
# dat_label
# 
# plot <- ggplot(mean_proportions, aes(x = year, y = total_proportion, color = as.factor(FAp53block), group = FAp53block)) +
#   geom_line(size = 2) +
#   facet_wrap(~corrBlock, scales = "fixed", nrow = 1, labeller = as_labeller(corrBlock_labels)) +
#   labs(x = "Year", y = "Tissue proportion", color = "Persistence coefficient") +
#   scale_color_manual(values = blue_shades, labels = custom_labels) +
#   theme_minimal() +
#   #scale_x_continuous(limits = c(1, 46), breaks = custom_breaks) +
#   scale_x_continuous(breaks = seq(1, 50, 5)) +
#   scale_y_continuous(limits = c(0, 0.25), breaks = seq(0,0.25, by = 0.05)) +
#   #theme(legend.position = c(0.2, 0.55),
#   theme(legend.position = "none",
#         legend.title = element_text(face = "italic")) +
#   geom_richtext(data = dat_label, aes(label = custom_labels[as.character(FAp53block)]), 
#                 color = "black", label.size = NA, fill = NA,
#                 label.margin = unit(4, "pt"),
#                 label.padding = unit(3, "pt"),
#                 hjust = 0, show.legend = FALSE) +
#   theme(
#     strip.text = element_text(size = 8, color='black'),
#     axis.text = element_text(size = 8, color='black'),
#     axis.title = element_text(size = 8, color='black'), 
#     legend.text = element_text(size = 8, color='black'), 
#     legend.title = element_text(size = 8, color='black'),
#     axis.line = element_line(color="black"),
#     axis.ticks = element_line(color = "black")
#   )
# show(plot)
# 
# if(cluster){
#   ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/tp53_tissueProportion.png"), plot, width = 4.5, height = 2, units = "in", dpi = 300)
# }else{
#   ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/tp53_tissueProportion.png"), plot, width = 4.5, height = 2, units = "in", dpi = 300)
# }
# 
# q()
# 
# 
# 
# 
# 
# alldat
# 
# tissue_size = 4900 # 70x70 basal layer = 4900
# 
# ## counts the number of cells that make up each tp53 clone within each replicate
# ## change to mut == 68 for only tp53 + FA
# tp53clones <- alldat %>%
#   filter(!(corrBlock == 0 & corrTime == 1)) %>% ## plots no correction case with corrBlock == 0
#   filter(!(year < 1)) %>%
#   filter(mut %in% c(68, 777)) %>% # tp53 + FA or tp53 + correction
#   group_by(replicate, corrBlock, year, FAp53block) %>%
#   summarize(count = n(), .groups = 'drop') %>%
#   filter(count >= threshold_vector) %>% ## lower limit of detection
#   mutate(proportion = count/tissue_size)
# tp53clones
# 
# ## plot all trajectories
# ggplot(tp53clones, aes(x = year, y = proportion, color = as.factor(replicate), group = replicate)) +
#   geom_line() +
#   labs(x = "Year", y = "Tissue proportion", color = "Replicate") +
#   facet_wrap(~ corrBlock) +
#   theme_minimal() +
#   theme(
#     strip.text = element_text(size = 12),
#     axis.text = element_text(size = 10),
#     axis.title = element_text(size = 12), 
#     legend.position = "none"
#   )
# 
# ## calculate the mean and IQR of tp53 muts
# mean_proportions <- tp53clones %>%
#   group_by(corrBlock, year, FAp53block) %>%
#   summarize(
#     mean_proportion = mean(proportion), 
#     Q1 = quantile(proportion, 0.25),  # Lower quartile
#     Q3 = quantile(proportion, 0.75),  # Upper quartile
#     .groups = 'drop'
#   )
# mean_proportions
# ## if no mutations at a given timepoint, set proportion to 0
# all_years <- seq(1, 46, by = 5)
# 
# combinations <- expand.grid(
#   corrBlock = unique(tp53clones$corrBlock),
#   FAp53block = unique(tp53clones$FAp53block),
#   year = all_years
# )
# 
# mean_proportions <- combinations %>%
#   left_join(mean_proportions, by = c("corrBlock", "year", "FAp53block")) %>%
#   mutate(
#     mean_proportion = ifelse(is.na(mean_proportion), 0, mean_proportion),
#     Q1 = ifelse(is.na(Q1), 0, Q1),
#     Q3 = ifelse(is.na(Q3), 0, Q3)
#   )
# 
# mean_proportions
# ## plot the mean proportion
# muted <- c('#bbe4dd', '#91d3c8', '#67c2b3', '#44aa99', '#338073', '#22564d', '#122c28')
# names(muted) <- c("0", "0.001", "0.01", "0.1", "0.2", "0.5", "1")
# custom_breaks <- seq(0,46,by=10)
# 
# plot <- ggplot(mean_proportions, aes(x = year, y = mean_proportion, color = as.factor(corrBlock), group = corrBlock)) +
#   geom_line() +
#   #geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.05, position = position_dodge(width = 0.7)) +  # Add error bars for IQR
#   labs(x = "Year", y = "Tissue proportion", color = "Persistence coefficient") +
#   facet_wrap(~ FAp53block) +
#   scale_color_manual("p", values = muted) +
#   theme_minimal() +
#   scale_x_continuous(limits = c(1, 46), breaks = custom_breaks) +
#   scale_y_continuous(limits = c(0, 1), breaks = seq(0,1, by = 0.2)) +
#   theme(legend.position = c(0.2, 0.55),
#         legend.title = element_text(face = "italic")) +
#   theme(
#     strip.text = element_text(size = 8),
#     axis.text = element_text(size = 8),
#     axis.title = element_text(size = 8), 
#     legend.text = element_text(size = 8), 
#     legend.title = element_text(size = 8)
#   )
# show(plot)
# 
# if(cluster){
#   ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-07-24-competition_sameMut_changeP/results/tp53_tissueProportion.png"), plot, width = 2, height = 2, units = "in", dpi = 300)
# }else{
#   ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-07-24-competition_sameMut_changeP/results/tp53_tissueProportion.png"), plot, width = 2, height = 2, units = "in", dpi = 300)
# }
# 
# q()
# 
# 
# 
# ## filter for replicates that were not lost
# ## find n replicates for each tp53Block/year
# completed_reps <- alldat %>%
#   group_by(corrBlock, FAp53block, corrTime, replicate) %>%
#   mutate(max_year = max(year)) %>%
#   filter(max_year == 46) %>%
#   group_by(corrBlock, year, FAp53block, corrTime) %>%
#   summarise(completed_replicates = n_distinct(replicate)) %>%
#   filter(!(year<1)) %>%
#   filter(!(corrBlock == 0 & corrTime == 1)) %>%
#   ungroup()
# completed_reps %>% print(n=Inf)
# 
# ## counts the number of cells that make up each tp53 clone within each replicate
# ## change to mut == 68 for only tp53 + FA
# tp53clones <- alldat %>%
#   filter(!(corrBlock == 0 & corrTime == 1)) %>% ## plots no correction case with corrBlock == 0
#   filter(!(year < 1)) %>%
#   dplyr::select(tp53Clone, replicate, corrBlock, year, mut, FAp53block) %>%
#   filter(mut %in% c(68, 777)) %>% # tp53 + FA or tp53 + correction
#   #filter(mut ==68) %>% # tp53 + FA only (no correction)
#   group_by(tp53Clone, replicate, corrBlock, year, mut, FAp53block) %>%
#   summarise(cells_count = n())
# tp53clones
# 
# ## filter data based on detection thresholds and combine all dataframes
# df_list <- list()
# for (val in threshold_vector) {
#   duplicated_df <- tp53clones %>%
#     filter(cells_count >= val) %>%
#     mutate(threshold = val)
#   df_list <- append(df_list, list(duplicated_df))
# }
# tp53_thresholds <- do.call(rbind, df_list)
# 
# ## counts the unique number of tp53 mutations within each replicate at each time
# unique_tp53Clone <- tp53_thresholds %>%
#   group_by(corrBlock, year, replicate, threshold, FAp53block) %>%
#   summarise(unique_count = n_distinct(tp53Clone)) %>%
#   ungroup() %>%
#   group_by(corrBlock, year, threshold, FAp53block) %>%
#   summarise(total_unique_count = sum(unique_count)) %>%
#   ungroup()
# unique_tp53Clone %>% print(n=Inf)
# 
# ## remove columns from completed dataframe
# # completed_reps <- completed_reps %>%
# #   dplyr::select(-mutRate, -corrTime)
# # completed_reps %>% print(n=Inf)
# 
# ## combine tp53 mut counts with the number of replicates completed in experiment
# unique_tp53Clone_reps <- completed_reps %>%
#   left_join(unique_tp53Clone, by = c("corrBlock", "year", "FAp53block"))
# unique_tp53Clone_reps
# ## average mutations per replicate
# ave_unique_count <- unique_tp53Clone_reps %>%
#   mutate(ave_unique_count = total_unique_count/completed_replicates) %>%
#   drop_na()
# print(ave_unique_count)
# 
# ## label mutation modifier values to read as an "X increase"
# # custom_labels <- c("1" = "4X", "2" = "2X", "2.666" = "1.5X")
# # blue_shades <- c("1" = "#03045e", "2" = "#0077b6", "2.666" = "#48cae4")
# # 
# # ## extra shades
# # custom_labels <- c("0.375" = "10.6X", "0.5" = "8X", "1" = "4X", "2" = "2X", "2.666" = "1.5X")
# # blue_shades <- c("0.375" = "#012A4A", "0.5" = "#01497c", "1" = "#2a6f97", "2" = "#468faf", "2.666" = "#89c2d9")
# 
# custom_labels <- c("0.01" = "0.01", "0.015" = "0.015", "0.02" = "0.02")
# blue_shades <- c("0.01" = "#03045e", "0.015" = "#0077b6", "0.02" = "#48cae4")
# 
# ## construct a label df of the final timepoint for line labels
# dat_label <- ave_unique_count %>%
#   group_by(corrBlock, FAp53block) %>%
#   filter(year == max(year)) %>%
#   ungroup()
# dat_label
# 
# ## facet labels
# corrBlock_labels <- c("0" = "No correction", "0.01" = "0.01", "0.1" = "0.1")
# 
# ## Plot using ggplot
# plot <- ggplot(ave_unique_count, aes(x = year, y = ave_unique_count, color = factor(FAp53block))) +
#   geom_line(size = 2) +
#   labs(x = "Year", y = "TP53 mutation frequency", color = "FA p53 persistence") +
#   theme_minimal() + 
#   #theme(aspect.ratio = 0.75) + 
#   scale_x_continuous(breaks = seq(1, 50, 5)) +
#   facet_wrap(~corrBlock, scales = "fixed", nrow = 1, labeller = as_labeller(corrBlock_labels)) +
#   theme(text = element_text(size = 20), panel.spacing = unit(2,'lines')) +
#   scale_color_manual(values = blue_shades, labels = custom_labels) +
#   #theme(legend.position = c(0.8, 0.8))
#   theme(legend.position = "none") +
#   geom_richtext(data = dat_label, aes(label = custom_labels[as.character(FAp53block)]), 
#                 color = "black", label.size = NA, fill = NA,
#                 label.margin = unit(4, "pt"),
#                 label.padding = unit(3, "pt"),
#                 hjust = 0, show.legend = FALSE) +
#   coord_cartesian(clip = "off") +
#   theme(plot.margin = margin(5.5, 35, 5.5, 5.5))
#   
# show(plot)
# ## save plot to file
# if(cluster){
#   #ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/tp53_mutations_per_simulation_thresholds.png"), plot, width = 15, height = 6, units = "in", dpi = 300)
#   ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-07-24-competition_sameMut_changeP/results/tp53_mutations_per_simulation_thresholds.png"), plot, width = 15, height = 6, units = "in", dpi = 300)
#   
# }else{
#   ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/tp53_mutations_per_simulation_thresholds.png"), plot, width = 15, height = 6, units = "in", dpi = 300)
# }
# 
# 
# ############
# ##### Fill missing values
# ############
# 
# ### Determine the missing replicates to be added.
# summary_data <- alldat %>%
#   filter(!(year < 1)) %>%
#   filter(!(corrBlock == 0 & corrTime == 1)) %>%
#   group_by(corrBlock, year, replicate, mutRate) %>%
#   summarize(observations = n())
# summary_data %>% print()
# print(unique(summary_data$mutRate))
# 
# 
# max_replicates = n_distinct(unique(summary_data$replicate[summary_data$corrBlock ==0]))
# 
# #summary_data$replicate[summary_data$corrBlock == 0]
# 
# missing_replicates <- lapply(years_total, function(year_value) {
#   summary_data %>%
#     filter(year == year_value) %>%
#     group_by(corrBlock, mutRate) %>%
#     print() %>%
#     summarize(available_replicates = n_distinct(replicate)) %>%
#     mutate(missing_replicates = max_replicates - available_replicates) %>%
#     mutate(missing_replicates = ifelse(missing_replicates < 0, 0, missing_replicates)) %>%
#     mutate(year = year_value)
# })
# 
# missing_replicates <- do.call(rbind, missing_replicates)
# print(missing_replicates, n=Inf)
# missing_replicates <- missing_replicates %>% filter(missing_replicates >0)
# 
# tp53clones <- alldat %>%
#   filter(!(corrBlock == 0 & corrTime == 1)) %>% ## plots no correction case with corrBlock == 0
#   filter(!(year < 1)) %>%
#   dplyr::select(tp53Clone, replicate, corrBlock, year, mut, mutRate) %>%
#   filter(mut %in% c(68, 777)) %>% # tp53 + FA or tp53 + correction
#   #filter(mut ==68) %>% # tp53 + FA only (no correction)
#   group_by(tp53Clone, replicate, corrBlock, year, mut, mutRate) %>%
#   summarise(cells_count = n())
# tp53clones
# tp53clones %>% print(n=Inf)
# tp53clones
# 
# missing_replicates
# for (i in 1:nrow(missing_replicates)){
#   n_missing <- missing_replicates$missing_replicates[i]
#   corrBlock_replace <- missing_replicates$corrBlock[i]
#   year <- missing_replicates$year[i]
#   mutRate <- missing_replicates$mutRate[i]
#   for (j in 1:n_missing){
#     n_reps = 100
#     random_rep = sample(0:n_reps-1, size = 1)
#     sampled_rep <- tp53clones %>%
#       filter(corrBlock==0, year == year, replicate == random_rep, mutRate == mutRate)
#     #print(sampled_rep)
#     
#     sampled_rep <- mutate(sampled_rep, corrBlock = corrBlock_replace)
#     #print(corrBlock_replace)
#     #print(sampled_rep)
#     tp53clones <- bind_rows(tp53clones, sampled_rep)
#     
#   }
# }
# tp53clones
# 
# unique_reps <- tp53clones %>%
#   group_by(corrBlock, year, mutRate) %>%
#   summarize(unique_reps = n_distinct(replicate))
# unique_reps
# 
# 
# ## filter for replicates that were not lost
# ## find n replicates for each tp53Block/year
# completed_reps <- alldat %>%
#   group_by(corrBlock, mutRate, corrTime, replicate) %>%
#   mutate(max_year = max(year)) %>%
#   filter(max_year == 46) %>%
#   group_by(corrBlock, year, mutRate, corrTime) %>%
#   summarise(completed_replicates = n_distinct(replicate)) %>%
#   filter(!(year<1)) %>%
#   filter(!(corrBlock == 0 & corrTime == 1)) %>%
#   ungroup()
# completed_reps %>% print(n=Inf)
# ### APPLY ONLY IF RANDOM SAMPLING
# completed_reps <- mutate(completed_reps, completed_replicates = 100)
# 
# 
# ## counts the number of cells that make up each tp53 clone within each replicate
# ## change to mut == 68 for only tp53 + FA
# # tp53clones <- alldat %>%
# #   filter(!(corrBlock == 0 & corrTime == 1)) %>% ## plots no correction case with corrBlock == 0
# #   filter(!(year < 1)) %>%
# #   select(tp53Clone, replicate, corrBlock, year, mut) %>%
# #   filter(mut %in% c(68, 777)) %>% # tp53 + FA or tp53 + correction
# #   #filter(mut ==68) %>% # tp53 + FA only (no correction)
# #   group_by(tp53Clone, replicate, corrBlock, year, mut) %>%
# #   summarise(cells_count = n())
# # tp53clones %>% print(n=500)
# 
# 
# ## filter data based on detection thresholds and combine all dataframes
# df_list <- list()
# for (val in threshold_vector) {
#   duplicated_df <- tp53clones %>%
#     filter(cells_count >= val) %>%
#     mutate(threshold = val)
#   df_list <- append(df_list, list(duplicated_df))
# }
# tp53_thresholds <- do.call(rbind, df_list)
# 
# ## counts the unique number of tp53 mutations within each replicate at each time
# unique_tp53Clone <- tp53_thresholds %>%
#   group_by(corrBlock, year, replicate, threshold, mutRate) %>%
#   summarise(unique_count = n_distinct(tp53Clone)) %>%
#   ungroup() %>%
#   group_by(corrBlock, year, threshold, mutRate) %>%
#   summarise(total_unique_count = sum(unique_count)) %>%
#   ungroup()
# unique_tp53Clone %>% print(n=Inf)
# 
# ## remove columns from completed dataframe
# # completed_reps <- completed_reps %>%
# #   dplyr::select(-mutRate, -corrTime)
# # completed_reps %>% print(n=Inf)
# 
# ## combine tp53 mut counts with the number of replicates completed in experiment
# unique_tp53Clone_reps <- completed_reps %>%
#   left_join(unique_tp53Clone, by = c("corrBlock", "year", "mutRate"))
# 
# ## average mutations per replicate
# ave_unique_count <- unique_tp53Clone_reps %>%
#   mutate(ave_unique_count = total_unique_count/completed_replicates) %>%
#   drop_na()
# print(ave_unique_count)
# ave_unique_count <- ave_unique_count %>%
#   filter(threshold == 108)
# ave_unique_count
# 
# #### Add 0's to data to show up on plot even if no mutations
# #### must manually check data first if 0s can be added.
# # start_years <- c(1, 6, 11)
# # unique_corrBlock <- unique(ave_unique_count$corrBlock)
# # unique_threshold <- unique(ave_unique_count$threshold)
# # 
# # # Create a dataframe with the additional rows
# # additional_rows <- expand.grid(year = start_years, 
# #                                ave_unique_count = 0, 
# #                                corrBlock = unique_corrBlock, 
# #                                threshold = unique_threshold)
# # 
# # # Combine with the original dataframe
# # ave_unique_count_extended <- bind_rows(ave_unique_count, additional_rows) %>%
# #   arrange(year, corrBlock, threshold)
# # 
# # print(ave_unique_count_extended, n=Inf)
# 
# 
# ## Plot using ggplot
# num_colors <- length(unique(ave_unique_count$corrBlock))
# my_palette <- rev(viridis(num_colors))
# my_palette[1] <- "#BADE28FF"
# plot <- ggplot(ave_unique_count, aes(x = year, y = ave_unique_count, color = factor(mutRate))) +
#   geom_line(size = 3) +
#   ylim(0, 0.45) +
#   #labs(x = "Year", y = "TP53 mutation frequency", color = "FA mutant mutation rate") +
#   labs(x = "Year", y = expression(paste(italic("TP53"), " mutation frequency")), color = "p") +
#   theme_minimal() + 
#   scale_color_manual(values = my_palette, 
#                      breaks = unique(ave_unique_count$corrBlock), 
#                      labels = unique(ave_unique_count$corrBlock)) +
#   #theme(aspect.ratio = 0.75) + 
#   scale_x_continuous(breaks = seq(1, 50, 5)) +
#   #facet_wrap(~corrBlock, scales = "fixed", nrow = 1, labeller = as_labeller(corrBlock_labels)) +
#   facet_wrap(~corrBlock, scales = "fixed", nrow = 1) +
#   scale_color_manual(values = blue_shades, labels = custom_labels) +
#   theme(legend.position = c(0.8, 0.8)) +
#   theme(
#     axis.line = element_line(color="black"), 
#     axis.text = element_text(color="black", size=28),
#     axis.title.x = element_text(color="black", size=28),
#     axis.title.y = element_text(color="black", size=28),
#     axis.ticks = element_line(color = "black"),
#     legend.text = element_text(color="black", size=28), 
#     legend.title = element_text(color="black", face="italic", size=28),
#     #strip.background = element_blank(),
#     #strip.text.x = element_blank(),
#   )
# 
# 
# show(plot)
# 
# ## save plot to file
# if(cluster){
#   ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/tp53_mutations_per_simulation_thresholds.png"), plot, width = 20, height = 10, units = "in", dpi = 300)
# }else{
#   ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/tp53_mutations_per_simulation_thresholds_poster.png"), plot, width = 10, height = 6, units = "in", dpi = 300)
# }
# 
# 
