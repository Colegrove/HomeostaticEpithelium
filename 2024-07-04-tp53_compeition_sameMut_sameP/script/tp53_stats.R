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
require(viridis)
require(ggtext)

cluster = TRUE
set.seed(223)
#threshold_vector = c(10,50,100,400)
threshold_vector = c(108)
project_path = "2024-07-04-tp53_competition_sameMut_sameP"

## read in data from compileVisFiles.R
if(cluster){
  alldat <- read_csv(glue("/net/feder/vol1/project/HomeostaticEpithelium/dat/{project_path}/position_data_All.csv"))
}else{
  #alldat <- read_csv(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/position_data_All.csv"))
  alldat <- read_csv(glue("/Volumes/feder-vol1/project/HomeostaticEpithelium/dat/{project_path}/position_data_All.csv"))
}
alldat
#alldat1 <- read_csv(glue("/Volumes/feder-vol1/project/HomeostaticEpithelium/dat/2024-06-13-tp53_competition_sameMut_sameP/position_data_All.csv"))
#alldat2 <- read_csv(glue("/Volumes/feder-vol1/project/HomeostaticEpithelium/dat/2024-07-03-tp53_competition_sameMut_sameP/position_data_All.csv"))
alldat1 <- read_csv(glue("/net/feder/vol1/project/HomeostaticEpithelium/dat/2024-06-13-tp53_competition_sameMut_sameP/position_data_All.csv"))
alldat1
alldat2 <- read_csv(glue("/net/feder/vol1/project/HomeostaticEpithelium/dat/2024-07-03-tp53_competition_sameMut_sameP/position_data_All.csv"))
alldat2
alldat1 <- alldat1 %>% mutate(replicate = replicate + 100)
alldat2 <- alldat2 %>% mutate(replicate = replicate + 200)
alldat <- bind_rows(alldat, alldat1, alldat2)

rm(alldat1)
rm(alldat2)
###############
#### Proportion of basal layer occupied by tp53 mutations
#### Incluces tp53 + FA or tp53 + correction
#### change to mut == 68 for only tp53 + FA
###############


## filter for replicates that were not lost
## find n replicates for each tp53Block/year
completed_reps <- alldat %>%
  group_by(corrBlock, mutRate, corrTime, replicate) %>%
  mutate(max_year = max(year)) %>%
  filter(max_year == 46) %>%
  group_by(corrBlock, year, mutRate, corrTime) %>%
  summarise(completed_replicates = n_distinct(replicate)) %>%
  filter(!(year<1)) %>%
  filter(!(corrBlock == 0 & corrTime == 1)) %>%
  ungroup()
completed_reps %>% print(n=Inf)

total_reps = completed_reps %>% filter(corrBlock == 0 & year == 46) %>%
  pull(completed_replicates)

tissue_size = 4900 # 70x70 basal layer = 4900
total_size = tissue_size*total_reps

## counts the number of tp53 cells within each replicate/year
## change to mut == 68 for only tp53 + FA
tp53clones <- alldat %>%
  filter(!(corrBlock == 0 & corrTime == 1)) %>% ## plots no correction case with corrBlock == 0
  filter(!(year < 1)) %>%
  filter(mut %in% c(68, 777)) %>% # tp53 + FA or tp53 + correction
  group_by(replicate, corrBlock, year, mutRate, tp53Clone) %>%
  #mutate(count = n()) %>%
  summarize(count = n(), .groups = 'drop') %>%
  filter(count >= threshold_vector) %>% ## lower limit of detection
  mutate(proportion = count/tissue_size)

### Fill in 0's for simulations that acquired no tp53 mutations 
all_years <- seq(1, 46, by = 5)
tp53clones_noCorrection <- tp53clones %>% filter(corrBlock == 0)

unique(tp53clones_noCorrection$mutRate)

expanded_tp53clones <- expand_grid(replicate = seq(0,total_reps-1), 
                                   corrBlock = 0, 
                                   year = all_years, 
                                   mutRate = unique(tp53clones_noCorrection$mutRate))

# Merge with the original data and fill missing values with 0
tp53clones_full <- expanded_tp53clones %>%
  left_join(tp53clones_noCorrection, by = c("replicate", "corrBlock", "year", "mutRate")) %>%
  replace_na(list(count = 0, proportion = 0, tp53Clone = -1))

### randomly sample lost simulations from the no correction case
missing_reps <- completed_reps %>% filter(corrBlock != 0) %>%
  mutate(missing_replicates = total_reps - completed_replicates)


results <- list()
for (i in 1:nrow(missing_reps)){
  filter_corrBlock <- missing_reps$corrBlock[i]
  filter_year <- missing_reps$year[i]
  missing_replicates <- missing_reps$missing_replicates[i]
  
  sampled_data <- tp53clones_full %>%
    filter(year == filter_year) %>%
    select(replicate, tp53Clone, count)
    #sample(missing_replicates, replace=TRUE)
  
  sampled_replicates <- sample(0:(total_reps-1), missing_replicates, replace = TRUE)
    
  expanded_data <- sampled_data %>%
    filter(replicate %in% sampled_replicates) %>%
    group_by(replicate) %>%
    mutate(
      corrBlock = filter_corrBlock,
      year = filter_year
    ) %>%
    ungroup()
  
  # Append results to the list
  results[[i]] <- expanded_data
  
  
    # Append results to the list
  # results[[i]] <- data.frame(
  #   corrBlock = rep(filter_corrBlock, missing_replicates),
  #   year = rep(filter_year, missing_replicates),
  #   count = sampled_count
  # )
}

sampled_dataframe <- bind_rows(results)
tp53clones_sampled <- bind_rows(tp53clones, sampled_dataframe)

write.csv(tp53clones_sampled, glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/tp53clones_sampled_baseline.csv"), row.names = FALSE)

cat('tp53clones_sampled')
print(tp53clones_sampled, n = Inf)

tp53clones_sampled_46 <- tp53clones_sampled %>% filter(year == 46) %>%
  filter(tp53Clone != -1) %>%
  group_by(corrBlock) %>%
  summarise(total_tp53Clones = n(),
            total_count = sum(count, na.rm = TRUE), 
            total_proportion = total_count/total_size,
            .groups = 'drop')

cat("tp53clones_sampled_46")
print(tp53clones_sampled_46)

write.csv(tp53clones_sampled_46, glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/tp53clones_sampled_46_baseline.csv"), row.names = FALSE)

tp53clones_sampled_summarized <- tp53clones_sampled %>%
  filter(tp53Clone != -1) %>%
  group_by(corrBlock, year) %>%
  summarize(total_count = sum(count),
            total_tp53Clones = n(), 
            total_proportion = total_count/total_size, .groups = 'drop') %>%
  mutate(total_proportion = total_count/total_size)

write.csv(tp53clones_sampled_summarized, glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/tp53clones_sampled_summarized_baseline.csv"), row.names = FALSE)

q()







# tp53clones_sampled_summarized <- tp53clones_sampled %>%
#   group_by(corrBlock, year) %>%
#   summarize(total_count = sum(count), .groups = 'drop') %>%
#   mutate(total_proportion = total_count/total_size)
# tp53clones_sampled_summarized
# 
# ## if no mutations at a given timepoint, set proportion to 0
# 
# combinations <- expand.grid(
#   corrBlock = unique(tp53clones_sampled_summarized$corrBlock),
#   #mutRate = unique(tp53clones_sampled_summarized$mutRate),
#   year = all_years
# )
# 
# mean_proportions <- combinations %>%
#   left_join(tp53clones_sampled_summarized, by = c("corrBlock", "year")) %>%
#   mutate(
#     total_proportion = ifelse(is.na(total_proportion), 0, total_proportion),
#     total_count = ifelse(is.na(total_count), 0, total_count),
#   )
# mean_proportions
# 
# ## plot the mean proportion
# muted <- c('#bbe4dd', '#91d3c8', '#67c2b3', '#44aa99', '#338073', '#22564d', '#122c28')
# names(muted) <- c("0", "0.001", "0.01", "0.1", "0.2", "0.5", "1")
# custom_breaks <- seq(0,46,by=10)
# corrBlock_labels <- c(
#   "0" = "No correction",
#   "0.01" = "0.01",
#   "0.1" = "0.1"
# )
# legend_labels <- c("0" = "No correction", "0.01" = "0.01", "0.1" = "0.1")
# 
# df_labs <- data.frame(
#   corrBlock = c(0, 0.01, 0.1),  # The unique values of corrBlock
#   #corrBlock_label = c("No~correction", "italic(p) == 0.01", "italic(p) == 0.1")  # Math expressions
#   corrBlock_label = c("No~correction", "italic(p)[italic(FANC)^'+']^italic(TP53)^'+' == 0.01", "italic(p)[italic(FANC)^'+']^italic(TP53)^'+' == 0.1")
# )
# mean_proportions <- merge(mean_proportions, df_labs, by = "corrBlock", all.x = TRUE)
# mean_proportions$corrBlock_label <- factor(mean_proportions$corrBlock_label,
#                                            levels = c("No~correction", "italic(p)[italic(FANC)^'+']^italic(TP53)^'+' == 0.01", "italic(p)[italic(FANC)^'+']^italic(TP53)^'+' == 0.1"))
# 
# 
# plot <- ggplot(mean_proportions, aes(x = year, y = total_proportion, color = as.factor(corrBlock), group = corrBlock)) +
#   geom_line(size = 1) +
#   facet_wrap(~corrBlock_label, scales = "fixed", nrow = 1, labeller = label_parsed) +
#   labs(x = "Year", y = "Tissue proportion", color = bquote(italic(p) [italic(FANC)^"+"] ^italic(TP53)^"+")) +
#   scale_color_manual(values = muted, labels = legend_labels) +
#   theme_minimal() +
#   #scale_x_continuous(limits = c(1, 46), breaks = custom_breaks) +
#   #scale_x_continuous(breaks = seq(1, 50, 5)) +
#   scale_x_continuous(breaks = c(1, 10, 20, 30, 40, 50)) +
#   scale_y_continuous(limits = c(0, 0.4), breaks = seq(0,0.4, by = 0.1)) +
#   #theme(legend.position = c(0.2, 0.55),
#   theme(legend.position = c(0.1, 0.55),
#   legend.title = element_text(face = "italic", size=8)) +
#   theme(
#     #legend.text=element_text(size=8),
#     legend.key.size = unit(0.5, "lines"),
#     legend.key.height = unit(0.5, "lines"),
#     legend.margin = margin(0, 0, 0, 0),
#     legend.key.spacing.y = unit(0.0001, "cm"),
#     legend.key.spacing.x = unit(0.0001, "cm"),
#     
#     strip.text = element_text(size = 8, color='black'),
#     axis.text = element_text(size = 8, color='black'),
#     axis.title = element_text(size = 8, color='black'), 
#     legend.text = element_text(size = 8, color='black'), 
#     legend.title = element_text(size = 8, color='black'),
#     axis.line = element_line(color="black"),
#     axis.ticks = element_line(color = "black"),
#     plot.margin = margin(0.5, 25.5, 0.5, 5.5), 
#     panel.spacing = unit(1.5, "lines")
#   )
# show(plot)
# 
# if(cluster){
#   ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/tp53_tissueProportion.png"), plot, width = 6, height = 1.35, units = "in", dpi = 300)
# }else{
#   ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/tp53_tissueProportion.png"), plot, width = 6, height = 1.35, units = "in", dpi = 300)
# }
# 
# 
# 
# 
# 
# 
# 
# 
# q()
# ### histogram figure showing persistence of each experiment
# 
# ## experiment 1 and 2:
# p_data <- tibble(
#   TP53 = factor(c("TP53+", "TP53+", "TP53-", "TP53-", "TP53+", "TP53+", "TP53-", "TP53-"),
#                 levels = c("TP53+", "TP53-")),
#   FA = c("FANC-", "FANC-", "FANC-", "FANC-", "FANC+", "FANC+", "FANC+", "FANC+"),
#   Category = c("p_correction", "p_tp53_FA", "p_correction", "p_tp53_FA", "p_correction", "p_tp53_FA", "p_correction", "p_tp53_corr"),
#   Value = c(0, 0, 0, 0.01, 0.1, 0, 0.1, 0.01),
#   Position = c(0, 0, 0, 0, 0, 0.1, 0, 0.1)
# )
# 
# hist1 <- ggplot(p_data, aes(x = TP53, y = Value, fill = Category)) +
#   geom_bar(stat = "identity", position = "stack") +
#   facet_wrap(~ FA, nrow = 1, strip.position = "bottom") +
#   scale_fill_manual(values = c("p_correction" = "#67c2b3", "p_tp53_FA" = "black", "p_tp53_corr" = "black")) +
#   labs(
#     x = NULL,
#     y = "Persistence"
#   ) +
#   theme_minimal() +
#   ylim(0,0.15) +
#   theme(
#     strip.background = element_blank(),
#     axis.line = element_line(color="black"),
#     strip.text = element_text(size = 8, face = "italic", margin = margin(b=0.2)),
#     axis.text.x = element_text(angle = 45, hjust = 1, face = 'italic'),
#     axis.title.y = element_text(size = 8),
#     axis.title.x = element_blank(),
#     panel.spacing = unit(0.1, "lines"),
#     strip.placement = "outside", 
#     legend.position = "none", 
#     axis.text = element_text(size = 8, color='black')
#   )  +
#   geom_text(
#     data = subset(p_data, Category == "p_correction" & FA == "FANC+"), 
#     aes(label = "p", y = Value / 2),  # Adjust 'Label' and positioning as needed
#     color = "white",
#     size = 8/.pt,
#     vjust = 0.5,
#     fontface = "italic"
#   )
# show(hist1)
# ## experiment 3:
# p_data <- tibble(
#   TP53 = factor(c("TP53+", "TP53+", "TP53-", "TP53-", "TP53-", "TP53+", "TP53+", "TP53-", "TP53-"),
#                 levels = c("TP53+", "TP53-")),
#   FA = c("FANC-", "FANC-", "FANC-", "FANC-", "FANC-", "FANC+", "FANC+", "FANC+", "FANC+"),
#   Category = factor(c("p_correction", "p_tp53", "p_correction", "p_tp53_FA", "p_tp53", "p_correction", "p_tp53", "p_correction", "p_tp53"),
#                 levels = c("p_tp53_FA", "p_correction", "p_tp53")),  # Define the stacking order
#   Value = c(0, 0, 0, 0.01, 0.01, 0.1, 0, 0.1, 0.01),
#   Position = c(0, 0, 0, 0, 0.01, 0, 0.1, 0, 0.1)
# )
# p_data
# hist2 <- ggplot(p_data, aes(x = TP53, y = Value, fill = Category)) +
#   geom_bar(stat = "identity", position = "stack") +
#   facet_wrap(~ FA, nrow = 1, strip.position = "bottom") +
#   scale_fill_manual(values = c( "p_tp53_FA" = "#CC6677", "p_tp53" = "black", "p_correction" = "#67c2b3")) +
#   labs(
#     x = NULL,
#     y = "Persistence"
#   ) +
#   theme_minimal() +
#   ylim(0,0.15) +
#   theme(
#     strip.background = element_blank(),
#     axis.line = element_line(color="black"),
#     strip.text = element_text(size = 8, face = "italic", margin = margin(b=0.2)),
#     axis.text.x = element_text(size = 8, angle = 45, hjust = 1, face = 'italic'),
#     axis.title.y = element_text(size = 8),
#     axis.title.x = element_blank(),
#     panel.spacing = unit(0.1, "lines"),
#     strip.placement = "outside", 
#     legend.position = "none", 
#     axis.text = element_text(size = 8, color='black')
#   ) +
#   geom_text(
#     data = subset(p_data, Category == "p_correction" & FA == "FANC+"), 
#     aes(label = "p", y = Value / 2),  # Adjust 'Label' and positioning as needed
#     color = "white",
#     size = 8/.pt,
#     vjust = 0.5,
#     fontface = "italic"
#   )
# show(hist2)
# height = 1.25
# width = 1.5
# if(cluster){
#   ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/histogram_1.png"), hist1, width = width, height = height, units = "in", dpi = 300)
#   ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/histogram_2.png"), hist2, width = width, height = height, units = "in", dpi = 300)
#   }else{
#   ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/histogram_1.png"), hist1, width = width, height = height, units = "in", dpi = 300)
#   ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/histogram_2.png"), hist2, width = width, height = height, units = "in", dpi = 300)
#   }
# 
# 
# 
# 
# 
# # ## filter data based on detection thresholds and combine all dataframes
# # df_list <- list()
# # for (val in threshold_vector) {
# #   duplicated_df <- tp53clones %>%
# #     filter(cells_count >= val) %>%
# #     mutate(threshold = val)
# #   df_list <- append(df_list, list(duplicated_df))
# # }
# # tp53_thresholds <- do.call(rbind, df_list)
# 
# 
# ## combine tp53 mut counts with the number of replicates completed in experiment
# unique_tp53Clone_reps <- completed_reps %>%
#   left_join(unique_tp53Clone, by = c("corrBlock", "year"))
# unique_tp53Clone_reps
# ## average mutations per replicate
# ave_unique_count <- unique_tp53Clone_reps %>%
#   mutate(ave_unique_count = total_unique_count/completed_replicates) %>%
#   drop_na()
# print(ave_unique_count)
# ## Plot using ggplot
# plot <- ggplot(ave_unique_count, aes(x = year, y = ave_unique_count, color = factor(corrBlock))) +
#   geom_line() +
#   labs(x = "Year", y = "Average number of tp53 mutations per replicate", color = "Blocking Probability b") +
#   theme_minimal() + 
#   #theme(aspect.ratio = 0.75) + 
#   scale_x_continuous(breaks = seq(1, 50, 5)) +
#   facet_wrap(~threshold, scales = "fixed", nrow = 1) +
#   theme(text = element_text(size = 20))
#   #theme(strip.text = element_text(size = 20))
# 
# show(plot)
# ## save plot to file
# if(cluster){
#   ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/tp53_mutations_per_simulation_thresholds.png"), plot, width = 20, height = 8, units = "in", dpi = 300)
# }else{
#   ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/tp53_mutations_per_simulation_thresholds.png"), plot, width = 20, height = 8, units = "in", dpi = 300)
# }
# 
# 
# ############
# ##### Fill missing values
# ############
# 
# reps_num = 100
# reps_num = 300 # use if combining all experiments
# 
# ### Determine the missing replicates to be added.
# summary_data <- alldat %>%
#   filter(!(year < 1)) %>%
#   filter(!(corrBlock == 0 & corrTime == 1)) %>%
#   group_by(corrBlock, year, replicate) %>%
#   summarize(observations = n())
# #summary_data %>% print(n=1000)
# 
# years_total = unique(summary_data$year)
# 
# max_replicates = n_distinct(unique(summary_data$replicate[summary_data$corrBlock ==0]))
# 
# missing_replicates <- lapply(years_total, function(year_value) {
#   summary_data %>%
#     filter(year == year_value) %>%
#     group_by(corrBlock) %>%
#     print() %>%
#     summarize(available_replicates = n_distinct(replicate)) %>%
#     mutate(missing_replicates = max_replicates - available_replicates) %>%
#     mutate(missing_replicates = ifelse(missing_replicates < 0, 0, missing_replicates)) %>%
#     mutate(year = year_value)
# })
# missing_replicates <- do.call(rbind, missing_replicates)
# print(missing_replicates, n=Inf)
# missing_replicates
# missing_replicates <- missing_replicates %>% filter(missing_replicates >0)
# 
# tp53clones <- alldat %>%
#   filter(!(corrBlock == 0 & corrTime == 1)) %>% ## plots no correction case with corrBlock == 0
#   filter(!(year < 1)) %>%
#   dplyr::select(tp53Clone, replicate, corrBlock, year, mut) %>%
#   filter(mut %in% c(68, 777)) %>% # tp53 + FA or tp53 + correction
#   #filter(mut ==68) %>% # tp53 + FA only (no correction)
#   group_by(tp53Clone, replicate, corrBlock, year, mut) %>%
#   summarise(cells_count = n())
# tp53clones %>% print(n=Inf)
# tp53clones
# 
# missing_replicates
# for (i in 1:nrow(missing_replicates)){
#   n_missing <- missing_replicates$missing_replicates[i]
#   corrBlock_replace <- missing_replicates$corrBlock[i]
#   year <- missing_replicates$year[i]
#   for (j in 1:n_missing){
#     n_reps = reps_num
#     random_rep = sample(0:n_reps-1, size = 1)
#     sampled_rep <- tp53clones %>%
#       filter(corrBlock==0, year == year, replicate == random_rep)
#     sampled_rep <- mutate(sampled_rep, corrBlock = corrBlock_replace)
# 
#     tp53clones <- bind_rows(tp53clones, sampled_rep)
#     
#   }
# }
# tp53clones
# 
# unique_reps <- tp53clones %>%
#   group_by(corrBlock, year) %>%
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
# completed_reps <- mutate(completed_reps, completed_replicates = reps_num)
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
#   group_by(corrBlock, year, replicate, threshold) %>%
#   summarise(unique_count = n_distinct(tp53Clone)) %>%
#   ungroup() %>%
#   group_by(corrBlock, year, threshold) %>%
#   summarise(total_unique_count = sum(unique_count)) %>%
#   ungroup()
# unique_tp53Clone %>% print(n=Inf)
# 
# ## remove columns from completed dataframe
# completed_reps <- completed_reps %>%
#   dplyr::select(-mutRate, -corrTime)
# completed_reps %>% print(n=Inf)
# 
# ## combine tp53 mut counts with the number of replicates completed in experiment
# unique_tp53Clone_reps <- completed_reps %>%
#   left_join(unique_tp53Clone, by = c("corrBlock", "year"))
# 
# ## average mutations per replicate
# ave_unique_count <- unique_tp53Clone_reps %>%
#   mutate(ave_unique_count = total_unique_count/completed_replicates) %>%
#   drop_na()
# print(ave_unique_count)
# ave_unique_count <- ave_unique_count %>%
#   filter(threshold == 108)
# 
# #### Add 0's to data to show up on plot even if no mutations
# #### must manually check data first if 0s can be added.
# start_years <- c(1, 6, 11)
# start_years <- c(1)
# unique_corrBlock <- unique(ave_unique_count$corrBlock)
# unique_threshold <- unique(ave_unique_count$threshold)
# 
#  # Create a dataframe with the additional rows
# additional_rows <- expand.grid(year = start_years,
#                                 ave_unique_count = 0,
#                                 corrBlock = unique_corrBlock,
#                                 threshold = unique_threshold)
# 
# #Combine with the original dataframe
# ave_unique_count_extended <- bind_rows(ave_unique_count, additional_rows) %>%
#   arrange(year, corrBlock, threshold)
# print(ave_unique_count_extended, n=Inf)
# 
# ## Plot using ggplot
# num_colors <- length(unique(ave_unique_count$corrBlock))
# my_palette <- rev(viridis(num_colors))
# my_palette[1] <- "#BADE28FF"
# plot <- ggplot(ave_unique_count_extended, aes(x = year, y = ave_unique_count, color = factor(corrBlock))) +
#   geom_line(size = 3) +
#   ylim(0, 0.45) +
#   #labs(x = "Year", y = "TP53 mutation frequency", color = "p") +
#   labs(x = "Year", y = expression(paste(italic("TP53"), " mutation frequency")), color = "p") +
#   theme_minimal() + 
#   scale_color_manual(values = my_palette, 
#                      breaks = unique(ave_unique_count$corrBlock), 
#                      labels = unique(ave_unique_count$corrBlock)) +
#   #theme(aspect.ratio = 0.75) + 
#   scale_x_continuous(breaks = seq(1, 50, 5)) +
#   facet_wrap(~threshold, scales = "fixed", nrow = 1) +
#   #theme(legend.position = "top") +
#   theme(legend.position = c(0.1, 0.75)) +
#   theme(
#     axis.line = element_line(color="black"), 
#     axis.text = element_text(color="black", size=28),
#     axis.title.x = element_text(color="black", size=28),
#     axis.title.y = element_text(color="black", size=28),
#     axis.ticks = element_line(color = "black"),
#     legend.text = element_text(color="black", size=28), 
#     legend.title = element_text(color="black", face="italic", size=28),
#     strip.background = element_blank(),
#     strip.text.x = element_blank(),
#   )
# 
# 
# show(plot)
# 
# ## save plot to file
# if(cluster){
#   ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/tp53_mutations_per_simulation_thresholds_color.png"), plot, width = 10, height = 6, units = "in", dpi = 300)
# }else{
#   ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/tp53_mutations_per_simulation_thresholds_poster.png"), plot, width = 10, height = 6, units = "in", dpi = 300)
# }
# 
# 
