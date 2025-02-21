## Hunter Colegrove
## FA dynamics
##
## Reads in clone data from tp53_stats.R
## 


require(foreach)
require(tidyverse)
require(glue)
library(gridExtra)
require(grid)
library(imager)
library(gtools)
library(png)
library(ggtext)


###############################################################################
####### total proportion
###############################################################################

#########################
####### baseline
#########################

project_path <- "2024-07-04-tp53_competition_sameMut_sameP"
tp53clones_sampled_summarized <- mutRate_df <- read_csv("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-07-04-tp53_competition_sameMut_sameP/dat/tp53clones_sampled_summarized_baseline.csv")


combinations <- expand.grid(
  corrBlock = unique(tp53clones_sampled_summarized$corrBlock),
  #mutRate = unique(tp53clones_sampled_summarized$mutRate),
  year = all_years
)

mean_proportions <- combinations %>%
  left_join(tp53clones_sampled_summarized, by = c("corrBlock", "year")) %>%
  mutate(
    total_proportion = ifelse(is.na(total_proportion), 0, total_proportion),
    total_count = ifelse(is.na(total_count), 0, total_count),
  )
mean_proportions

## plot the mean proportion
muted <- c('#bbe4dd', '#91d3c8', '#67c2b3', '#44aa99', '#338073', '#22564d', '#122c28')
names(muted) <- c("0", "0.001", "0.01", "0.1", "0.2", "0.5", "1")
custom_breaks <- seq(0,46,by=10)
corrBlock_labels <- c(
  "0" = "No correction",
  "0.01" = "0.01",
  "0.1" = "0.1"
)
legend_labels <- c("0" = "No correction", "0.01" = "0.01", "0.1" = "0.1")

df_labs <- data.frame(
  corrBlock = c(0, 0.01, 0.1),  # The unique values of corrBlock
  #corrBlock_label = c("No~correction", "italic(p) == 0.01", "italic(p) == 0.1")  # Math expressions
  corrBlock_label = c("No~correction", "italic(p)[italic(FANC)^'+']^italic(TP53)^'+' == 0.01", "italic(p)[italic(FANC)^'+']^italic(TP53)^'+' == 0.1")
)
mean_proportions <- merge(mean_proportions, df_labs, by = "corrBlock", all.x = TRUE)
mean_proportions$corrBlock_label <- factor(mean_proportions$corrBlock_label,
                                           levels = c("No~correction", "italic(p)[italic(FANC)^'+']^italic(TP53)^'+' == 0.01", "italic(p)[italic(FANC)^'+']^italic(TP53)^'+' == 0.1"))

plot <- ggplot(mean_proportions, aes(x = year, y = total_proportion, color = as.factor(corrBlock), group = corrBlock)) +
  geom_line(size = 1) +
  facet_wrap(~corrBlock_label, scales = "fixed", nrow = 1, labeller = label_parsed) +
  #labs(x = "Year", y = "Tissue proportion", color = bquote(italic(p) [italic(FANC)^"+"] ^italic(TP53)^"+")) +
  labs(x = "Year", y = "Tissue proportion", color = bquote(italic(p)[corr])) +
  scale_color_manual(values = muted, labels = legend_labels) +
  theme_minimal() +
  #scale_x_continuous(limits = c(1, 46), breaks = custom_breaks) +
  #scale_x_continuous(breaks = seq(1, 50, 5)) +
  scale_x_continuous(breaks = c(1, 10, 20, 30, 40, 50)) +
  scale_y_continuous(limits = c(0, 0.4), breaks = seq(0,0.4, by = 0.1)) +
  #theme(legend.position = c(0.2, 0.55),
  theme(legend.position = c(0.1, 0.6),
        legend.title = element_text(face = "italic", size=8)) +
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
    plot.margin = margin(0.5, 2, 0.5, 5.5), 
    panel.spacing = unit(1.5, "lines")
  )
show(plot)

plot_noLabels <- plot + 
  theme(strip.text = element_blank(), 
        axis.title.y = element_text(hjust = 0.75),
        plot.margin = margin(10, 2, 0.5, 5.5), 
        legend.title = element_text(margin = margin(b = 0))
  )

if(cluster){
  ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/tp53_tissueProportion_stats.png"), plot, width = 5.8, height = 1.35, units = "in", dpi = 300)
}else{
  ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/tp53_tissueProportion_stats.png"), plot, width = 5.8, height = 1.35, units = "in", dpi = 300)
  ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/tp53_tissueProportion_noFacetLabels_stats.png"), plot_noLabels, width = 5.8, height = 1, units = "in", dpi = 300)
  
}


#########################
####### mutRate
#########################

project_path <- "2024-06-13-tp53_competition_changeMut_sameP"
tp53clones_sampled_summarized <- mutRate_df <- read_csv("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-06-13-tp53_competition_changeMut_sameP/dat/tp53clones_sampled_summarized_mutRate.csv")


combinations <- expand.grid(
  corrBlock = unique(tp53clones_sampled_summarized$corrBlock),
  mutRate = unique(tp53clones_sampled_summarized$mutRate),
  year = all_years
)

mean_proportions <- combinations %>%
  left_join(tp53clones_sampled_summarized, by = c("corrBlock", "year", "mutRate")) %>%
  mutate(
    total_proportion = ifelse(is.na(total_proportion), 0, total_proportion),
    total_count = ifelse(is.na(total_count), 0, total_count),
  )
#mean_proportions

cat("mean_proportions\n")
print(mean_proportions)

## plot the mean proportion
muted <- c('#bbe4dd', '#91d3c8', '#67c2b3', '#44aa99', '#338073', '#22564d', '#122c28')
names(muted) <- c("0", "0.001", "0.01", "0.1", "0.2", "0.5", "1")

custom_labels <- c("0.375" = "10.6X", "0.5" = "8X", "1" = "4X", "2" = "2X", "2.666" = "1.5X")
blue_shades <- c("0.375" = "#012A4A", "0.5" = "#01497c", "1" = "#2a6f97", "2" = "#468faf", "2.666" = "#89c2d9")

custom_breaks <- seq(0,46,by=10)
# corrBlock_labels <- c(
#   "0" = "No correction",
#   "0.01" = bquote(italic(p) ~ "= 0.01"),
#   "0.1" = bquote(italic(p) ~ "= 0.1")
# )


## for line labels
dat_label <- mean_proportions %>%
  group_by(corrBlock, mutRate) %>%
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

plot <- ggplot(mean_proportions, aes(x = year, y = total_proportion, color = as.factor(mutRate), group = mutRate)) +
  geom_line(size = 1) +
  facet_wrap(~corrBlock_label, scales = "fixed", nrow = 1, labeller = label_parsed) +  # Using label_parsed
  labs(x = "Year", y = "Tissue proportion", color = bquote(mu * phantom(x) [italic(FANC)^"-"] ~ "increase")) +
  scale_color_manual(values = blue_shades, labels = custom_labels) +
  theme_minimal() +
  scale_x_continuous(breaks = c(1, 10, 20, 30, 40, 50)) +
  scale_y_continuous(limits = c(0, 0.4), breaks = seq(0,0.4, by = 0.1)) +
  theme(legend.position = c(0.1, 0.7))+
  coord_cartesian(clip = "off") +
  theme(
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
    plot.margin = margin(0.5, 2, 0.5, 5.5), 
    panel.spacing = unit(1.5, "lines")
  )

show(plot)

plot_noLabels <- plot + 
  theme(strip.text = element_blank(), 
        axis.title.y = element_text(hjust = 0.75),
        plot.margin = margin(10, 2, 0.5, 5.5), 
        legend.title = element_text(margin = margin(b = 0))
  )

if(cluster){
  ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/tp53_tissueProportion_stat.png"), plot, width = 5.8, height = 1.35, units = "in", dpi = 300)
  ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/tp53_tissueProportion_noFacetLabels_stat.png"), plot_noLabels, width = 6, height = 1.15, units = "in", dpi = 300)
  
}else{
  ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/tp53_tissueProportion_stat.png"), plot, width = 5.8, height = 1.35, units = "in", dpi = 300)
  ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/tp53_tissueProportion_noFacetLabels_stat.png"), plot_noLabels, width = 5.8, height = 1, units = "in", dpi = 300)
  
}

#########################
####### tp53 persistence
#########################

tp53clones_sampled_summarized <- mutRate_df <- read_csv("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-07-24-competition_sameMut_changeP/dat/tp53clones_sampled_summarized_p53Block.csv")
tp53clones_sampled_summarized

combinations <- expand.grid(
  corrBlock = unique(tp53clones_sampled_summarized$corrBlock),
  FAp53block = unique(tp53clones_sampled_summarized$FAp53block),
  year = all_years
)

mean_proportions <- combinations %>%
  left_join(tp53clones_sampled_summarized, by = c("corrBlock", "year", "FAp53block")) %>%
  mutate(
    total_proportion = ifelse(is.na(total_proportion), 0, total_proportion),
    total_count = ifelse(is.na(total_count), 0, total_count),
  )

custom_labels <- c( "0.04" = "4X", "0.02" = "2X", "0.015" = "1.5X")
red_shades <- c("0.04" = "#882255", "0.02" = "#CC6677", "0.015" = "#EE99AA")

custom_breaks <- seq(0,46,by=10)

## for line labels
dat_label <- mean_proportions %>%
  group_by(corrBlock, FAp53block) %>%
  filter(year == max(year)) %>%
  ungroup()

df_labs <- data.frame(
  corrBlock = c(0, 0.01, 0.1),  # The unique values of corrBlock
  corrBlock_label = c("No~correction", "italic(p)[italic(FANC)^'+']^italic(TP53)^'+' == 0.01", "italic(p)[italic(FANC)^'+']^italic(TP53)^'+' == 0.1")  # Math expressions
)
mean_proportions <- merge(mean_proportions, df_labs, by = "corrBlock", all.x = TRUE)
mean_proportions$corrBlock_label <- factor(mean_proportions$corrBlock_label,
                                           levels = c("No~correction", "italic(p)[italic(FANC)^'+']^italic(TP53)^'+' == 0.01", "italic(p)[italic(FANC)^'+']^italic(TP53)^'+' == 0.1"))
mean_proportions$FAp53block <- factor(mean_proportions$FAp53block, 
                                      levels = c("0.04", "0.02", "0.015"))  # Order by 4X, 2X, 1.5X

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
    plot.margin = margin(0.5, 2, 0.5, 5.5), 
    panel.spacing = unit(1.5, "lines")
  )
show(plot)

plot_noLabels <- plot + 
  theme(strip.text = element_blank(), 
        axis.title.y = element_text(hjust = 0.75),
        plot.margin = margin(10, 2, 0.5, 5.5), 
        legend.title = element_text(margin = margin(b = 0))
  )

show(plot_noLabels)
if(cluster){
  ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-07-24-competition_sameMut_changeP/results/tp53_tissueProportion.png"), plot, width = 5.8, height = 1.35, units = "in", dpi = 300)
}else{
  ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-07-24-competition_sameMut_changeP/results/tp53_tissueProportion.png"), plot, width = 5.8, height = 1.35, units = "in", dpi = 300)
  ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-07-24-competition_sameMut_changeP/results/tp53_tissueProportion_noFacetLabels.png"), plot_noLabels, width = 5.8, height = 1, units = "in", dpi = 300)
  
}

###############################################################################
####### In text values
###############################################################################

baseline_df <- read_csv("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-07-04-tp53_competition_sameMut_sameP/dat/tp53clones_sampled_46_baseline.csv")
mutRate_df <- read_csv("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-06-13-tp53_competition_changeMut_sameP/dat/tp53clones_sampled_46_mutRate.csv")
tp53Block_df <- read_csv("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-07-24-competition_sameMut_changeP/dat/tp53clones_sampled_46_p53Block.csv")

## mutRate changes at 46 years
## tissue proportion change from baseline with 8x increase

# In the absence of gene therapy, elevated mutation rates in FANC- tissue 
# significantly increased the tissue proportion with a TP53- clone after 
# 46 years from XXX with background mutation rate to YYY with an 8x 
# elevated mutation rate (Figure 4F).

baseline_noCorrection <- baseline_df %>% 
  filter(corrBlock == 0) %>% 
  pull(total_proportion)

mutRate_noCorrection <- mutRate_df %>%
  filter(corrBlock == 0) %>%
  filter(mutRate == 0.5) %>%
  pull(total_proportion)

# Introducing gene correction with a small persistence coefficient 
# (pcorr = 0.01, Figure 4G) reduced TP53- tissue proportion moderately 
# (~XXX% reduction at an 8x mutation rate) compared to no correction.

mutRate_0.01 <- mutRate_df %>%
  filter(corrBlock == 0.01) %>%
  filter(mutRate == 0.5) %>%
  pull(total_proportion)

percent_change_mut <- (mutRate_noCorrection - mutRate_0.01) / mutRate_noCorrection * 100


# TP53 proliferative changes at 46 years
# In the absence of gene therapy, when pFANC-TP53- was increased by a factor of 
# four, this increased persistence of TP53- cells and led to more extensive 
# tissue coverage: TP53- mutants covered XXX% of tissue on average at 46 years 
# compared to only YYY% of tissue under background persistence advantage.

tp53Block_noCorrection <- tp53Block_df %>%
  filter(corrBlock == 0) %>%
  filter(FAp53block == 0.04) %>%
  pull(total_proportion)

baseline_noCorrection

# When the persistence advantage of gene correction was small 
# (pcorr = 0.01, Figure 4J), reduced average TP53- tissue proportion 
# (~XXX% at a 4x proliferative advantage) compared to no correction.

tp53Block_0.01 <- tp53Block_df %>%
  filter(corrBlock == 0.01) %>%
  filter(FAp53block == 0.04) %>%
  pull(total_proportion)

tp53Block_0.01

percent_change_block <- (tp53Block_noCorrection - tp53Block_0.01) / tp53Block_noCorrection * 100
percent_change_block


fileConn<-file("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/tp53_figure_stats.txt")

writeLines(c(
  "## mutRate changes at 46 years",
  paste("In the absence of gene therapy, elevated mutation rates in FANC- tissue significantly increased the tissue proportion with a TP53- clone after 46 years from", 
        sprintf("%.3f", baseline_noCorrection), 
        "with background mutation rate to", 
        sprintf("%.3f", mutRate_noCorrection), 
        "with an 8x elevated mutation rate (Figure 4F)."),
  paste("Introducing gene correction with small persistence coefficients (p = 0.01) reduced TP53- prevalence moderately (~", 
        sprintf("%.2f%%", percent_change_mut), 
        " reduction at an 8x mutation rate) compared to no correction."),
  "",
  "## tp53Block changes at 46 years",
  paste("In the absence of gene therapy, when pFANC-TP53- was increased by a factor of four, this increased persistence of TP53- cells and led to more extensive tissue coverage: TP53- mutants covered", 
        sprintf("%.1f%%", tp53Block_noCorrection * 100), 
        "of tissue on average at 46 years compared to only",
        sprintf("%.3f", baseline_noCorrection), 
        " of tissue under background persistence advantage."),
  paste("When the persistence advantage of gene correction was small (pcorr = 0.01, Figure 4J), reduced average TP53- tissue proportion (~", 
        sprintf("%.3f%%", percent_change_block), 
        "when TP53 mutations in FANC- backgrounds had no added persistence advantage.")
), fileConn)

close(fileConn)

###############################################################################
####### Supplemental Figure
###############################################################################

baseline_df <- read_csv("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-07-04-tp53_competition_sameMut_sameP/dat/tp53clones_sampled_46_baseline.csv")
mutRate_df <- read_csv("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-06-13-tp53_competition_changeMut_sameP/dat/tp53clones_sampled_46_mutRate.csv")
tp53Block_df <- read_csv("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-07-24-competition_sameMut_changeP/dat/tp53clones_sampled_46_p53Block.csv")

## filter and combine dfs
mutRate_df
baseline_df <- baseline_df %>% 
  mutate(FAp53block = 0.01) %>% 
  mutate(mutRate = 4) %>%
  mutate(clones_per_section = total_tp53Clones/300) %>%
  mutate(ave_clone_size = total_count/total_tp53Clones)

mutRate_df <- mutRate_df %>% 
  mutate(FAp53block = 0.01) %>%
  mutate(clones_per_section = total_tp53Clones/100)%>%
  mutate(ave_clone_size = total_count/total_tp53Clones)

tp53Block_df <- tp53Block_df %>% 
  mutate(mutRate = 4) %>%
  mutate(clones_per_section = total_tp53Clones/100)%>%
  mutate(ave_clone_size = total_count/total_tp53Clones)

all_clones <- baseline_df %>% bind_rows(mutRate_df) %>% bind_rows(tp53Block_df)

df_labs <- data.frame(
  corrBlock = c(0, 0.01, 0.1),
  corrBlock_label = c("No~correction", "italic(p)[corr] == 0.01", "italic(p)[corr] == 0.1")
)
all_clones <- merge(all_clones, df_labs, by = "corrBlock", all.x = TRUE)
all_clones$corrBlock_label <- factor(all_clones$corrBlock_label,
                                           levels = c("No~correction", "italic(p)[corr] == 0.01", "italic(p)[corr] == 0.1"))


all_clones

## set common themes and colors
text8pt_theme <- theme_minimal() + 
  theme(
    strip.text = element_text(size = 8, color = 'black'),
    axis.text = element_text(size = 8, color = 'black'),
    axis.title = element_text(size = 8, color = 'black'), 
    legend.text = element_text(size = 8, color = 'black'), 
    legend.title = element_text(size = 8, color = 'black'),
    legend.key.size = unit(0.3, "cm"),
    #legend.position = c(.85,.7),
    legend.position = 'none',
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
    plot.margin = margin(0.5, 2, 0.5, 5.5)
  )
p53_labels <- c("0.01" = "1X", "0.04" = "4X", "0.02" = "2X", "0.015" = "1.5X")
red_shades <- c("0.01" = "#44aa99", "0.04" = "#882255", "0.02" = "#CC6677", "0.015" = "#EE99AA")

mutation_labels <- c("4" = "1X", "0.5" = "8X", "1" = "4X", "2" = "2X", "2.666" = "1.5X")
blue_shades <- c("4" = "#44aa99", "0.5" = "#01497c", "1" = "#2a6f97", "2" = "#468faf", "2.666" = "#89c2d9")

#POS_WIDTH = 5
WIDTH = 0.75

### compare the number of clones based on mutation rate
compare_mut <- all_clones %>% filter(FAp53block == 0.01)

p1 <- ggplot(compare_mut, aes(x = factor(mutRate, levels = c(4, 2.666, 2, 1, 0.5)), y = clones_per_section, fill = factor(mutRate))) +
  geom_bar(stat = "identity", width=WIDTH) +
  facet_wrap(~ corrBlock_label, labeller = label_parsed) +
  scale_x_discrete(labels = mutation_labels) +
  scale_fill_manual(values = blue_shades, labels = mutation_labels) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(x = bquote(mu * phantom(x) [italic(FANC)^"-"] ~ "increase"), y = "Ave. # clones/sample", fill = bquote(mu * phantom(x) [italic(FANC)^"-"] ~ "increase")) +
  text8pt_theme
p1

### compare the number of clones based on tp53 persistence coefficient
compare_p53 <- all_clones %>% 
  filter(mutRate == 4) %>%
  mutate(FAp53block = factor(FAp53block, levels = rev(c(0.01, 0.015, 0.02, 0.04))))


p2 <- ggplot(compare_p53, aes(x = factor(FAp53block, levels = c(.01,0.015,0.02,0.04)), y = clones_per_section, fill = factor(FAp53block))) +
  geom_bar(stat = "identity", width=WIDTH) +
  facet_wrap(~ corrBlock_label, labeller = label_parsed) +
  scale_x_discrete(labels = p53_labels) +
  scale_fill_manual(values = red_shades, labels = p53_labels) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(x = bquote(italic(p) [italic(FANC)^"-"] ^italic(TP53)^"-" ~ "increase"), y = "Ave. # clones/sample", fill = bquote(italic(p) [italic(FANC)^"-"] ^italic(TP53)^"-" ~ "increase")) +
  text8pt_theme


### compare the average clone size based on mutation rate
compare_mut <- all_clones %>% filter(FAp53block == 0.01)

p3 <- ggplot(compare_mut, aes(x = factor(mutRate, levels = c(4, 2.666, 2, 1, 0.5)), y = ave_clone_size, fill = factor(mutRate))) +
  geom_bar(stat = "identity", width=WIDTH) +
  facet_wrap(~ corrBlock_label, labeller = label_parsed) +
  scale_x_discrete(labels = mutation_labels) +
  scale_fill_manual(values = blue_shades, labels = mutation_labels) +
  scale_y_continuous(limits = c(0, 4000), breaks = seq(0, 4000, by = 1000)) +
  labs(x = bquote(mu * phantom(x) [italic(FANC)^"-"] ~ "increase"), y = "Ave. clone size", fill = bquote(mu * phantom(x) [italic(FANC)^"-"] ~ "increase")) +
  text8pt_theme

### compare the average clone size based on tp53 persistence coefficient
compare_p53 <- all_clones %>% 
  filter(mutRate == 4) %>%
  mutate(FAp53block = factor(FAp53block, levels = rev(c(0.01, 0.015, 0.02, 0.04))))

p4 <- ggplot(compare_p53, aes(x = factor(FAp53block, levels = c(.01,0.015,0.02,0.04)), y = ave_clone_size, fill = factor(FAp53block))) +
  geom_bar(stat = "identity", width=WIDTH) +
  facet_wrap(~ corrBlock_label, labeller = label_parsed) +
  scale_x_discrete(labels = p53_labels) +
  scale_fill_manual(values = red_shades, labels = p53_labels) +
  scale_y_continuous(limits = c(0, 4000), breaks = seq(0, 4000, by = 1000)) +
  labs(x = bquote(italic(p) [italic(FANC)^"-"] ^italic(TP53)^"-" ~ "increase"), y = "Ave. clone size", fill = bquote(italic(p) [italic(FANC)^"-"] ^italic(TP53)^"-" ~ "increase")) +
  text8pt_theme

supplement <- plot_grid(p1,p2,p3,p4, labels = c("A","C","B","D"), label_size=14, label_x = 0.05)
show(supplement)
ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/tp53_supplement.png"), supplement, width = 5.5, height = 3, units = "in", dpi = 300)

