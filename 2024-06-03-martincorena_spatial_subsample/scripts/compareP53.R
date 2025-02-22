## Hunter Colegrove
## FA dynamics
## 
## Plot mutation frequency per over time and compare to Martincorena esophagus

## First, run subsectEsophagus.R to generate random subsections of martincorena data
## Run, countClones_tp53.R to generate tp53 clone sizes across a range of detection thresholds


cluster = FALSE
vis_mutRate = 4
lower_detection_limit = 108 # cells

require(tidyverse)
require(viridis)
require(glue)
require(fitdistrplus)
require(entropy)
require(cowplot)


###### project path of subsampled reconstructed tissue
subsample_path = "2024-06-03-martincorena_spatial_subsample"

###### project path of simulations
project_path = "2024-05-03-tp53_martincorena_mutRate_54years"


if(cluster){
  ## subsamples
  subsamplefilepath = glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{subsample_path}/dat/subsamples/mutations_subsamples.csv")
  
  ## tp53 unique clones across detection thresholds
  simfilepath = glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/countClones_tp53_R.csv")
  
  ## VAF path
  vaffilepath = glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/VAF_tp53_R.csv")
}else{
  ## subsamples
  subsamplefilepath = glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{subsample_path}/dat/subsamples/mutations_subsamples.csv")
  
  ## tp53 unique clones across detection thresholds
  simfilepath = glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/countClones_tp53_R.csv")
  
  ## VAF path
  vaffilepath = glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/VAF_tp53_R.csv")
}



###### read in files
tp53sim = read_delim(simfilepath, delim = ",")
tp53subsample = read_delim(subsamplefilepath, delim = ",")
#tp53vafSim = read_delim(vaffilepath, delim = ",")


## number of samples collected for each age range pulled from Martincorena data
#nSamples <- c(84, 93, 96, 94, 95, 98, 95, 95, 94)
nSamples <- c(84, 93, 96, 94, 95, 98)
sum(nSamples)
## calculate the average number of mutations per tissue sample
tp53_persample <- tp53subsample %>%
  group_by(age) %>%
  summarise(MutationCount = n()) %>%
  mutate(nSamples = nSamples) %>%
  mutate(AverageMutationCount = MutationCount / (nSamples)) %>%
  mutate(blockProb = "IM_data")

## line plot of martincorena data. Mutation frequency per sample by age range.
ggplot(tp53_persample, aes(x = age, y = AverageMutationCount, group = 1)) +
  geom_line() +
  labs(x = "Age Range", y = "Mutation frequency per subsample") +
  theme_minimal()

#######################
##### simulation data

## adjust age ranges to match martincorena data
tp53simAge <- tp53sim %>%
  #rename(age = time) %>% ## toggle on if using python
  mutate(age = case_when(
    age == 21.5 ~ "20-23",
    age == 25.5 ~ "24-27",
    age == 37.5 ~ "36-39",
    age == 45.5 ~ "44-47",
    age == 49.5 ~ "48-51",
    age == 53.5 ~ "52-55",
    ## below age conversion uses schenk default time
    # age == 4.5 ~ "20-23",
    # age == 5.5 ~ "24-27",
    # age == 8.5 ~ "36-39",
    # age == 10 ~ "44-47",
    # age == 11 ~ "48-51",
    # age == 12 ~ "52-55",
    # age == 13 ~ "56-59",
    # age == 15.5 ~ "68-71",
    # age == 16.5 ~ "72-75",
    TRUE ~ as.character(age)  # Keep the original value if none of the conditions match
  )) %>%
  mutate(blockProb = as.character(blockProb))


########
#### mutation rate adjustment
########

## filter the mutation rate
tp53simAge1 <- tp53simAge %>% 
  filter(mutRate == vis_mutRate)

## apply threshold for counting mutations
threshold_vector1 = unique(tp53simAge1$threshold)
df_list <- list()
for (val in threshold_vector1) {
  duplicated_df <- tp53_persample
  duplicated_df$threshold <- val
  df_list <- append(df_list, list(duplicated_df))
}
tp53_persample1 <- do.call(rbind, df_list)

## combine the data between martincorena and simulated tissue
combined_tibble1 <- bind_rows(tp53simAge1, tp53_persample1)
combined_tibble1 <- combined_tibble1 %>%
  mutate(blockProb = as.numeric(blockProb))


## lineplot of average mutation counts per tissue section - martincorena and simulated tissue

combined_tibble_plot1 <- combined_tibble1 %>%
  mutate(blockProb = replace_na(as.character(blockProb), "Empirical"))
palette <- c("0" = "#fcae91",
             "0.007" = "#fb6a4a",
             "0.01" = "#ef3b2c",
             "0.0125" = "#cb181d",
             "0.015" = "#99000d",
             "Empirical" = "#636363")

numPlot1 <- ggplot(combined_tibble_plot1, aes(x = age, y = AverageMutationCount, group = as.factor(blockProb), color = as.factor(blockProb))) +
  geom_point(size = 5) +
  geom_line() +
  scale_color_manual(values = palette) +
  labs(x = "Age Range", y = "Mutations per Sample", color = "p") +
  theme_minimal() +
  #facet_wrap(~threshold, scales = "free", nrow = 2) +
  theme(strip.text = element_text(size = 20),
        legend.position = c(0.15,0.7)
        ) + 
  ylim(c(0,1))


show(numPlot1)
combined_tibble1 %>% print(n=Inf)

########
#### default mutation rate
########

## filter for HomeostaticEpidermis default mutation rate
tp53simAge2 <- tp53simAge %>% 
  filter(mutRate == 1)

## apply threshold for counting mutations
threshold_vector2 = unique(tp53simAge2$threshold)
df_list <- list()
for (val in threshold_vector2) {
  duplicated_df <- tp53_persample
  duplicated_df$threshold <- val
  df_list <- append(df_list, list(duplicated_df))
}
tp53_persample2 <- do.call(rbind, df_list)

## combine the data between martincorena and simulated tissue
combined_tibble2 <- bind_rows(tp53simAge2, tp53_persample2)
combined_tibble2 <- combined_tibble2 %>%
  mutate(blockProb = as.numeric(blockProb))


## lineplot of average mutation counts per tissue section - martincorena and simulated tissue

combined_tibble_plot2 <- combined_tibble2 %>%
  mutate(blockProb = replace_na(as.character(blockProb), "Empirical"))
### undo
#combined_tibble_plot2 <- combined_tibble_plot2 %>% filter(blockProb != 0)
palette <- c("0" = "#fcae91",
             "0.007" = "#fb6a4a",
             "0.01" = "#ef3b2c",
             "0.0125" = "#cb181d",
             "0.015" = "#99000d",
             "Empirical" = "#636363")

numPlot2 <- ggplot(combined_tibble_plot2, aes(x = age, y = AverageMutationCount, group = as.factor(blockProb), color = as.factor(blockProb))) +
  geom_point(size = 5) +
  geom_line() +
  scale_color_manual(values = palette) +
  labs(x = "Age Range", y = "Mutations per Sample", color = "p") +
  theme_minimal() +
  #facet_wrap(~threshold, scales = "free", nrow = 2) +
  theme(strip.text = element_text(size = 20), 
        legend.position = c(0.15,0.7)
        ) + 
  ylim(c(0,1))
show(numPlot2)
combined_tibble2 %>% print(n=Inf)



## save plot to file
if(cluster){
  ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{subsample_path}/dat/martincorena_comparison_thresholds_mutRate_{vis_mutRate}.png"), plot, width = 12, height = 6, units = "in", dpi = 300)
}else{
  ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{subsample_path}/dat/martincorena_comparison_thresholds.png"), plot, width = 12, height = 6, units = "in", dpi = 300)
}


############### Mean Squared Error
## Mean mutation count per sample of IM data and of simulations
## Calculate difference of each mean and square it
## Average across time points

IM_df <- combined_tibble1 %>%
  filter(is.na(blockProb)) %>%
  dplyr::select(age, IM_averageMutCount = AverageMutationCount)

mse_df <- combined_tibble1 %>%
  filter(!is.na(blockProb)) %>%
  left_join(IM_df, by = "age") %>%
  mutate(
    squared_error = (AverageMutationCount - IM_averageMutCount) ^ 2
  ) %>%
  #filter(blockProb != 0) %>%
  group_by(blockProb, age)
  
mse_df

ggplot(mse_df, aes(x = age, y = squared_error, group = blockProb, color = factor(blockProb))) +
  geom_line() +
  labs(x = "Age", y = "MSE", color = "p")

## average across all timepoints
avg_mse <- mse_df %>%
  group_by(blockProb) %>%
  summarize(avg_MSE = mean(squared_error)) %>%
  rename(p = blockProb)
avg_mse

## average across timepoints excluding age 44-47

mse_df_44 <- mse_df %>% filter(age!="44-47")
ggplot(mse_df_44, aes(x = age, y = squared_error, group = blockProb, color = factor(blockProb))) +
  geom_line() +
  labs(x = "Age", y = "MSE", color = "p")

avg_mse <- mse_df_44 %>%
  group_by(blockProb) %>%
  summarize(avg_MSE = mean(squared_error)) %>%
  rename(p = blockProb)
avg_mse

palette <- c("0" = "#fcae91",
             "0.007" = "#fb6a4a",
             "0.01" = "#ef3b2c",
             "0.0125" = "#cb181d",
             "0.015" = "#99000d",
             "Empirical" = "#636363")

MSE <- ggplot(avg_mse, aes(x=factor(p), y=avg_MSE, fill=factor(p))) + 
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = palette) +
  labs(x = "Persistence coefficient (p)", y = "Mean Squared Error") +
  theme_minimal() + 
  geom_text(aes(label = round(avg_MSE, 4)), vjust = -0.3, size = 3.5) +
  theme(legend.position = 'none')
show(MSE)

ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{subsample_path}/dat/martincorena_comparison_thresholds_mse.png"), MSE, width = 6, height = 4, units = "in", dpi = 300)





############################################################################
############################################################################
##### VAF 
############################################################################
############################################################################



## read in files
#tp53sim = read_delim(simfilepath, delim = ",")
tp53subsample = read_delim(subsamplefilepath, delim = ",")
tp53vafSim = read_delim(vaffilepath, delim = ",")

########
### VAF of martincorena SNVs
########

#### plot martincorena vaf values only from subsected samples
# ggplot(tp53subsect, aes(x = age, y = vaf)) +
#   geom_point() + #position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 2) +
#   labs(x = "Age",
#        y = "VAF") +
#   theme_minimal()
# print(tp53subsect)

######### VAF combined simulations and martincorena 

## re-assign age ranges based on martincorena data
tp53vafSim <- tp53vafSim %>%
  mutate(age = case_when(
    age == 21.5 ~ "20-23",
    age == 25.5 ~ "24-27",
    age == 37.5 ~ "36-39",
    age == 45.5 ~ "44-47",
    age == 49.5 ~ "48-51",
    age == 53.5 ~ "52-55",
    # age == 4.5 ~ "20-23",
    # age == 5.5 ~ "24-27",
    # age == 8.5 ~ "36-39",
    # age == 10 ~ "44-47",
    # age == 11 ~ "48-51",
    # age == 12 ~ "52-55",
    # age == 13 ~ "56-59",
    # age == 15.5 ~ "68-71",
    # age == 16.5 ~ "72-75",
    TRUE ~ as.character(age)  # Keep the original value if none of the conditions match
  ))


## only grab needed columns and assign IM label to martincorena data
tp53vaf_subsample = tp53subsample %>%
  dplyr::select(vaf, age, size) %>%
  mutate(blockProb = "Empirical") %>%
  rename(count = size)

## order simulation data from weakest to strongest block and filter LLOD
tp53vafSim <- tp53vafSim %>%
  arrange(desc(blockProb)) %>%
  mutate(blockProb = as.factor(blockProb)) %>%
  rename(vaf = VAF) %>%
  filter(count >= lower_detection_limit)

## combine simulation data and IM data
tp53_allVAF <- bind_rows(tp53vafSim, tp53vaf_subsample)

## filter for mutRate
tp53_allVAF <- tp53_allVAF %>% filter(mutRate == vis_mutRate | is.na(mutRate)) %>%
  #filter(blockProb != 0) %>%
  mutate(blockProb = if_else(blockProb == "Empirical", "Empirical", blockProb))

### plot VAF of simulations and martincorena
# Define color
palette <- c("0" = "#fcae91",
             "0.007" = "#fb6a4a",
             "0.01" = "#ef3b2c",
             "0.0125" = "#cb181d",
             "0.015" = "#99000d",
             "Empirical" = "#636363")
# tp53_allVAF <- tp53_allVAF %>%
#   filter(!(blockProb == 0.0125 | blockProb == 0.015))

## plot
raw_VAF <- ggplot(tp53_allVAF, aes(x = age, y = vaf, color = blockProb)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.01, dodge.width = .9), size = 2) +
  labs(x = "Age",
       y = "VAF",
       color = "TP53 Block Probability") +
  scale_color_manual(values = palette) +
  #breaks = c(NA, factor_levels),  # Set the breaks in numerical order
  #labels = c("IM_data", factor_levels)) +  # Set the labels in numerical order
  theme_minimal() +
  theme(aspect.ratio = 0.5)
show(raw_VAF)

## save plot to file
if(cluster){
  ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{subsample_path}/dat/martincorena_comparison_thresholds_VAF.png"), raw_VAF, width = 10, height = 5, units = "in", dpi = 300)
}else{
  ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{subsample_path}/dat/martincorena_comparison_thresholds_VAF.png"), raw_VAF, width = 10, height = 5, units = "in", dpi = 300)
}


########## plot distributions of mutation VAFs

plot_vaf_distribution <- function(data, p_list, bin_width = bin_width, my_palette) {
  for (p in p_list) {
    # Filter data for the specified blockProb value
    filtered_data <- data[data$blockProb == p | data$blockProb == "Empirical", ]
    print(filtered_data)
    # Create the plot
    plot <- ggplot(filtered_data, aes(x = vaf, fill = blockProb)) +
      geom_histogram(binwidth = bin_width, position = "identity", alpha = 0.5) +
      geom_density(alpha = 0.5) + 
      labs(x = "VAF",
           y = "Count",
           fill = "p",
           title = paste("Distribution of VAF for blockProb =", p, "or blockProb = IM_data")) +
      scale_fill_manual(values = my_palette) +
      theme_minimal() +
      facet_wrap(~ age, ncol = 3) +  # Separate by age with 3 columns
      theme(aspect.ratio = 0.5, legend.position = "top")
    
    
    print(plot)
  }
}


p_list <- c(0.007, 0.01, 0.0125, 0.015)
bin_width <- 0.005
plot_vaf_distribution(tp53_allVAF, p_list, bin_width, palette)

#### fit an exponential distribution to the data

# analyze_vaf <- function(data, p_list, lower_detection_limit, ages){
#   
#   # initialize empty tibble
#   results <- tibble(age = character(), p = numeric(), loglike = numeric(), Sim_count = numeric(), IM_count = numeric())
#   
#   # iterate through each persistence coefficient 
#   for (p in p_list){ 
#     
#     filtered_data_p <- data[data$blockProb == p, ]
#     filtered_data_IM <- data[data$blockProb == "IM data", ]
# 
#     for (i in ages){
#       ## filter simulated data by age
#       filtered_data_p_age <- filtered_data_p %>%
#         filter(age == i)
#       vafs_p <- filtered_data_p_age$vaf
#       
#       ## filter martincorena data by age
#       filtered_data_IM_age <- filtered_data_IM %>%
#         filter(age == i)
#       vafs_IM <- filtered_data_IM_age$vaf
#       
#       ## subtract off limit of detection
#       offset <- lower_detection_limit_VAF
#       
#       ## find the best fit distribution for simulated data
#       best_fit_dist <- fitdist(vafs_p - offset, "exp")
#       
#       print(best_fit_dist$estimate)
#       ## log-likelihood of martincorena data under the fitted exponential
#       loglike = sum(dexp(vafs_IM - offset, best_fit_dist$estimate, log = TRUE))
#       
#       ## how many martincorena values are used to find total log-likelhood
#       IM_count <- length(vafs_IM)
#       Sim_count <- length(vafs_p)
#       
#       # Add resutls to tibble
#       results <- results %>% add_row(age = i, p = p, loglike = loglike, Sim_count = Sim_count, IM_count = IM_count)
#     }
#   }
#   return(results)
# }


analyze_vaf_integral <- function(data, p_list, lower_detection_limit, ages, eps){
  
  # initialize empty tibble
  results <- tibble(age = character(), p = numeric(), loglike = numeric(), Sim_count = numeric(), IM_count = numeric(), lambda = numeric())
  
  # iterate through each persistence coefficient 
  for (p in p_list){ 
    
    filtered_data_p <- data[data$blockProb == p, ]
    filtered_data_IM <- data[data$blockProb == "Empirical", ]
    
    for (i in ages){
      ## filter simulated data by age
      filtered_data_p_age <- filtered_data_p %>%
        filter(age == i)
      vafs_p <- filtered_data_p_age$vaf
      
      ## filter martincorena data by age
      filtered_data_IM_age <- filtered_data_IM %>%
        filter(age == i)
      vafs_IM <- filtered_data_IM_age$vaf
      
      ## subtract off limit of detection
      offset <- lower_detection_limit_VAF
      print(vafs_p)
      ## find the best fit distribution for simulated data
      best_fit_dist <- fitdist(vafs_p - offset, "exp")
      
      print(min(vafs_IM - offset))
      print(best_fit_dist$estimate)
      
      vafs_IM_offset <- vafs_IM - offset
      lower <- pexp(vafs_IM_offset - eps, best_fit_dist$estimate, lower.tail=TRUE)
      upper <- pexp(vafs_IM_offset + eps, best_fit_dist$estimate, lower.tail=TRUE)
      
      integral_density <- upper - lower

      loglike <- -sum(log(integral_density))
      cat('loglike')
      print(loglike)
      
      ## how many martincorena values are used to find total log-likelhood
      IM_count <- length(vafs_IM)
      Sim_count <- length(vafs_p)
      
      # Add resutls to tibble
      results <- results %>% add_row(age = i, p = p, loglike = loglike, Sim_count = Sim_count, IM_count = IM_count, lambda = best_fit_dist$estimate)
    }
  }
  return(results)
}


p_list <- c(0.007, 0.01, 0.0125, 0.015)
ages <- c("20-23", "24-27", "36-39", "44-47", "48-51", "52-55")
bin_width <- 0.005
eps <- 0.0002
lower_detection_limit_VAF <- lower_detection_limit /(15000*2*(70*70/15000))

results <- analyze_vaf_integral(tp53_allVAF, p_list, lower_detection_limit, ages, eps)

options(pillar.sigfig = 5)

methods_chart_filtering <- results %>% 
  filter(age != "20-23") %>%
  dplyr::select(IM_count, lambda, loglike, Sim_count,  lambda, age, p) %>%
  mutate(
    lambda = round(lambda, 2),
    loglike = round(loglike, 2),
  )

write.csv(methods_chart_filtering, glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{subsample_path}/dat/methods_chart.csv"))




results_exclude <- results %>% filter(age != "20-23") %>% group_by(p) %>% summarise(sum_loglike = sum(loglike))
print(results_exclude)
results %>% group_by(p) %>% summarise(sum(loglike))

palette <- c("0" = "#fcae91",
             "0.007" = "#fb6a4a",
             "0.01" = "#ef3b2c",
             "0.0125" = "#cb181d",
             "0.015" = "#99000d",
             "Empirical" = "#636363")

loglike <- ggplot(results_exclude, aes(x=factor(p), y=sum_loglike, color=factor(p))) +
  #geom_bar(stat = "identity", width = 0.7) +
  #scale_fill_manual(values = palette) +
  geom_point(size = 2) +
  scale_color_manual(values = palette) +
  xlab("Persistence coefficient (p)") +
  ylab("-Log-likelihood sum") +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_text(aes(label = round(sum_loglike, 2)), vjust = -0.3, size = 3.5)
show(loglike)





ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{subsample_path}/dat/martincorena_comparison_thresholds_log_like.png"), loglike, width = 6, height = 4, units = "in", dpi = 300)






############################
####### Supplemental figure
############################


library("ggpubr")

PLOT_MARGIN = 1

palette <- c("0" = "#fcae91",
             "0.007" = "#fb6a4a",
             "0.01" = "#ef3b2c",
             "0.0125" = "#cb181d",
             "0.015" = "#99000d",
             "Empirical" = "#636363")

#### plot 1 default mutation rate
numPlot1 <- ggplot(combined_tibble_plot1, aes(x = age, y = AverageMutationCount, group = as.factor(blockProb), color = as.factor(blockProb))) +
  geom_point(size = 1) +
  geom_line() +
  scale_color_manual(values = palette) +
  labs(x = "Age Range", y = "Mutations per sample", color = "p") +
  theme_minimal() +
  #facet_wrap(~threshold, scales = "free", nrow = 2) +
  theme(
    legend.position = c(0.15,0.7),
    legend.key.size = unit(0.5, "lines"),
    legend.key.height = unit(0.5, "lines"),
    legend.margin = margin(0, 0, 0, 0),
    legend.key.spacing.y = unit(0.0001, "cm"),
    legend.key.spacing.x = unit(0.1, "cm"),
    
    strip.text = element_text(size = 8, color='black'),
    axis.text = element_text(size = 8, color='black'),
    axis.title = element_text(size = 8, color='black'), 
    legend.text = element_text(size = 8, color='black', margin = margin(l = 1)), 
    legend.title = element_text(size = 8, color='black', margin = margin(r = 2)),
    axis.line = element_line(color="black"),
    axis.ticks = element_line(color = "black"),
    #plot.margin = margin(0.5, 25.5, 0.5, 5.5), 
    plot.margin = margin(PLOT_MARGIN, 0, PLOT_MARGIN, 0),  
    panel.spacing = unit(1.5, "lines")
  )+
  # theme(strip.text = element_text(size = 20),
  #       legend.position = c(0.15,0.7),
  #       axis.title.x = element_text(size = 8),
  #       axis.title.y = element_text(size = 8),
  #       axis.text.x = element_text(size = 8),
  #       axis.text.y = element_text(size = 8)
  # ) + 
  ylim(c(0,1))


#### plot 2 mutation rate/4
numPlot2 <- ggplot(combined_tibble_plot2, aes(x = age, y = AverageMutationCount, group = as.factor(blockProb), color = as.factor(blockProb))) +
  geom_point(size = 1) +
  geom_line() +
  scale_color_manual(values = palette) +
  labs(x = "Age Range", y = "Mutations per sample", color = "p") +
  theme_minimal() +
  #facet_wrap(~threshold, scales = "free", nrow = 2) +
  theme(
    legend.position = c(0.15,0.7),
    legend.key.size = unit(0.5, "lines"),
    legend.key.height = unit(0.5, "lines"),
    legend.margin = margin(0, 0, 0, 0),
    legend.key.spacing.y = unit(0.0001, "cm"),
    legend.key.spacing.x = unit(0.1, "cm"),
    
    strip.text = element_text(size = 8, color='black'),
    axis.text = element_text(size = 8, color='black'),
    axis.title = element_text(size = 8, color='black'), 
    legend.text = element_text(size = 8, color='black', margin = margin(l = 1)), 
    legend.title = element_text(size = 8, color='black', margin = margin(r = 2)),
    axis.line = element_line(color="black"),
    axis.ticks = element_line(color = "black"),
    #plot.margin = margin(0.5, 25.5, 0.5, 5.5), 
    plot.margin = margin(PLOT_MARGIN, 0, PLOT_MARGIN, 0),  
    panel.spacing = unit(1.5, "lines")
  )+
  # theme(strip.text = element_text(size = 20), 
  #       legend.position = c(0.15,0.7),
  #       axis.title.x = element_text(size = 8),
  #       axis.title.y = element_text(size = 8),
  #       axis.text.x = element_text(size = 8),
  #       axis.text.y = element_text(size = 8)
  # ) + 
  ylim(c(0,1))
show(numPlot2)
combined_tibble2 %>% print(n=Inf)


#### plot 3 raw VAF
raw_VAF <- ggplot(tp53_allVAF, aes(x = age, y = vaf, color = blockProb)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.01, dodge.width = .9), size = 1) +
  labs(x = "Age Range",
       y = "Variant allele frequency",
       color = "p") +
  scale_color_manual(values = palette) +
  #breaks = c(NA, factor_levels),  # Set the breaks in numerical order
  #labels = c("IM_data", factor_levels)) +  # Set the labels in numerical order
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.key.size = unit(0.5, "lines"),
    legend.key.height = unit(0.5, "lines"),
    legend.margin = margin(0, 0, 0, 0),
    legend.key.spacing.y = unit(0.0001, "cm"),
    legend.key.spacing.x = unit(0.1, "cm"),
    strip.text = element_text(size = 8, color='black'),
    axis.text = element_text(size = 8, color='black'),
    axis.title = element_text(size = 8, color='black'), 
    legend.text = element_text(size = 8, color='black', margin = margin(l = 1)), 
    legend.title = element_text(size = 8, color='black', margin = margin(r = 2)),
    axis.line = element_line(color="black"),
    axis.ticks = element_line(color = "black"),
    #plot.margin = margin(0.5, 25.5, 0.5, 5.5), 
    plot.margin = margin(PLOT_MARGIN, 0, PLOT_MARGIN, 0),  
    panel.spacing = unit(1.5, "lines")
  )
  # theme(strip.text = element_text(size = 20), 
  #       #legend.position = c(0.08,0.5),
  #       legend.position = "top",
  #       axis.title.x = element_text(size = 8),
  #       axis.title.y = element_text(size = 8),
  #       axis.text.x = element_text(size = 8),
  #       axis.text.y = element_text(size = 8)
  # ) #+ 
  #theme(aspect.ratio = 0.5)
show(raw_VAF)


#### plot 4 MSE
MSE <- ggplot(avg_mse, aes(x=factor(p), y=avg_MSE, fill=factor(p))) + 
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = palette) +
  labs(x = expression(italic(p)[italic(FANC)^'+' ]^italic(TP53)^'−'), y = "Mean squared error") +
  theme_minimal() + 
  geom_text(aes(label = round(avg_MSE, 4)), vjust = -0.3, size = 8/.pt) +
  theme(
    legend.position = 'none',
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
    #plot.margin = margin(0.5, 25.5, 0.5, 5.5), 
    plot.margin = margin(PLOT_MARGIN, 0, PLOT_MARGIN, 0),  
    panel.spacing = unit(1.5, "lines")
  )+
  # theme(
  #   legend.position = 'none',
  #   axis.title.x = element_text(size = 8),
  #   axis.title.y = element_text(size = 8),
  #   axis.text.x = element_text(size = 8),
  #   axis.text.y = element_text(size = 8),
  #   plot.margin = unit(c(1, 1, 1, 1), "cm")
  # ) +
  ylim(0, max(avg_mse$avg_MSE) * 1.175)


#### plot 5 Log likelihood
# Create the loglike plot
loglike <- ggplot(results_exclude, aes(x=factor(p), y=sum_loglike, color=factor(p))) +
  #geom_bar(stat = "identity", width = 0.7) +
  #scale_fill_manual(values = palette) +
  geom_point(size =2)+ 
  scale_color_manual(values = palette) +
  labs(x = expression(italic(p)[italic(FANC)^'+' ]^italic(TP53)^'−'), y = "-Log-likelihood sum") +
  theme_minimal() +
  geom_text(aes(label = round(sum_loglike, 2)), vjust = -1, size = 8/.pt) +
  theme(
    legend.position = 'none',
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
    #plot.margin = margin(0.5, 25.5, 0.5, 5.5), 
    plot.margin = margin(PLOT_MARGIN, 0, PLOT_MARGIN, 0), 
    panel.spacing = unit(1.5, "lines")
  ) + 
  ylim(550, 580) 
  #ylim(0, max(results_exclude$sum_loglike) * 1.25) 

# Combine the plots into one
#supp_plot <- plot_grid(numPlot2, numPlot1, raw_VAF, MSE, loglike, ncol = 3, rel_heights = c(1, 0.5), labels = "AUTO")

# numPlot2 <- numPlot2 + theme(legend.position = "none")
# numPlot1 <- numPlot1 + theme(legend.position = "none")
# raw_VAF <- raw_VAF + theme(legend.position = "none")

numPlot2 <- numPlot2 + labs(color = expression(italic(p)[italic(FANC)^'+']^italic(TP53)^'−'))
numPlot1 <- numPlot1 + labs(color = expression(italic(p)[italic(FANC)^'+']^italic(TP53)^'−'))
raw_VAF <- raw_VAF + labs(color = expression(italic(p)[italic(FANC)^'+']^italic(TP53)^'−'))

# supp_plot <- plot_grid(
#   numPlot2, numPlot1, raw_VAF,                 # First row
#   NULL, MSE, loglike,                          # Second row, with NULL to align MSE and loglike
#   ncol = 3,                                    # Ensure a 3-column layout
#   rel_heights = c(1, 0.65),
#   align = 'v'                                  # Align vertically
# )


# supp_plot <- ggarrange(
#   numPlot2, numPlot1, raw_VAF,    # First row of plots
#   NULL, MSE, loglike,             # Second row of plots
#   labels = c("A", "B", "C", "", "D", "E"),
#   ncol = 3,                       # Set 3 columns to match your layout
#   nrow = 2,                       # Set 2 rows to match layout
#   common.legend = TRUE,           # Single shared legend
#   heights = c(1.5, 1),
#   legend = "top",
#   font.label = list(size = 14, family = "Arial", face = "bold"),  # Set labels to bold Arial 14
#   label.y = 1.05 
# )
# 
# # Show the final combined plot with the shared legend on top
# supp_plot
# 
# show(supp_plot)


top_row <- ggarrange(
  numPlot2, numPlot1, raw_VAF,
  labels = c("A", "B", "C"),
  ncol = 3,
  font.label = list(size = 14, family = "Arial", face = "bold"),
  label.y = 1.1,
  label.x = -0.025,
  common.legend = TRUE,
  legend = "top" 
)

show(top_row)
# Arrange the bottom row with a different label.y value
bottom_row <- ggarrange(
  NULL, MSE, loglike,
  labels = c("", "D", "E"),
  ncol = 3,
  font.label = list(size = 14, family = "Arial", face = "bold"),
  label.y = 1.175,  # Adjust label position for the bottom row
  label.x = -0.025
)
show(bottom_row)
# Combine the top and bottom rows, and add the shared legend at the top
supp_plot <- ggarrange(
  top_row, bottom_row,
  ncol = 1,
  heights = c(1.6, 1)      # Adjust row heights
)

show(supp_plot)


ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{subsample_path}/dat/martincorena_comparison_p53.png"), supp_plot, width = 7, height = 3.5, units = "in", dpi = 300)
#ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{subsample_path}/dat/martincorena_comparison_p53.png"), supp_plot)

