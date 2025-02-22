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

## project path of subsampled reconstructed tissue
subsample_path = "2024-06-03-martincorena_spatial_subsample"

## project path of simulations
project_path = "2024-05-03-tp53_martincorena_mutRate_54years"

if(cluster){
  ## old subsect way
  #subsectfilepath = glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/subsect_clones_tp53.csv")
  ## new subsample way
  subsamplefilepath = glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{subsample_path}/dat/subsamples/mutations_subsamples.csv")
  ## tp53 unique clones across detection thresholds
  simfilepath = glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/countClones_tp53_R.csv")
  ## VAF path
  vaffilepath = glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/VAF_tp53_R.csv")
}else{
  ## old subsect way
  #subsectfilepath = glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/subsect_clones_tp53.csv")
  ## new subsample way
  subsamplefilepath = glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{subsample_path}/dat/subsamples/mutations_subsamples.csv")
  ## tp53 unique clones across detection thresholds
  simfilepath = glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/countClones_tp53_R.csv")
  ## VAF path
  vaffilepath = glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/VAF_tp53_R.csv")
}


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
  mutate(blockProb = "IM_data") %>%
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
tp53_allVAF <- tp53_allVAF %>% filter(mutRate == vis_mutRate | is.na(mutRate))

### plot VAF of simulations and martincorena
factor_levels <- c("0.007", "0.01", "0.0125", "0.015")

# Define custom color palette using a colorblind-friendly palette
my_palette <- c("IM_data" = "black",
                setNames(c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"), factor_levels))

# tp53_allVAF <- tp53_allVAF %>%
#   filter(!(blockProb == 0.0125 | blockProb == 0.015))

## plot
plot <- ggplot(tp53_allVAF, aes(x = age, y = vaf, color = blockProb)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.01, dodge.width = .9), size = 2) +
  labs(x = "Age",
       y = "VAF",
       color = "TP53 Block Probability") +
  scale_color_manual(values = my_palette) +
                     #breaks = c(NA, factor_levels),  # Set the breaks in numerical order
                     #labels = c("IM_data", factor_levels)) +  # Set the labels in numerical order
  theme_minimal() +
  theme(aspect.ratio = 0.5)
show(plot)

## save plot to file
if(cluster){
  ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{subsample_path}/dat/martincorena_comparison_thresholds_VAF.png"), plot, width = 10, height = 5, units = "in", dpi = 300)
}else{
  ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{subsample_path}/dat/martincorena_comparison_thresholds_VAF.png"), plot, width = 10, height = 5, units = "in", dpi = 300)
}


########## plot distributions of mutation VAFs

plot_vaf_distribution <- function(data, p_list, bin_width = bin_width, my_palette) {
  for (p in p_list) {
    # Filter data for the specified blockProb value
    filtered_data <- data[data$blockProb == p | data$blockProb == "IM_data", ]
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
plot_vaf_distribution(tp53_allVAF, p_list, bin_width, my_palette)




#### fit an exponential distribution to the data



analyze_vaf <- function(data, p_list, lower_detection_limit, ages, bin_width){
  
  results <- tibble(age = character(), p = numeric(), loglike = numeric(), Sim_count = numeric(), IM_count = numeric())
  
  for (p in p_list){ 
    filtered_data_p <- data[data$blockProb == p, ]
    filtered_data_IM <- data[data$blockProb == "IM_data", ]
    for (i in ages){
      ## timepoint VAFs from simulation
      filtered_data_p_age <- filtered_data_p %>%
        filter(age == i)
      vafs_p <- filtered_data_p_age$vaf
      
      ## timepoint VAFs from martincorena
      filtered_data_IM_age <- filtered_data_IM %>%
        filter(age == i)
      vafs_IM <- filtered_data_IM_age$vaf
      
      ## subtract off whatever is limit of detection
      offset <- lower_detection_limit_VAF
      
      ## find the best fit distribution for simulated data
      best_fit_dist <- fitdist(vafs_p - offset, "exp")
      
      ## log-likelihood of martincorena data under the fitted exponential
      loglike = sum(dexp(vafs_IM, best_fit_dist$estimate, log = TRUE))
      
      ## how many martincorena values are used to find total log-likelhood
      IM_count <- length(vafs_IM)
      Sim_count <- length(vafs_p)
      
      ###### KL divergence
      # min_vaf <- min(vafs_p, vafs_IM)
      # max_vaf <- max(vafs_p, vafs_IM)
      # breaks <- seq(min_vaf, max_vaf + bin_width, by = bin_width)
      # print(breaks)
      # 
      # hist_p <- hist(vafs_p, breaks = breaks, plot = TRUE)
      # hist_IM <- hist(vafs_IM, breaks = breaks, plot = FALSE)
      # prob_p <- hist_p$counts / sum(hist_p$counts)
      # prob_IM <- hist_IM$counts / sum(hist_IM$counts)
      # 
      # kl_divergence <- KL.plugin(prob_IM, prob_p, unit = "log")
      
      
      # Add resutls to tibble
      results <- results %>% add_row(age = i, p = p, loglike = loglike, Sim_count = Sim_count, IM_count = IM_count)
    }
  }
  return(results)
}



p_list <- c(0.007, 0.01, 0.0125, 0.015)
ages <- c("20-23", "24-27", "36-39", "44-47", "48-51", "52-55")
bin_width <- 0.005
lower_detection_limit_VAF <- lower_detection_limit /(15000*2*(70*70/15000))

results <- analyze_vaf(tp53_allVAF, p_list, lower_detection_limit, ages, bin_width)
print(results, n=Inf)


results_exclude <- results %>% filter(age != "20-23") %>% group_by(p) %>% summarise(sum_loglike = sum(loglike))
print(results_exclude)
results %>% group_by(p) %>% summarise(sum(loglike))

palette <- c("0.0" = "#fcae91",
             "0.007" = "#fb6a4a",
             "0.01" = "#ef3b2c",
             "0.0125" = "#cb181d",
             "0.015" = "#99000d",
             "IM data" = "#636363")

loglike <- ggplot(results_exclude, aes(x=factor(p), y=sum_loglike, fill=factor(p))) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = palette) +
  xlab("Persistence coefficient (p)") +
  ylab("Log-likelihood sum") +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_text(aes(label = round(sum_loglike, 2)), vjust = -0.3, size = 3.5)

ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{subsample_path}/dat/martincorena_comparison_thresholds_log_like.png"), loglike, width = 6, height = 4, units = "in", dpi = 300)






GeneLengths <- c(3191618082, 2463, 3666, 1617, 1443, 1444, 13692, 6205, 5506, 1212,
                 4569, 2301, 1440, 888, 2649, 514, 3279, 7329, 2307, 5767, 3633,
                 2931, 3695, 4023, 3927, 2239, 13767, 14944, 2124, 2616, 2614, 
                 2443, 7176, 4391, 2781, 633, 1404, 1976, 1875, 687, 4185, 
                 13482, 8520, 1771, 7668, 7416, 6966, 570, 7351, 1299, 3207, 
                 3750, 3369, 4821, 4337, 1212, 4383, 2787, 2951, 3985, 5376, 
                 6026, 7695, 1659, 2364, 954, 5103, 1455, 1289, 2041, 7098)

mutation_rate <- 3.2e-9
# Calculate expected mutations for each gene
ExpectedMuts <- GeneLengths * mutation_rate

# Create a data frame to display the results
gene_data <- data.frame(
  Gene = 0:(length(GeneLengths) - 1),
  GeneLength = GeneLengths,
  ExpectedMutations = ExpectedMuts
)

# Print the results
print(gene_data)