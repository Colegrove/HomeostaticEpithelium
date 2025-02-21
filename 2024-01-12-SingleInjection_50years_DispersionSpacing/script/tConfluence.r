## FA dynamics
## clone size distribution plot
## 18 Jan 2024

require(tidyverse)
require(ggplot2)
require(dplyr)
require(tidyr)
require(readr)
require(tibble)
require(cowplot)
require(khroma)
require(RColorBrewer)
require(grid)
require(png)
#require(viridis)

## read in data csv as tibble
#outcomes <- tibble(read.csv("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/dat/timeToConfluence.csv"))
outcomes <- tibble(read.csv("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/dat/timeToConfluence.csv"))


########
### Confluence/Loss filtering
########

# proportionSummaryLoss <- outcomes %>% 
#     select(-tConfluence) %>% 
#     group_by(Block_prob, Sigma, Dose) %>% 
#     summarize(pLoss = sum(Confluence == "Loss") /n())

# Calculate tConfluence median
dataSummaries <- outcomes %>%
  group_by(Confluence, Block_prob, Sigma, Dose) %>%
  summarize(Q1 = quantile(tConfluence, 0.25, na.rm = TRUE),
            Q3 = quantile(tConfluence, 0.75, na.rm = TRUE),
            tConfluence = median(tConfluence)
            )

# Filter data for Confluence and Loss separately
dataSummariesConfluence <- dataSummaries %>%
  filter(Confluence == "Confluence")

dataSummariesLoss <- dataSummaries %>%
  filter(Confluence == "Loss")

violet_shades <- c("2" = "#C2A5CF", "10" = "#9970AB", "20" = "#762A83")

LINE_SIZE = 0.5

### Images of basal layer injections
#dose3sigma2 <-
#dose3sigma10 <- readPNG("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/results/basal_sigma_10_dose_3.png")
#dose3sigma20
#dose10sigma2
#dose10sigma10 <- readPNG("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/results/basal_sigma_10_dose_10.png")
#dose10sigma20
#dose30sigma2
#dose30sigma10 <- readPNG("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/results/basal_sigma_10_dose_30.png")
#dose30sigma20

######## 
### Facet by Dose
########

block_levels <- c(0, 0.001)
sigma_levels <- c(2, 10, 20)
dose_levels <- c(3, 10, 30)

#### Null data for plotting style
null_data <- expand.grid(Confluence = "Confluence",
                         Block_prob = block_levels,
                         Sigma = sigma_levels,
                         Dose = dose_levels)
null_data$tConfluence <- NA
dataSummariesConfluence <- bind_rows(dataSummariesConfluence, null_data)

## include only dose = 10
dataSummariesConfluenceDose <- dataSummariesConfluence %>% 
    filter(Dose == 10)
proportionSummaryLossDose <- proportionSummaryLoss %>%
    filter(Dose == 10)

spacing_labels <- c("low", "med", "high")
## Create line plot for tConfluence - facet by Dose
linePlotConfluence <- dataSummariesConfluenceDose %>%
  ggplot(aes(x = as.factor(Block_prob), y = tConfluence, color = as.factor(Sigma), group = interaction(Sigma, Dose))) +
  geom_line(linewidth = LINE_SIZE) +
  geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.2) +
  facet_wrap(~ Dose, scales = "free", ncol = 3, drop = FALSE) +
  theme_classic() +
  labs(x = expression("Persistence coefficient (" * italic(p) * ")"), y = "Years to confluence") +
  scale_color_manual(expression("Transgene diffusion (" * italic(D) * ")"), values = violet_shades, labels=spacing_labels) +
  scale_y_continuous(limits = c(0,50), breaks = seq(0, 50, by = 10)) +
  theme(strip.text = element_blank(), 
    # legend.position = c(0.8, 0.65),
    legend.position = c(0.25, 0.4),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    legend.background = element_blank(), 
    text=element_text(size=8), 
    legend.key.size = unit(0.5, "lines"),
    legend.key.height = unit(0.5, "lines"),
    legend.margin = margin(0, 0, 0, 0),
    legend.key.spacing.y = unit(0.0001, "cm"),
    legend.key.spacing.x = unit(0.0001, "cm"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8, color = 'black'),
    legend.text=element_text(size=8)
    )

confluence_spacing_dose_10 <- linePlotConfluence
show(confluence_spacing_dose_10)

#ggsave("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/results/tConfluence_spacing_dose_10.png", plot = confluence_spacing_dose_10, height = 4, width = 4)

#### presentation style

# linePlotConfluence <- linePlotConfluence + theme(
#   panel.background = element_rect(fill = "transparent",
#                                   colour = NA_character_), # necessary to avoid drawing panel outline
#   panel.grid.major = element_blank(), # get rid of major grid
#   panel.grid.minor = element_blank(), # get rid of minor grid
#   plot.background = element_rect(fill = "transparent",
#                                  colour = NA_character_), # necessary to avoid drawing plot outline
#   legend.background = element_rect(fill = "transparent"),
#   legend.box.background = element_blank(),
#   #legend.box.background = element_rect(color="#DDDDDD",fill = "transparent"),
#   legend.key = element_rect(fill = "transparent"),
#   axis.line = element_line(color="#DDDDDD"), 
#   axis.text = element_text(color="#DDDDDD", size=14),
#   axis.title.x = element_text(color="#DDDDDD", size=14, face="bold"),
#   axis.title.y = element_text(color="#DDDDDD", size=14, face="bold"),
#   axis.ticks = element_line(color = "#DDDDDD"),
#   legend.text = element_text(color="#DDDDDD", size=14, face="bold"), 
#   legend.title = element_text(color="#DDDDDD", size=14, face="bold")
# )

#ggsave("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/results/tConfluence_lines_facetDose.png", plot = linePlotConfluence, height = 4, width = 4)
#ggsave("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/results/tConfluence_lines_facetDose.png", plot = linePlotConfluence, height = 4, width = 4)


########
### Create line plot for pLoss - Facet by Dose
########

linePlotLoss <- proportionSummaryLossDose %>%
  ggplot(aes(x = as.factor(Block_prob), y = pLoss, color = as.factor(Sigma), group = interaction(Sigma, Dose))) +
  geom_line(linewidth = LINE_SIZE) +
  facet_wrap(~ Dose, scales = "free", ncol = 3, drop = FALSE) +
  theme_classic() +
  labs(x = expression("Persistence coefficient (" * italic(p) * ")"), y = "Loss proportion") +
  scale_color_manual(expression("Transgene diffusion (" * italic(D) * ")"), values = violet_shades, labels=spacing_labels) +
  scale_y_continuous(breaks = seq(0, 1, by = .2)) +
  theme(strip.text = element_blank(), 
  # legend.position = c(0.8, 0.65),
  legend.position = c(0.25, 0.4),
  axis.text.x = element_text(angle = 45, hjust = 1),
  axis.title.x = element_blank(),
  legend.background = element_blank(),
  text=element_text(size=8), 
  legend.key.size = unit(0.5, "lines"),
  legend.key.height = unit(0.5, "lines"),
  legend.margin = margin(0, 0, 0, 0),
  legend.key.spacing.y = unit(0.0001, "cm"),
  legend.key.spacing.x = unit(0.0001, "cm"),
  axis.title = element_text(size = 8),
  axis.text = element_text(size = 8, color = 'black'),
  legend.text=element_text(size=8)
  )
loss_spacing_dose_10 <- linePlotLoss
show(loss_spacing_dose_10)
#ggsave("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/results/tLoss_spacing_dose_10.png", plot = loss_spacing_dose_10, height = 4, width = 4)


####
## poster styling
####
# linePlotLoss_poster <- linePlotLoss
# linePlotLoss_poster <- linePlotLoss_poster + theme(
#   axis.line = element_line(color="black"), 
#   axis.text = element_text(color="black", size=28),
#   axis.title.x = element_text(color="black", size=28),
#   axis.title.y = element_text(color="black", size=28),
#   axis.ticks = element_line(color = "black"),
#   legend.text = element_text(color="black", size=28), 
#   legend.title = element_text(color="black", size=28)
# )
# show(linePlotLoss_poster)
#ggsave("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/results/pLoss_lines_facetDose_poster.png", plot = linePlotLoss_poster, height = 5, width = 10)


####
## Presentation styling
####
# linePlotLoss <- linePlotLoss + theme(
#   panel.background = element_rect(fill = "transparent",
#                                   colour = NA_character_), # necessary to avoid drawing panel outline
#   panel.grid.major = element_blank(), # get rid of major grid
#   panel.grid.minor = element_blank(), # get rid of minor grid
#   plot.background = element_rect(fill = "transparent",
#                                  colour = NA_character_), # necessary to avoid drawing plot outline
#   legend.background = element_rect(fill = "transparent"),
#   legend.box.background = element_blank(),
#   #legend.box.background = element_rect(color="#DDDDDD",fill = "transparent"),
#   legend.key = element_rect(fill = "transparent"),
#   axis.line = element_line(color="#DDDDDD"), 
#   axis.text = element_text(color="#DDDDDD", size=14),
#   axis.title.x = element_text(color="#DDDDDD", size=14, face="bold"),
#   axis.title.y = element_text(color="#DDDDDD", size=14, face="bold"),
#   axis.ticks = element_line(color = "#DDDDDD"),
#   legend.text = element_text(color="#DDDDDD", size=14, face="bold"), 
#   legend.title = element_text(color="#DDDDDD", size=14, face="bold")
# )
# show(linePlotLoss)
#ggsave("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/results/pLoss_lines_facetDose.png", plot = linePlotLoss, height = 4, width = 4)
#ggsave("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/results/pLoss_lines_facetDose.png", plot = linePlotLoss, height = 4, width = 4)


red_shades <- c("3" = "#F4A582", "10" = "#D6604D", "30" = "#B2182B")


########
### Facet by Sigma (dispersion)
########

block_levels <- c(0, 0.001)
sigma_levels <- c(2, 10, 20)
dose_levels <- c(3, 10, 30)
    
#### Null data for plotting style
null_data <- expand.grid(Confluence = "Confluence",
                         Block_prob = block_levels,
                         Sigma = sigma_levels,
                         Dose = dose_levels)
null_data$tConfluence <- NA
dataSummariesConfluence <- bind_rows(dataSummariesConfluence, null_data)


## include only sigma = 10
dataSummariesConfluenceSigma <- dataSummariesConfluence %>% 
    filter(Sigma == 10)
proportionSummaryLossSigma <- proportionSummaryLoss %>%
    filter(Sigma == 10)

# Create line plot for tConfluence - Facet by Sigma
linePlotConfluence <- dataSummariesConfluenceSigma %>%
  ggplot(aes(x = as.factor(Block_prob), y = tConfluence, color = as.factor(Dose), group = interaction(Sigma, Dose))) +
  geom_line(linewidth = LINE_SIZE) +
  geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.2) +
  facet_wrap(~ Sigma, scales = "free", ncol = 3, drop = FALSE) +
  theme_classic() +
  labs(x = expression("Persistence coefficient (" * italic(p) * ")"), y = "Years to confluence") +
  scale_color_manual(expression("# of corrected cells (" * italic(k) * ")"), values = red_shades) +
  scale_x_discrete(drop=FALSE) +  # Set x-axis ticks manually
  scale_y_continuous(limits = c(0,50), breaks = seq(0, 50, by = 10)) +
  theme(strip.text = element_blank(), 
    legend.position = c(0.82, 0.65),
    #legend.title.align = 1,
    legend.background=element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    text=element_text(size=8), 
    legend.title.align=0.5,
    legend.key.size = unit(0.5, "lines"),
    legend.key.height = unit(0.5, "lines"),
    legend.margin = margin(0, 0, 0, 0),
    legend.key.spacing.y = unit(0.0001, "cm"),
    legend.key.spacing.x = unit(0.0001, "cm"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8, color = 'black'),
    legend.text=element_text(size=8)
    )

confluence_dose_spacing_10 <- linePlotConfluence
show(confluence_dose_spacing_10)
#ggsave("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/results/confluence_dose_spacing_10.png", plot = confluence_dose_spacing_10, height = 4, width = 4)


########
### Presentation styling
########
# linePlotConfluence <- linePlotConfluence + theme(
#   panel.background = element_rect(fill = "transparent",
#                                   colour = NA_character_), # necessary to avoid drawing panel outline
#   panel.grid.major = element_blank(), # get rid of major grid
#   panel.grid.minor = element_blank(), # get rid of minor grid
#   plot.background = element_rect(fill = "transparent",
#                                  colour = NA_character_), # necessary to avoid drawing plot outline
#   #legend.background = element_rect(fill = "transparent"),
#   #legend.box.background = element_blank(),
#   #legend.box.background = element_rect(color="#DDDDDD",fill = "transparent"),
#   legend.key = element_rect(fill = "transparent"),
#   axis.line = element_line(color="#DDDDDD"), 
#   axis.text = element_text(color="#DDDDDD", size=14),
#   axis.title.x = element_text(color="#DDDDDD", size=14, face="bold"),
#   axis.title.y = element_text(color="#DDDDDD", size=14, face="bold"),
#   axis.ticks = element_line(color = "#DDDDDD"),
#   #legend.text = element_text(color="#DDDDDD", size=14, face="bold"), 
#   #legend.title = element_text(color="#DDDDDD", size=14, face="bold")
# )
# show(linePlotConfluence)
#ggsave("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/results/tConfluence_multidose_midSigma.png", plot = linePlotConfluence, height = 4, width = 4)
#ggsave("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/results/tConfluence_multidose_midSigma.png", plot = linePlotConfluence, height = 4, width = 6)

########
### Loss
########
proportionSummaryLoss
proportionSummaryLossSigma10 <- proportionSummaryLoss %>%
  filter(Sigma == 10) 
  #filter(Dose == 10 | Dose == 30)
proportionSummaryLossSigma10

#### calculate the IQR

proportionSummaryLossSigma10 %>% print(n=Inf)
# Create line plot for pLoss - Facet by Simga
linePlotLoss <- proportionSummaryLossSigma10 %>%
  ggplot(aes(x = as.factor(Block_prob), y = pLoss, color = as.factor(Dose), group = interaction(Sigma, Dose))) +
  geom_line(linewidth = LINE_SIZE) +
  facet_wrap(~ Sigma, scales = "free", ncol = 3, drop = FALSE) +
  theme_classic() +
  labs(x = expression("Persistence coefficient (" * italic(p) * ")"), y = "Loss proportion") +
  scale_color_manual(expression("# of corrected cells (" * italic(k) * ")"), values = red_shades) +
  scale_y_continuous(breaks = seq(0, 1, by = .2)) +
  theme(strip.text = element_blank(), 
    legend.position = c(0.82, 0.65),
    #legend.title.align = 1,
    legend.background=element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    text=element_text(size=8), 
    legend.title.align=0.5,
    legend.key.size = unit(0.5, "lines"),
    legend.key.height = unit(0.5, "lines"),
    legend.margin = margin(0, 0, 0, 0),
    legend.key.spacing.y = unit(0.0001, "cm"),
    legend.key.spacing.x = unit(0.0001, "cm"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8, color = 'black'),
    legend.text=element_text(size=8)
    )
show(linePlotLoss)
loss_dose_spacing_10 <- linePlotLoss
#ggsave("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/results/loss_dose_spacing_10.png", plot = loss_dose_spacing_10, height = 4, width = 4)




########
### poster styling
########
# linePlotLoss_poster <- linePlotLoss
# linePlotLoss_poster <- linePlotLoss_poster + theme(
#   axis.line = element_line(color="black"), 
#   axis.text = element_text(color="black", size=28),
#   axis.title.x = element_text(color="black", size=28),
#   axis.title.y = element_text(color="black", size=28),
#   axis.ticks = element_line(color = "black"),
#   legend.text = element_text(color="black", size=28), 
#   legend.title = element_text(color="black", size=28)
# )
#show(linePlotLoss_poster)
#ggsave("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/results/pLoss_lines_facetSigma_poster.png", plot = linePlotLoss_poster, height = 5, width = 10)

########
### Research reports styling
########
# linePlotLoss <- linePlotLoss + theme(
#   panel.background = element_rect(fill = "transparent",
#                                   colour = NA_character_), # necessary to avoid drawing panel outline
#   panel.grid.major = element_blank(), # get rid of major grid
#   panel.grid.minor = element_blank(), # get rid of minor grid
#   plot.background = element_rect(fill = "transparent",
#                                  colour = NA_character_), # necessary to avoid drawing plot outline
#   legend.background = element_rect(fill = "transparent"),
#   legend.box.background = element_blank(),
#   #legend.box.background = element_rect(color="#DDDDDD",fill = "transparent"),
#   legend.key = element_rect(fill = "transparent"),
#   axis.line = element_line(color="#DDDDDD"), 
#   axis.text = element_text(color="#DDDDDD", size=14),
#   axis.title.x = element_text(color="#DDDDDD", size=14, face="bold"),
#   axis.title.y = element_text(color="#DDDDDD", size=14, face="bold"),
#   axis.ticks = element_line(color = "#DDDDDD"),
#   legend.text = element_text(color="#DDDDDD", size=14, face="bold"), 
#   legend.title = element_text(color="#DDDDDD", size=14, face="bold")
# )
###########
#show(linePlotLoss)

#ggsave("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/results/pLoss_multiDose_midSigma.png", plot = linePlotLoss, height = 4, width = 6)
#ggsave("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/results/pLoss_multiDose_midSigma_dose10_30.png", plot = linePlotLoss, height = 4, width = 6)



###########################
##### Compile figure panels
###########################

loss_dose_spacing_10 <- loss_dose_spacing_10 + 
  theme(legend.position = "top", 
  legend.key.spacing.x = unit(0.15, "cm")) +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5))

confluence_dose_spacing_10 <- confluence_dose_spacing_10 + 
  theme(legend.position = "top", 
  legend.key.spacing.x = unit(0.15, "cm")) +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5))

loss_spacing_dose_10 <- loss_spacing_dose_10 + 
  theme(legend.position = "top", 
  legend.key.spacing.x = unit(0.15, "cm")) +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5))

confluence_spacing_dose_10 <- confluence_spacing_dose_10 + 
  theme(legend.position = "top", 
  legend.key.spacing.x = unit(0.15, "cm")) + 
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5))


shared_x_label <- textGrob(expression("Persistence coefficient (" * italic(p)[corr] * ")"), gp = gpar(fontsize = 8), vjust = -.2)
top_plots <- plot_grid(
  loss_dose_spacing_10, confluence_dose_spacing_10,
  ncol = 2, labels = c("B", "C"), label_size = 12
)

# Combine bottom plots
bottom_plots <- plot_grid(
  loss_spacing_dose_10, confluence_spacing_dose_10,
  ncol = 2, labels = c("E", "F"), label_size = 12
)

# Add shared x-axis labels using `plot_grid()`
top_with_label <- plot_grid(
  top_plots, 
  ggdraw() + draw_grob(shared_x_label), 
  ncol = 1, rel_heights = c(1, 0.1)
)

bottom_with_label <- plot_grid(
  bottom_plots, 
  ggdraw() + draw_grob(shared_x_label), 
  ncol = 1, rel_heights = c(1, 0.1)
)

# Combine the top and bottom sections
f2_panB_C_D_E <- plot_grid(
  top_with_label, bottom_with_label,
  ncol = 1, rel_heights = c(1, 1)
)

# f2_panB_C_D_E <- plot_grid(loss_dose_spacing_10, confluence_dose_spacing_10, 
#                            loss_spacing_dose_10, confluence_spacing_dose_10, 
#                            ncol = 2, labels = c("B", "C", "E", "F"), 
#                            label_size = 12, rel_heights = c(1,1), 
#                            label_x = c(0, 0, 0, 0)
#                            )
f2_panB_C_D_E

ggsave("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/results/f2_panB_C_D_E.png", plot = f2_panB_C_D_E, height = 3.245, width = 3.25)
