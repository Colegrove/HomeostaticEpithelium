## Hunter Colegrove
## FA dynamics
## 29 Mar 2024

## Speed calculation of the spread of corrected clones.
## Uses corrCellCounts.csv file generated from corrected_cellCounts.R
## outputs averageAreaSpread.csv to plot all figure panels together

require(glue)
require(foreach)
require(tidyverse)
require(khroma)
require(utils)
library(broom)

## Project information and parameters used in experiment (enter from snakefile)
project_path = "2024-04-10-SingleInjection_Confluence_blockProb"
sigmas = c(2)
doses = c(10)
block_values = c(0, 0.001, 0.01, 0.1, 0.2, 0.5, 1.0)
#block_values = c(0, 1.0)
#p53block = 0.005
corrTimes = c(1)
nreplicates = 100
replicate = 0:nreplicates


#write.csv(alldat, glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/corrCellCounts.csv"), row.names = FALSE)

corrCellCounts <- read_csv(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/corrCellCounts.csv"))
#corrCellCounts <- read_csv(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/corrCellCounts.csv"))

cell_density = 15000 # cells/mm^2
area_per_cell_sq_mm <- 1 / cell_density
#area_per_cell_sq_mm <- 10*10 / 1e6  # Convert 10 square micrometers to square millimeters
#total_area_sq_mm <- number_of_cells * area_per_cell_sq_mm

########### cell count numbers convert to area

round_to_nearest <- function(x, step){
  print(x)
  print(step)
  print(round(x/step)*step)
  return(round(x/step)*step)
}

countOnly <- corrCellCounts %>%
  mutate(year = round_to_nearest(year, 0.5)) %>%
  group_by(replicate, corrBlock) %>%
  mutate(maxCount = max(cellCount)) %>%
  filter(maxCount >= 8000) %>% ## filter out non-confluence
  mutate(cellChange = cellCount - lag(cellCount)) %>%
  mutate(change_in_time = year - lag(year)) %>%
  mutate(speed_sq_mm_year = cellChange * area_per_cell_sq_mm / (change_in_time)) %>%
  filter(!is.na(speed_sq_mm_year)) %>%
  ungroup() %>%
  filter(change_in_time > 0)
countOnly %>% print()


## average across replicates
average_countOnly <- countOnly %>%
  group_by(corrBlock, year) %>%
  summarize(
    avg_cellCount = mean(cellCount, na.rm = TRUE),
    avg_speed = mean(speed_sq_mm_year, na.rm = TRUE)
  )
print(average_countOnly, n = 500)

average_countOnly <- average_countOnly %>%
  mutate(avg_tissueArea = avg_cellCount * area_per_cell_sq_mm)
average_countOnly

## output file to plot all figure panels together
write.csv(average_countOnly, glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/averageAreaSpread.csv"), row.names = FALSE)

#### single plot
muted <- color("muted")(9)[c(1, 3, 5:9)] # skip green, keep teal
names(muted) <- c("0.01", "0.1", "0.2", "0.5", "1")


ave_countOnly <- ggplot(average_countOnly, aes(x = year, y = avg_tissueArea, color = as.factor(corrBlock))) +
  geom_line(size = 1.6) +
  scale_color_manual("Blocking Probability (b)", values = muted) +
  labs(x = "Year", y = expression(paste('Area (mm' ^2*')' )), color = "CorrBlock") +
  theme_minimal() +
  theme(legend.position = "none")  # Remove legend

ave_countOnly <- ave_countOnly + theme(
  panel.background = element_rect(fill = "transparent",
                                  colour = NA_character_), # necessary to avoid drawing panel outline
  panel.grid.major = element_blank(), # get rid of major grid
  panel.grid.minor = element_blank(), # get rid of minor grid
  plot.background = element_rect(fill = "transparent",
                                 colour = NA_character_), # necessary to avoid drawing plot outline
  legend.background = element_rect(fill = "transparent"),
  legend.box.background = element_blank(),
  #legend.box.background = element_rect(color="#DDDDDD",fill = "transparent"),
  legend.key = element_rect(fill = "transparent"),
  axis.line = element_line(color="#DDDDDD"), 
  axis.text = element_text(color="#DDDDDD", size=16, face="bold"),
  axis.title.x = element_text(color="#DDDDDD", size=18),
  axis.title.y = element_text(color="#DDDDDD", size=18, face="bold"),
  axis.ticks = element_line(color = "#DDDDDD"),
  legend.text = element_text(color="#DDDDDD", size=16, face="bold"), 
  legend.title = element_text(color="#DDDDDD", size=16, face="bold")
)
show(ave_countOnly)
ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/areaSpread_time.png"), plot = ave_countOnly, height = 4, width = 6)


#################
#### Calculate average speed of each persistence coefficient
#################

### averages of replicates

average_countOnly
fit_line <- average_countOnly %>%
  group_by(corrBlock) %>%
  do(model = lm(avg_tissueArea ~ year, data = .))  # Fit linear model
fitted <- fit_line$model[[1]]
fitted$coefficients[2]

ggplot(average_countOnly, aes(x = year, y = avg_tissueArea)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, aes(color = corrBlock), show.legend = TRUE) +
  facet_wrap(~ corrBlock) +
  labs(
    x = "Year",
    y = "Avg Area (sq mm)"
  ) +
  theme_minimal()

slope_by_block <- average_countOnly %>%
  group_by(corrBlock) %>%
  do(tidy(lm(avg_tissueArea ~ year, data = .))) %>%  # Fit model and tidy output
  filter(term == "year")
  
print(slope_by_block)


### all replicates

countOnly_plot <- countOnly %>%
  mutate(tissue_area = cellCount * area_per_cell_sq_mm)

slope_by_block <- countOnly_plot %>%
  group_by(corrBlock) %>%
  do(tidy(lm(tissue_area ~ year, data = .))) %>%
  filter(term == "year") #%>%
  #select(corrBlock, slope = estimate, std_error = std.error)

print(slope_by_block)

# Plot individual linear regressions for each corrBlock
ggplot(countOnly_plot, aes(x = year, y = tissue_area)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, aes(color = corrBlock), show.legend = TRUE) +
  facet_wrap(~ corrBlock) +
  labs(
    x = "Year",
    y = "Area (sq mm)"
  ) +
  theme_minimal()




############################################################
########## Spacing vs time to convergence of adjacent clones
############################################################
average_radius <- average_countOnly %>%
  mutate(avg_radius_mm = sqrt(avg_tissueArea / pi)) %>%
  mutate(avg_radius_um = avg_radius_mm * 1000) %>%
  mutate(spacing_interaction = avg_radius_um *2) %>%
  print(n=500)

average_radius <- average_radius %>%
  group_by(corrBlock) %>%
  mutate(min_year_8000 = min(ifelse(avg_cellCount >= 8000, year, NA), na.rm = TRUE)) %>%
  filter(year <= min_year_8000 | avg_cellCount <8000) %>%
  select(-min_year_8000)

print(average_radius, n=200)

ave_radius <- ggplot(average_radius, aes(x = spacing_interaction, y = year, color = as.factor(corrBlock))) +
  geom_line() +
  scale_color_manual("Blocking Probability (b)", values = muted) +
  labs(x = expression(paste("Microneedle spacing (", mu, "m)")), y = "Time of clone intersection (year)", color = "CorrBlock") +
  xlim(200, 1010) +
  scale_y_continuous(breaks = seq(0, 55, by = 5)) +
  scale_x_continuous(breaks = seq(0, 1000, by = 100)) +
  theme_minimal() + 
  theme(legend.position = "none")  # Remove legend


ave_radius <- ave_radius + theme(
  panel.background = element_rect(fill = "transparent",
                                  colour = NA_character_), # necessary to avoid drawing panel outline
  panel.grid.major = element_blank(), # get rid of major grid
  panel.grid.minor = element_blank(), # get rid of minor grid
  plot.background = element_rect(fill = "transparent",
                                 colour = NA_character_), # necessary to avoid drawing plot outline
  legend.background = element_rect(fill = "transparent"),
  legend.box.background = element_blank(),
  #legend.box.background = element_rect(color="#DDDDDD",fill = "transparent"),
  legend.key = element_rect(fill = "transparent"),
  axis.line = element_line(color="#DDDDDD"), 
  axis.text = element_text(color="#DDDDDD", size=13),
  axis.title.x = element_text(color="#DDDDDD", size=18),
  axis.title.y = element_text(color="#DDDDDD", size=18),
  axis.ticks = element_line(color = "#DDDDDD"),
  legend.text = element_text(color="#DDDDDD", size=16, face="bold"), 
  legend.title = element_text(color="#DDDDDD", size=16, face="bold")
)

show(ave_radius)
ggsave(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/spacing_intersection_time.png"), plot = ave_radius, height = 5, width = 6)


########### convert diameter to area
remove_loss <- alldat %>%
  group_by(corrBlock, replicate) %>%
  filter(!any(horizontal_diameter == -Inf | vertical_diameter == -Inf | ave_diameter == -Inf)) %>%
  ungroup()
print(remove_loss)

diameter_change <- remove_loss %>%
  group_by(replicate) %>%
  mutate(change_in_diameter = if_else(year == 2/364/4.5, NA_real_, ave_diameter - lag(ave_diameter)),
         change_in_cells = if_else(year == 2/364/4.5, NA_real_, (pi*(ave_diameter/2)^2) - (pi*(lag(ave_diameter)/2)^2)), 
         change_in_time = if_else(year == 2/364/4.5, NA_real_, year - lag(year)),
         speed_sq_mm_year = change_in_cells * area_per_cell_sq_mm / (change_in_time)
        )

print(diameter_change, n=Inf)

# Drop the first row within each replicate as it doesn't have a previous time point for comparison
diameter_change <- diameter_change %>% filter(!is.na(change_in_diameter))
print(diameter_change)

diameter_change %>% filter(corrBlock == 0.01) %>% print(n=Inf)

average_diameter_change <- diameter_change %>%
  group_by(corrBlock, year) %>%
  summarize(
    avg_change_in_diameter = mean(change_in_diameter, na.rm = TRUE),
    avg_change_in_cells = mean(change_in_cells, na.rm = TRUE),
    avg_change_in_time = mean(change_in_time, na.rm = TRUE),
    avg_speed = mean(speed_sq_mm_year, na.rm = TRUE)
  )
print(average_diameter_change, n = Inf)

pi*(35/2)^2
pi*(50/2)^2

write.csv(average_diameter_change, glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/average_diameter_change_speed.csv"), row.names = FALSE)
average_diameter_change <- read_csv(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/average_diameter_change_speed.csv"))

#write.csv(average_diameter_change, glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/average_diameter_change_speed.csv"), row.names = FALSE)
#average_diameter_change <- read_csv(glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/average_diameter_change_speed.csv"))


ave_diameter_no_0.001 <- average_diameter_change %>% filter(!(corrBlock <= 0.001))
ave_diameter_no_0.001
ggplot(ave_diameter_no_0.001, aes(x = year, y = avg_speed, color = as.factor(corrBlock))) +
  geom_line() +
  labs(x = "Year", y = "mm^2/year", color = "CorrBlock") +
  theme_minimal()



# tempdat <- alldat %>%
#   filter(year == 2.5) %>%
#   filter(replicate == 0)
#   #filter(injSite == 0) 
# print(tempdat)
# 
# plotCorrectionEvents <- function(df){
#   #ggplot(df, aes(x = x, y = z, fill = factor(injSite))) +
#   ggplot(df, aes(x = x, y = z, fill = factor(mut))) +
#     geom_tile(width = 1, height = 1) +
#     #scale_fill_manual(values = c("0" = "#44AA99", "-1" = "#DDDDDD")) +
#     scale_fill_manual(values = c("666" = "#44AA99", "0" = "#DDDDDD", "68" = "red")) +
#     labs(x = NULL, y = NULL, fill = "Corrected") +
#     #labs(x = "X", y = "Y", fill = "Corrected") +
#     geom_hline(yintercept = seq(-0.5,99.5), color = "black", size = 0.5) +
#     geom_vline(xintercept = seq(-0.5,99.5), color = "black", size = 0.5) +
#     theme(legend.position="none", 
#           axis.text.x = element_blank(),
#           axis.ticks.x = element_blank(),
#           axis.text.y = element_blank(),
#           axis.ticks.y = element_blank(),
#           plot.margin = grid::unit(c(0,0,-1,-1), "mm"), 
#           plot.background=element_rect(fill="black", colour=NA)) +
#     #coord_cartesian(xlim = c(36, 58), ylim = c(36, 58)) # Setting plot limits
#     coord_cartesian(xlim = c(4, 95), ylim = c(4, 95)) # Setting plot limits
#   #geom_segment(data = segment_lines, aes(x = x, y = y, xend = xend, yend = yend))
# }
# plotCorrectionEvents(tempdat)