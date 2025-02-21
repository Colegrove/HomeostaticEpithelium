path = "/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-11-18-divisionRateComparison/dat/divFile.txt"

custom_columns <- c("x", "y", "t", "division")
file <- read_delim(path, delim = "\t", col_names = custom_columns)

average_divisions_per_position <- file %>%
  group_by(x, y) %>%                           # Group by position (x, y)
  summarise(total_divisions = sum(division),   # Sum divisions for each position across all timepoints
            .groups = "drop") %>%
  summarise(avg_divisions_per_position = mean(total_divisions)) # Compute the average across all positions


ave_divisions <- round(average_divisions_per_position/100*7, 1)


# Ave divisions per week per position
divisions_per_week <- file %>%
  group_by(x, y) %>% 
  summarise(total_divisions = sum(division), .groups = "drop") %>%
  mutate(divisions_per_week = total_divisions / 100 * 7)


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
    axis.text.x = element_text(angle = 00, hjust = 1, vjust = 1),
    plot.margin = margin(0.5, 2, 0.5, 5.5)
  )


divisions_per_week
plotdist <- ggplot(divisions_per_week, aes(x = divisions_per_week)) +
  geom_histogram(binwidth = 0.13, color = "#DDDDDD", fill = "black") +
  geom_vline(xintercept = 0.4, color = "red", linetype = "solid", size = 1) +
  annotate("text", x = 0.4, 
           y = max(hist(divisions_per_week$divisions_per_week, plot = FALSE)$counts) * 0.95, 
           label = "0.4", color = "red", hjust = -0.2, vjust = -0.6) +
  scale_y_continuous(
    limits = c(0, 600),            # Set y-axis limits
    breaks = seq(0, 600, by = 200) # Set breaks every 200
  ) +
  labs(
    x = "Divisions per week",
    y = "Count"
  ) +
  text8pt_theme
show(plotdist)
ggsave("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-11-18-divisionRateComparison/results/divisionRate_distribution.png", plotdist, width = 1.75, height = 1.75)
