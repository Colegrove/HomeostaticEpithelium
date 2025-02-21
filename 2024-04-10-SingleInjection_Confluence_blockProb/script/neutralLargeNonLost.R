## Hunter Colegrove
## 26 Sep 2024
## Find the largest non-lost clone size of a neutral clone


dir_path <- "/Volumes/feder-vol1/project/HomeostaticEpithelium/dat/2024-04-10-SingleInjection_Confluence_blockProb/results"
# list files
file_list <- list.files(path = dir_path)

# Filter files to keep only those that start with 'VisFile_block_0_'
filtered_list <- grep("^VisFile_block_0_", file_list, value = TRUE)

# Extract numeric values from the filenames for sorting
numbers <- as.numeric(gsub(".*\\.(\\d+)\\.txt$", "\\1", filtered_list))

# Find the highest integer value
max_value <- max(numbers)

# Find all files that correspond to this max_value
max_files <- filtered_list[numbers == max_value]

# Filter all timepoints of the non-lost replicate (adjusting for your non-lost condition)
replicate_nonLost <- grep("_(45|93|68)\\.[0-9]+\\.txt$", filtered_list, value = TRUE)

# Initialize a dataframe to hold the results
corrected_cell_counts <- data.frame(file = character(), count = integer(), stringsAsFactors = FALSE)

# Loop through each file in the filtered list and process the data
for (file in replicate_nonLost) {
  file_path <- file.path(dir_path, file)  
  column_names = c('x', 'z', 'y', 'h', 's', 'v', 'alpha', 'mut', 'injSite', 'tp53Clone')
  
  # Read the file
  data <- read_tsv(file_path, col_names = column_names)
  
  # Filter and count the cells
  df <- data %>%
    filter(y == 0) %>%
    filter(injSite != -1) %>%
    summarise(count = n())
  
  # Extract the count value
  value <- df$count[1]
  
  # Append the value to the results
  corrected_cell_counts <- rbind(corrected_cell_counts, data.frame(file = file, count = value))
}

# Sort the results by the count in descending order
corrected_cell_counts <- corrected_cell_counts %>%
  arrange(desc(count))

# Find all instances where the count is equal to the maximum count
max_count <- max(corrected_cell_counts$count)
max_count_files <- corrected_cell_counts %>%
  filter(count == max_count)

# Output the final filtered dataframe with all max count cases
max_count_files




######################
##### neutral largest lost clone
######################



dir_path <- "/Volumes/feder-vol1/project/HomeostaticEpithelium/dat/2024-04-10-SingleInjection_Confluence_blockProb/results"
# list files
# List all files
file_list <- list.files(path = dir_path)

# Filter files to keep only those that start with 'VisFile_block_0_'
filtered_list <- grep("^VisFile_block_0_", file_list, value = TRUE)

lost_simulations <- filtered_list[!grepl("_(45|93|68)\\.[0-9]+\\.txt$", filtered_list)]

# Extract numeric values from the filenames for sorting
numbers <- as.numeric(gsub(".*\\.(\\d+)\\.txt$", "\\1", lost_simulations))

# Find the highest integer value and corresponding file
max_index <- which.max(numbers)
max_value <- numbers[max_index]
max_file <- lost_simulations[max_index]

# Output the file with the maximum number
print(max_file)

# Process all timepoints of the lost simulations
corrected_cell_counts <- data.frame(file = character(), count = integer(), stringsAsFactors = FALSE)

for (file in lost_simulations) {
  file_path <- file.path(dir_path, file)  
  column_names = c('x', 'z', 'y', 'h', 's', 'v', 'alpha', 'mut', 'injSite', 'tp53Clone')
  
  # Read the file
  data <- read_tsv(file_path, col_names = column_names)
  
  # Filter the data as needed and count the cells
  df <- data %>%
    filter(y == 0) %>%
    filter(injSite != -1) %>%
    summarise(count = n())
  
  # Extract the count value
  value <- df$count[1]
  
  # Append the count to the results
  corrected_cell_counts <- rbind(corrected_cell_counts, data.frame(file = file, count = value))
}

# Sort the results by the count in descending order
corrected_cell_counts <- corrected_cell_counts %>%
  arrange(desc(count))

# Output the final sorted dataframe
corrected_cell_counts