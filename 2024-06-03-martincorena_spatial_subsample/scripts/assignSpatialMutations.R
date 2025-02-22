## Hunter Colegrove
## FA dynamics
## 24 May 2024

## Spatially reconstruct a tissue section from martincorena 2018 esophagus mutation frequencies 

library(tidyverse)
library(glue)
library(ggplot2)
library(RColorBrewer)

project_path = "2024-06-03-martincorena_spatial_subsample"



#####################
###### Process martincorena file
#####################

filepath = "/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-02-24-p53_fitnessMagnitude/martincorena_table2_esophagus.txt"
#filepath = "/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-02-24-p53_fitnessMagnitude/martincorena_table2_esophagus.txt"
tp53 = read_delim(filepath, delim = "\t")
cell_density = 15000 ## average based on Colom, Jones, & Doupe data

## number of samples collected for each age range
nSamples <- c(84, 93, 96, 94, 95, 98, 95, 95, 94)

print(unique(tp53$impact))
## filter for SNV only
tp53 <- tp53 %>%
  filter(impact == "Missense" | impact == "Nonsense" | impact == "Synonymous")

## calculate clone size area of each SNV
tp53 <- tp53 %>% 
  mutate(area = vaf*2*2) %>% # 2 copies per cell*vaf = proportion of area * 2mm^2 per sample 
  dplyr::select(sampleID, vaf, area) %>%
  group_by(sampleID) %>%
  mutate(mutation = row_number()) %>%
  arrange(sampleID, desc(vaf)) %>%
  mutate(size = area*cell_density) 

## add age ranges to data
tp53 <- tp53 %>%
  mutate(age = if_else(grepl("PD36806", sampleID), "20-23",
               if_else(grepl("PD36712", sampleID), "24-27",
               if_else(grepl("PD30272", sampleID), "36-39",
               if_else(grepl("PD30986", sampleID), "44-47",
               if_else(grepl("PD30987", sampleID), "48-51",
               if_else(grepl("PD30274", sampleID), "52-55",
               if_else(grepl("PD30988", sampleID), "56-59",
               if_else(grepl("PD30273", sampleID), "68-71",
               if_else(grepl("PD31182", sampleID), "72-75",
               NA_character_
               ))))))))))


## filter muts > 55 years old
tp53 <- tp53 %>% filter(!(age == "56-59" | age == "68-71" | age == "72-75"))

a <- print(unique(tp53$sampleID))
length(a)
#write.csv(sampled_mutations, glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/IM_tp53_area.csv"), row.names = FALSE)


#####################
###### Using martincorena mutation file, place mutations in spatial locations on a 2mm^2 grid
#####################

########
### Distance Function
########

precompute_neighbors <- function(size) {
  neighbors_list <- vector("list", size * size)
  for (x in 1:size) {
    for (y in 1:size) {
      idx <- (x - 1) * size + y
      directions <- tibble(
        x = c(-1, 1, 0, 0),
        y = c(0, 0, -1, 1)
      )
      neighbors <- tibble(
        x = x + directions$x,
        y = y + directions$y
      )
      neighbors <- neighbors %>% filter(x >= 1 & y >= 1 & x <= size & y <= size)
      neighbors_list[[idx]] <- neighbors
    }
  }
  return(neighbors_list)
}

get_precomputed_neighbors <- function(x, y, size, neighbors_list) {
  idx <- (x - 1) * size + y
  return(neighbors_list[[idx]])
}

# function to find moore neighbors 
# filter by positions within the grid

## von-nuemann
# get_neighbors <- function(x, y, size) {
#   directions <- tibble(
#     x = c(-1, 1, 0, 0),
#     y = c(0, 0, -1, 1)
#   )
#   neighbors <- tibble(
#     x = x + directions$x,
#     y = y + directions$y
#   )
#   neighbors <- neighbors %>% filter(x >= 1 & y >= 1 & x <= size & y <= size)
#   return(neighbors)
# }

########
### Function to place mutations on grid
########

place_mutations <- function(grid, mutation_df, size, seed) {
  #set.seed(seed)
  print(mutation_df)
  ## iterate through each mutation
  for(i in 1:nrow(mutation_df)){
    print(i)
    current_mutation <- mutation_df[i,]
    ## find empty cells on the grid
    empty_cells <- grid %>% filter(grid$Mut1 == 0)
    empty_cells_count <- nrow(empty_cells)
  
    ## if there aren't enough spaces, stack mutations 
    if(empty_cells_count < current_mutation$size){
      existing_muts <- grid %>% filter(grid$Mut1 != 0 & grid$Mut2 == 0)
      random_start <- existing_muts[sample(nrow(existing_muts), 1), ]
      placed_cells <- tibble(x = random_start$x, y = random_start$y)
      
      ## iterate through the rest of cells within this mutation and place it within a nearby cell containing the same Mut1 value
      for(j in 2:current_mutation$size){
        # find neighboring cells of all placed cells
        
        neighbor_cells <- do.call(rbind, lapply(1:nrow(placed_cells), function(k) {
          get_precomputed_neighbors(placed_cells$x[k], placed_cells$y[k], size, neighbors_list)
        }))
        
        neighbor_cells <- neighbor_cells %>%
          distinct() %>%
          left_join(grid, by = c("x", "y")) %>%
          filter(Mut1 == random_start$Mut1)
        
        # neighbor_cells <- placed_cells %>%
        #   rowwise() %>%
        #   do(get_neighbors(.$x, .$y, size)) %>%
        #   ungroup() %>%
        #   distinct()
        
        # filter for neighbors containing the same Mut1 value
        # neighbor_cells <- neighbor_cells %>%
        #   left_join(grid, by = c("x", "y")) %>%
        #   filter(Mut1 == random_start$Mut1)
  
        # If no Mut1 neighbors, throw error
        if(nrow(neighbor_cells) == 0) {
          cat("Error: no available neighboring cells containing the same Mut1 value\n")
          break
        }
        
        # select a random neighbor to place the next mutation
        selected_neighbor <- neighbor_cells[sample(nrow(neighbor_cells), 1), ]
        grid <- grid %>%
          mutate(Mut2 = case_when(x == selected_neighbor$x & y == selected_neighbor$y ~ current_mutation$mutation, TRUE ~ Mut2))
        
        # add to list of placed cells
        placed_cells <- placed_cells %>%
          add_row(x = selected_neighbor$x, y = selected_neighbor$y)
      }
    }
    
    
    ## if there is space, add mutation to empty locations
    else{
      random_start <- empty_cells[sample(nrow(empty_cells), 1), ]
      grid <- grid %>%
        mutate(Mut1 = case_when(x == random_start$x & y == random_start$y ~ current_mutation$mutation, TRUE ~ Mut1))
      placed_cells <- tibble(x = random_start$x, y = random_start$y)
      
      ## iterate through the rest of cells within this mutation and place it within a nearby cell of the current mutation
      for(j in 2:current_mutation$size){
        print(j)
        
        neighbor_cells <- do.call(rbind, lapply(1:nrow(placed_cells), function(k) {
          get_precomputed_neighbors(placed_cells$x[k], placed_cells$y[k], size, neighbors_list)
        }))
        
        # filter for empty neighbors
        neighbor_cells <- neighbor_cells %>%
          distinct() %>%
          left_join(grid, by = c("x", "y")) %>%
          filter(Mut1 == 0)
        
        # find neighboring cells of all placed cells
        # neighbor_cells <- placed_cells %>%
        #   rowwise() %>%
        #   do(get_neighbors(.$x, .$y, size)) %>%
        #   ungroup() %>%
        #   distinct()
        # 
        # # filter for empty neighbors
        # neighbor_cells <- neighbor_cells %>%
        #   left_join(grid, by = c("x", "y")) %>%
        #   filter(Mut1 == 0)
        
        # if no empty neighbors throw error
        if(nrow(neighbor_cells) == 0) {
          cat("Error: no available empty neighboring cells to place mutation\n")
          break
        }
        
        # select a random neighbor to place the next mutation
        selected_neighbor <- neighbor_cells[sample(nrow(neighbor_cells), 1), ]
        grid <- grid %>%
          mutate(Mut1 = case_when(x == selected_neighbor$x & y == selected_neighbor$y ~ current_mutation$mutation, TRUE ~ Mut1))
        
        # add to list of placed cells
        placed_cells <- placed_cells %>%
          add_row(x = selected_neighbor$x, y = selected_neighbor$y)
      }
    }
  }
  return(grid)

}


########
### Main
########

## Define grid size
grid_size <- 173
#grid_size <- 20

## set seed
seed_val <- 13
set.seed(seed_val)

## Generate empty grid
grid <- tibble(
  x = rep(1:grid_size, each = grid_size),
  y = rep(1:grid_size, times = grid_size),
  Mut1 = 0,
  Mut2 = 0
)

## precompute neighbors
neighbors_list <- precompute_neighbors(grid_size)

##### Mutations dataframe

## Martincorena
for(i in unique(tp53$sampleID)){
  mutation_df <- tp53 %>% filter(sampleID == i)
  id_age <- mutation_df$age
  grid_with_mutations <- place_mutations(grid, mutation_df, grid_size, seed = seed_val)
  print(grid_with_mutations)
  #write.csv(grid_with_mutations, glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/IM_tp53_area_{i}_{id_age}.csv"), row.names = FALSE)
  write.csv(grid_with_mutations, glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/reconstructed_tissue/IM_tp53_area_{i}_{id_age}.csv"), row.names = FALSE)

  
  #####################
  ###### Visualize
  #####################
  
  ## factor dataframe
  grid_with_mutations <- grid_with_mutations %>%
    mutate(Mut1 = as.factor(grid_with_mutations$Mut1)) %>%
    mutate(Mut2 = as.factor(grid_with_mutations$Mut2))
  
  # Mut1 plot
  mut1 <- ggplot(grid_with_mutations, aes(x = x, y = y, fill = Mut1)) +
    geom_tile() +
    scale_fill_brewer(palette = "Set3") +
    labs(title = "Mutation Distribution on 173x173 Grid",
         x = "Column",
         y = "Row",
         fill = "Mutation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    coord_fixed()
  
  # Mut2 plot
  mut2 <- ggplot(grid_with_mutations, aes(x = x, y = y, fill = Mut2)) +
    geom_tile() +
    scale_fill_brewer(palette = "Set3") +
    labs(title = "Mutation Distribution on 173x173 Grid",
         x = "Column",
         y = "Row",
         fill = "Mutation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    coord_fixed() 
  
  # Plot Mut1 + Mut2
  df0 <- grid_with_mutations
  df0_long <- pivot_longer(df0, cols = -c(x, y), names_to = "mut")
  df <- df0_long
  df1    <- df[!duplicated(interaction(df$x, df$y)),]
  df2    <- df[duplicated(interaction(df$x, df$y)),]
  df2    <- df[rep(seq(nrow(df)), each = 3),]
  df2$x1 <- as.numeric(as.factor(df2$x))
  df2$y1 <- as.numeric(as.factor(df2$y))
  df2$x1 <- df2$x1 + c(-0.5, 0.5, 0.5)
  df2$y1 <- df2$y1 + c(-0.5, -0.5, 0.5)
  df2$z  <- rep(seq(nrow(df2)/3), each = 3)
  
  mut12 <- ggplot(df1, aes(x = x, y = y, fill = value)) + 
    geom_tile() +
    scale_fill_brewer(palette = "Set3") +
    geom_polygon(data = df2, aes(x = x1, y = y1, group = z))
  
  
  ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/reconstructed_tissue/{i}_{id_age}_mut1_plot.png"), mut1, width = 10, height = 5, units = "in", dpi = 300)
  ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/reconstructed_tissue/{i}_{id_age}_mut2_plot.png"), mut2, width = 10, height = 5, units = "in", dpi = 300)
  ggsave(glue("/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/reconstructed_tissue/{i}_{id_age}_mut1_2_plot.png"), mut12, width = 10, height = 5, units = "in", dpi = 300)
  
  }


## Generate example dataframe of mutations and sizes
# mutation_df <- data.frame(
#   mutation = c(1, 2),
#   size = sample(1:300, 2, replace = TRUE)
# )

## Place mutations on the grid
#grid_with_mutations <- place_mutations(grid, mutation_df, grid_size, seed = seed_val)

