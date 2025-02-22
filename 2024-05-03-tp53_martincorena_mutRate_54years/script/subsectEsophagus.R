## Hunter Colegrove
## FA dynamics
## 21 Mar 2024

## 2mm^2 tissue sections from martincorena data randomly sampled into 9 subsections.


library(tidyverse)
library(glue)

project_path = "2024-05-03-tp53_martincorena_mutRate_54years"

random_sample_sections <- function(data, nsubsections) {
  ## initialize output lists
  sampleID = list()
  subsampleID = list()
  clone = list()
  size = list()
  vaf = list()
  ## Iterate through each sample in dataset
  unique_sampleID = unique(data$sampleID)
  #print(length(unique_sampleID))
  
  #count = 0
  for(sample_id in unique_sampleID){
    #print(glue("sample_id: {sample_id}"))
    #count = count + 1
    
    ## initalize key-value pair for subsection and its occupied area to be updated
    ## throughout each mutation
    key_value <- tibble(
      subsection = seq(1:nsubsections),
      occupied_area = rep(0,nsubsections)
    )
    #print(glue("sample_id: {sample_id}"))
    
    ## grab the ith sample
    sample_data = data[data$sampleID == sample_id, ]
    print(sample_data)
    ## Iterate through each mutation in sample 
    for (i in seq_len(nrow(sample_data))) {
      #print(glue("i: {i}"))
      mutation <- sample_data[i, ]
      mutation_area <- mutation$area ## of each mutation
      #print(glue("mutation_area: {mutation_area}"))
      
      ## sample a random subsection and place the mutation there
      while(mutation_area > 0){
        
        ## grab a random subsection
        random_subsection = sample(1:nsubsections,1)
        
        ## find area occupied by this subsection
        occupied_area = key_value$occupied_area[key_value$subsection == random_subsection]
        
        ## find how much area this subsection has available
        remaining_area = (2/nsubsections) - occupied_area #### 2mm^2/nsubsections
        #print(glue("remaining area: {remaining_area}"))
        
        if(remaining_area > 0){
          ## If space, allocate the entire mutation, else allocate what space is available
          allocate_area = min(mutation_area, remaining_area)
          #print(glue("allocate_area: {allocate_area}"))
          ## update the occupied area of this subsection
          key_value$occupied_area[random_subsection] <- occupied_area + allocate_area
          
          sampleID <- c(sampleID, sample_id)
          subsampleID <- c(subsampleID, random_subsection)
          clone <- c(clone, i)
          size <- c(size, allocate_area)
          vaf <- c(vaf, allocate_area/(4/nsubsections))

          ## Update remaining mutation area
          mutation_area = mutation_area - allocate_area
          #print(glue("available mutation_area: {mutation_area}"))
        }
        if(remaining_area == 0){
          if(sum(sum(key_value$occupied_area) > 1.99)){
            break
          }else{
            next
          }
          
        }
      }
    }
  }
  #print(glue("sampleID length: {length(sampleID)}"))
  #print(glue("muation length: {length(clone)}"))
  #print(glue("size length: {length(size)}"))
  result <- tibble(
    sampleID = sampleID,
    subsampleID = subsampleID,
    clone = clone,
    size = size, 
    vaf = vaf
  )
  return(result)

}

############
#### MAIN
############


# Function to randomly sample mutations into sections
filepath = "/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/2024-02-24-p53_fitnessMagnitude/martincorena_table2_esophagus.txt"
tp53 = read_delim(filepath, delim = "\t")
#print(tp53)

# number of samples collected for each age range
nSamples <- c(84, 93, 96, 94, 95, 98, 95, 95, 94)

#print(tp53)
## filter for SNV only
tp53 <- tp53 %>%
  filter(impact == "Missense" | impact == "Nonsense" | impact == "Synonymous")

## calculate clone size area of each SNV
tp53 <- tp53 %>% 
  mutate(area = vaf*2*2) %>% # 2 copies per cell*vaf = proportion of area * 2mm^2 per sample 
  select(sampleID, vaf, area) %>%
  group_by(sampleID)
print(tp53, n=Inf)

# Randomly sample mutations into sections
n_subsections = 9
sampled_mutations <- random_sample_sections(tp53, n_subsections)

# Output sampled mutations
#print(sampled_mutations)
sampled_mutations <- unnest(sampled_mutations, cols = c(sampleID, subsampleID, clone, size, vaf))

## add age ranges based on sample ID
sampled_mutations <- sampled_mutations %>%
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
print(tp53)
tp53 %>% filter(age == "44-47")

# Print the resulting dataframe
print(sampled_mutations)
sampled_mutations %>% filter(age == "44-47")
write.csv(sampled_mutations, glue("/Volumes/gs-vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/dat/subsect_clones_tp53.csv"), row.names = FALSE)


