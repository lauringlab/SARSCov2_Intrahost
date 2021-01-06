

# Make this into an R script and include into Snakefile. Do for both filtered and unfiltered.

# ========================= Modules ============================

library(tidyverse)

# ========================= Combine variant files ============================

file_list_1 <- Sys.glob("data/ivar_output/variants/*.tsv")

index <- 1
for (file in file_list_1)
{
  print(index)
  # if the merged dataset doesn't exist, create it
  if(!exists("variants"))
  {
    variants <- read.table(file, header = TRUE, sep = "\t")
    filename_only <- as.character(str_split(file, "/")[[1]][3])
    id <- str_split(filename_only, "[.]")[[1]][1]
    variants <- mutate(variants, ID = id)
  }
  
  # if the merged dataset does exist, append to it
  if(exists("variants"))
  {
    temp_dataset <- read.table(file, header = TRUE, sep = "\t")
    filename_only <- as.character(str_split(file, "/")[[1]][3])
    id <- str_split(filename_only, "[.]")[[1]][1]
    temp_dataset <- mutate(temp_dataset, ID = id)
    
    variants <- rbind(variants, temp_dataset)
    rm(temp_dataset)
  }
  
  index = index + 1
}

write.csv(variants, "data/ivar_output/all.variants.csv")

# ========================= Combine variant files ============================

file_list_2 <- Sys.glob("data/ivar_output/variants_final/*.tsv")

index <- 1
for (file in file_list_2)
{
  print(index)
  # if the merged dataset doesn't exist, create it
  if(!exists("variants_filtered"))
  {
    variants_filtered <- read.table(file, header = TRUE, sep = "\t")
    filename_only <- as.character(str_split(file, "/")[[1]][3])
    id <- str_split(filename_only, "[.]")[[1]][1]
    variants_filtered <- mutate(variants_filtered, ID = id)
  }
  
  # if the merged dataset does exist, append to it
  if(exists("variants_filtered"))
  {
    temp_dataset <- read.table(file, header = TRUE, sep = "\t")
    filename_only <- as.character(str_split(file, "/")[[1]][3])
    id <- str_split(filename_only, "[.]")[[1]][1]
    temp_dataset <- mutate(temp_dataset, ID = id)
    
    variants_filtered <- rbind(variants_filtered, temp_dataset)
    rm(temp_dataset)
  }
  
  index = index + 1
}

write.csv(variants_filtered, "data/ivar_output/all.variants.filtered.csv")

