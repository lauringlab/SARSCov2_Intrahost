

### Author: Andrew Valesano
### Purpose: Find iSNV in multiple samples and compare to global frequency.

# ========================= Load modules and data ============================

library(tidyverse)
library(lubridate)

metadata <- read.csv("data/metadata/sample_qpcr_dpso.csv", stringsAsFactors = FALSE)
variants <- read.csv("data/processed/processed.variants.csv", stringsAsFactors = FALSE)
dates <- read.csv("data/metadata/all_sample_dates.csv", stringsAsFactors = FALSE)

variants_dates <- left_join(variants, dates, by = "ID")

# =================================== Get iSNV found in multiple samples ======================================

variants %>%
  group_by(mutation) %>%
  summarize(count_all = n()) %>%
  mutate(position = as.numeric(substr(mutation, 2, nchar(mutation)-1)),
         ref = substr(mutation, 1, 1),
         alt = substr(mutation, nchar(mutation), nchar(mutation))) -> pos_counts

left_join(pos_counts, select(variants, mutation, type), by = "mutation") %>%
  unique() -> pos_counts_type

pos_counts_type %>%
  group_by(count_all) %>%
  summarize(count = n()) -> pos_counts_type_summary

# Location of mutations and their counts.
multiple.mutation.colors <- c("Synonymous" = "#575294", "Non-synonymous" = "coral2", "Noncoding" = "#989C97")
multiple.mutation.plot <- ggplot() +
  geom_bar(data = filter(pos_counts_type, count_all >= 2), aes(x = position, y = count_all, fill = type), stat = "identity", width = 80) +
  geom_point(data = filter(pos_counts_type, count_all >= 2), aes(x = position, y = count_all), size = 2) +
  theme_classic() +
  xlab("Genome Position") +
  ylab("Number of Samples") +
  scale_fill_manual(name = "", values = multiple.mutation.colors) +
  scale_y_continuous(breaks = seq(1, 15)) # PDF 4 by 8

mutations <- filter(pos_counts_type, count_all >= 3)

write.csv(mutations, "data/processed/multiple_mutations.csv", quote = FALSE, row.names = FALSE)

# ========================= Analyze counts of parallel iSNV in GISAID =======================

mutations <- read.csv("data/processed/multiple_mutations.csv", stringsAsFactors = FALSE)
global_ids <- read.csv("data/processed/multiple_mutations_gisaid/global_aln_ids.txt", stringsAsFactors = FALSE)
meta <- read.table("../SARSCoV2_Sequencing/data/gisaid/downloaded_data/20201111/metadata_2020-11-11_07-27.tsv", header = TRUE, fill = TRUE, quote = "", sep = "\t", stringsAsFactors = FALSE) %>% 
  filter(region %in% c("Africa", "North America", "South America", "Asia", "Europe", "Oceania")) %>%
  filter(strain %in% global_ids$ID) %>%
  filter(date > "2019-12-27") %>%
  mutate(days_relative = as.Date(date) - as.Date("2019-12-26"))


# Process the data

GetData <- function(mutation, meta)
{
  filename = paste0("data/processed/multiple_mutations_gisaid/", mutation, "_genomes.csv")
  genomes <- read.csv(filename, stringsAsFactors = FALSE, header = FALSE) %>%
    select(-V1) %>%
    filter(V2 != 0) %>%
    rename(ID = V2)
  
  meta <- mutate(meta, has_mutation = ifelse(strain %in% genomes$ID, TRUE, FALSE))
  meta <- mutate(meta, date = as.Date("2019-12-26") + days_relative)
  meta <- mutate(meta, week = week(date)) %>% filter(week < 50)
  meta <- mutate(meta, floor_week = floor_date(date, unit = "week"))
  
  
  filter(meta, has_mutation == TRUE) %>% group_by(floor_week) %>% summarise(num_true = n()) -> data_true
  filter(meta, has_mutation == FALSE) %>% group_by(floor_week) %>% summarise(num_false = n()) -> data_false
  left_join(data_false, data_true, by = "floor_week") -> data_full
  data_full <- mutate(data_full, num_true = ifelse(is.na(num_true), 0, num_true))
  data_full <- mutate(data_full, frequency = num_true / num_false)
  data_full <- mutate(data_full, total = num_true + num_false)
  
  data_full <- mutate(data_full, mutation = mutation)
  
  return(data_full)
}

data_complete <- data.frame()
for(r in 1:nrow(mutations))
{
  row <- mutations[r,]
  mut <- unique(row$mutation)
  mut_lookup <- substr(mut, 2, nchar(mut))
  data <- GetData(mut_lookup, meta)
  data_complete <- rbind(data_complete, data)
}


plot_colors <- c("11782G" = "#702082", "12331A" = "#FFCB05", "13914G" = "#D86018")
local.global.plot <- ggplot(filter(data_complete, total > 100), aes(x = floor_week, y = frequency, color = mutation)) +
  geom_line(size = 0.8) +
  theme_classic() +
  xlab("") +
  ylab("Frequency") +
  ylim(c(0, 0.1)) +
  theme(legend.position = "right") +
  scale_color_manual(name = "", values = plot_colors) +
  scale_x_date(date_breaks = "month", date_labels = "%m/%y") +
  geom_vline(xintercept = as.Date("2020-03-20"), linetype = "dotted", color = "#702082", size = 0.8) +
  geom_vline(xintercept = as.Date("2020-03-30"), linetype = "dotted", color = "#FFCB05", size = 0.8) +
  geom_vline(xintercept = as.Date("2020-03-25"), linetype = "dotted", color = "#D86018", size = 0.8) # PDF 3 by 5.
  
# Results: None rise above 1% frequency on a weekly basis from start of pandemic to 11/11.
# 11782 first found on 2020-03-20.
# 12331 first found on 2020-03-30.
# 13914 first found on 2020-03-25.

# None are found above 50% frequency.
