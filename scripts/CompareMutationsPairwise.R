
### Author: Andrew Valesano
### Purpose: Compare pairwise shared iSNV.

# ========================= Load modules and data ============================

library(tidyverse)
library(ape)
library(reshape2)
library(patchwork)

metadata <- read.csv("data/metadata/metadata_called_variants.csv", stringsAsFactors = FALSE)
variants <- read.csv("data/processed/processed.variants.csv", stringsAsFactors = FALSE) 
dates <- read.csv("data/metadata/all_sample_dates.csv", stringsAsFactors = FALSE)
clusters <- read.csv("data/metadata/employee_clusters.csv", stringsAsFactors = FALSE) %>%
  mutate(cluster = as.character(cluster_ID),
         ID = as.character(ID))

# Read in distance matrix
fasta <- read.dna("data/tree/data/all.consensus.wref2.aln.mask.fasta", format = "fasta")
dist <-  dist.dna(fasta, model = "N", as.matrix = TRUE)
dist_df <- melt(as.matrix(dist), varnames = c("row", "col"))
dist_df <- mutate(dist_df, is_self = ifelse(row == col, TRUE, FALSE)) %>%
  mutate(ID1 = gsub("/", "_", row), ID2 = gsub("/", "_", col)) %>%
  select(-row, -col)

# =============================== Distance metrics from reference =========================

dist_df_ref <- filter(dist_df, ID1 == "Wuhan-Hu-1" & !ID2 %in% c("Wuhan-Hu-1", "Wuhan_WH01_2019"))

ggplot(dist_df_ref, aes(value)) +
  geom_histogram(binwidth = 1)

median(dist_df_ref$value)
min(dist_df_ref$value)
max(dist_df_ref$value)
IQR(dist_df_ref$value)

# ================================= How far apart are the employee clusters? ============================

dist_df_clusters <- filter(dist_df, ID1 %in% clusters$ID & ID2 %in% clusters$ID)

dist_df_clusters <- mutate(dist_df_clusters, cluster1 = clusters$cluster[match(ID1, clusters$ID)])
dist_df_clusters <- mutate(dist_df_clusters, cluster2 = clusters$cluster[match(ID2, clusters$ID)])
dist_df_clusters <- filter(dist_df_clusters, ID1 != ID2) %>% filter(cluster1 == cluster2)

dist_df_clusters_distinct <- dist_df_clusters[!duplicated(t(apply(dist_df_clusters,1,sort))),] # Some clusters are related, but not others.

# ================================== What about multiple iSNV? ===============================

IDs_A11782G <- filter(variants, mutation == "A11782G")$ID
IDs_G12331A <- filter(variants, mutation == "G12331A")$ID
IDs_T13914G <- filter(variants, mutation == "T13914G")$ID

dist_df_A11782G <- filter(dist_df, ID1 %in% IDs_A11782G & ID2 %in% IDs_A11782G)
dist_df_G12331A <- filter(dist_df, ID1 %in% IDs_G12331A & ID2 %in% IDs_G12331A)
dist_df_T13914G <- filter(dist_df, ID1 %in% IDs_T13914G & ID2 %in% IDs_T13914G)

dist_df_A11782G_distinct <- dist_df_A11782G[!duplicated(t(apply(dist_df_A11782G,1,sort))),]
dist_df_G12331A_distinct <- dist_df_G12331A[!duplicated(t(apply(dist_df_G12331A,1,sort))),]
dist_df_T13914G_distinct <- dist_df_T13914G[!duplicated(t(apply(dist_df_T13914G,1,sort))),]

# ================================== Distribution of pairwise shared variants =============================

variants_analyze_pt <- filter(variants, str_detect(ID, "S_2006"))
variants_analyze_hcw <- filter(variants, !str_detect(ID, "S_2006"))

IDs_all <- unique(metadata$ID)
IDs_pt <- unique(variants_analyze_pt$ID)
IDs_hcw <- unique(variants_analyze_hcw$ID)

IDs_use <- IDs_all

expand.grid(IDs_use, IDs_use) %>%
  mutate(specimen1 = as.character(Var1), specimen2 = as.character(Var2)) %>%
  select(-Var1, -Var2) %>%
  filter(specimen1 != specimen2) -> all_spec_comparisons

all_spec_comparisons_final <- all_spec_comparisons[!duplicated(t(apply(all_spec_comparisons, 1, sort))),]

all_spec_comparisons_final <- mutate(all_spec_comparisons_final, shared_var = NA)
for(r in 1:nrow(all_spec_comparisons_final))
{
  row <- all_spec_comparisons_final[r, ]
  sample1 <- unique(row$specimen1)
  sample2 <- unique(row$specimen2)
  #print(r)
  
  #variants_final_use <- filter(variants_final, !mutation %in% c("G11083T")) # We have already filtered out this position
  variants_final_use <- variants
  variants_samples <- filter(variants_final_use, ID %in% c(sample1, sample2))
  
  variants_samples %>%
    group_by(mutation) %>%
    summarize(count = n()) -> counts
  
  num_shared <- nrow(filter(counts, count > 1))
  all_spec_comparisons_final[r, ]$shared_var <- num_shared
}

all_spec_comparisons_final <- mutate(all_spec_comparisons_final, specimen1_group = ifelse(str_detect(specimen1, "S_2006"), "patient", "employee"))
all_spec_comparisons_final <- mutate(all_spec_comparisons_final, specimen2_group = ifelse(str_detect(specimen2, "S_2006"), "patient", "employee"))

# Groups
all_spec_comparisons_final <- mutate(all_spec_comparisons_final, group = NA)
all_spec_comparisons_final <- mutate(all_spec_comparisons_final, group = ifelse(specimen1_group == "patient" & specimen2_group == "patient", "patient - patient", group))
all_spec_comparisons_final <- mutate(all_spec_comparisons_final, group = ifelse(specimen1_group == "employee" & specimen2_group == "employee", "employee - employee", group))

all_spec_comparisons_final <- mutate(all_spec_comparisons_final, cluster1 = clusters$cluster[match(specimen1, clusters$ID)])
all_spec_comparisons_final <- mutate(all_spec_comparisons_final, cluster2 = clusters$cluster[match(specimen2, clusters$ID)])
all_spec_comparisons_final <- mutate(all_spec_comparisons_final, group = ifelse((!is.na(cluster1) & !is.na(cluster2)) & cluster1 == cluster2, "employee cluster", group))

all_spec_comparisons_final <- mutate(all_spec_comparisons_final, group = ifelse(is.na(group), "patient - employee", group))

# Combining consensus SNV data from dist_df
all_spec_comparisons_final <- mutate(all_spec_comparisons_final, pair = paste0(specimen1, "-", specimen2))
dist_df <- mutate(dist_df, pair = paste0(ID1, "-", ID2))
all_spec_comparisons_final <- mutate(all_spec_comparisons_final, fixed_diffs = dist_df$value[match(pair, dist_df$pair)])
all_spec_comparisons_final <- mutate(all_spec_comparisons_final, distance_group = ifelse(fixed_diffs <= 1, "Near-Identical", "Not Identical"))

# Add time difference
all_spec_comparisons_final <- mutate(all_spec_comparisons_final, date1 = dates$collection_date[match(specimen1, dates$ID)])
all_spec_comparisons_final <- mutate(all_spec_comparisons_final, date2 = dates$collection_date[match(specimen2, dates$ID)])
all_spec_comparisons_final <- mutate(all_spec_comparisons_final, time_diff = abs(as.Date(date2) - as.Date(date1)))

time.diff.hist <- ggplot(all_spec_comparisons_final, aes(time_diff)) +
  geom_histogram(binwidth = 1)

# Writing for future use
#write.csv(all_spec_comparisons_final, "data/processed/all_spec_comparisons_final.csv")
#all_spec_comparisons_final <- read.csv("data/processed/all_spec_comparisons_final.csv", stringsAsFactors = FALSE)

### Compare pairs within one week of each other ###

# Change employee clusters back to employee - employee. Since there are no shared iSNV.
all_spec_comparisons_final_plot <- mutate(all_spec_comparisons_final, group = ifelse(group == "employee cluster", "employee - employee", group))

all_spec_comparisons_final_withinweek <- filter(all_spec_comparisons_final_plot, time_diff <= 7)
all_spec_comparisons_final_outsideweek <- filter(all_spec_comparisons_final_plot, time_diff > 7)


### Now plot ###
shared.by.fixed.colors = c("patient - patient" = "#702082", 
                           "patient - employee" = "coral2", 
                           "employee - employee" = "#575294",
                           "employee cluster" = "palevioletred3")

shared.by.fixed <- ggplot(all_spec_comparisons_final_withinweek, aes(y = shared_var, x = fixed_diffs, color = group)) +
  geom_jitter(width = 0.2, alpha = 0.8, height = 0.1) +
  theme_classic() +
  xlab("Consensus Differences") +
  ylab("Shared iSNV") +
  facet_wrap(~group) +
  scale_color_manual(name = "Pair Type", values = shared.by.fixed.colors) +
  theme() +
  scale_y_continuous(breaks = seq(0, 6, 1), limits = c(-0.1, 3)) +
  theme(text = element_text(size = 10), 
        strip.background = element_rect(colour = "white", fill = "white"), 
        strip.text.x = element_text(face = "bold"),
        legend.position = "none") +
  ggtitle("Within week")

shared.by.fixed.outsideweek <- ggplot(all_spec_comparisons_final_outsideweek, aes(y = shared_var, x = fixed_diffs, color = group)) +
  geom_jitter(width = 0.2, alpha = 0.8, height = 0.1) +
  theme_classic() +
  ggtitle("Outside week") +
  xlab("Consensus Differences") +
  ylab("Shared iSNV") +
  facet_wrap(~group) +
  scale_color_manual(name = "Pair Type", values = shared.by.fixed.colors) +
  theme() +
  scale_y_continuous(breaks = seq(0, 6, 1), limits = c(-0.1, 3)) +
  theme(text = element_text(size = 10), 
        strip.background = element_rect(colour = "white", fill = "white"), 
        strip.text.x = element_text(face = "bold"),
        legend.position = "none")

shared.by.fixed / shared.by.fixed.outsideweek # PDF 6 by 10

time.diff.hist <- ggplot(all_spec_comparisons_final_plot, aes(time_diff, fill = group)) +
  geom_histogram(binwidth = 1) +
  scale_color_manual(name = "Pair Type", values = shared.by.fixed.colors) +
  theme_classic() +
  xlab("Days Between Samples in Pair") +
  ylab("Count") +
  facet_wrap(~group)
  

fixed.hist <- ggplot(all_spec_comparisons_final_withinweek, aes(fixed_diffs, fill = group)) +
  geom_histogram(binwidth = 1, color = "white") +
  theme_classic() +
  facet_wrap(~group) +
  xlab("Consensus Differences") +
  ylab("Number of Pairs") +
  theme(text = element_text(size = 10), 
        strip.background = element_rect(colour = "white", fill = "white"), 
        strip.text.x = element_text(face = "bold"),
        legend.position = "none") +
  scale_fill_manual(name = "Pair Type", values = shared.by.fixed.colors)

