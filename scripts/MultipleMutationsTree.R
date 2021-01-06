
# ============================= Modules and data ==========================

library(tidyverse)
library(ggtree)
library(ape)
library(treeio)
library(phytools)
library(lubridate)

tree <- read.nexus("data/tree/data/treetime_output/divergence_tree.nexus")
dates <- read.table("data/tree/data/all_sample_dates.tsv", sep = "\t", header = TRUE)

variants <- read.csv("data/processed/processed.variants.csv", stringsAsFactors = FALSE)
multiple_mutations <- read.csv("data/processed/multiple_mutations.csv", stringsAsFactors = FALSE)

pangolin <- read.csv("data/metadata/all.pangolin.results.csv", stringsAsFactors = FALSE) %>% 
  mutate(name = Sequence.name)
clades <- select(pangolin, name, Lineage) %>% 
  mutate(Group = ifelse(str_detect(name, "S_2006"), "hospitalized", "employee")) %>%
  mutate(Group = ifelse(name == "Wuhan-Hu-1" | name == "Wuhan/WH01/2019", "reference", Group))

clusters <- read.csv("data/metadata/employee_clusters.csv", stringsAsFactors = FALSE) %>%
  rename(name = ID) %>%
  mutate(cluster = as.character(cluster_ID),
         name = as.character(name))

clades <- left_join(clades, clusters, by = "name")
clades$cluster[is.na(clades$cluster)] <- "None"
clades$cluster[clades$cluster %in% c("8", "26", "18", "14", "12", "11")] <- "None" # If we only sequenced one in a cluster, then let's not highlight here.

clades <- mutate(clades, has_11782 = ifelse(name %in% filter(variants, mutation == "A11782G")$ID, "A11782G", "None"))
clades <- mutate(clades, has_12331 = ifelse(name %in% filter(variants, mutation == "G12331A")$ID, "G12331A", "None"))
clades <- mutate(clades, has_13914 = ifelse(name %in% filter(variants, mutation == "T13914G")$ID, "T13914G", "None"))

rownames(clades) <- clades$name

# ========================== Plot the tree =========================

tree.plot <- ggtree(tree, color = "black", as.Date = TRUE, mrsd = max(dates$date), ladderize = TRUE)

tree_point_colors <- c("reference" = "black", "hospitalized" = "#00B2A9", "employee" = "#575294")
mutation_colors <- c("A11782G" = "#702082", "G12331A" = "#FFCB05", "T13914G" = "#D86018", "None" = "#989C97")

tree.plot.meta <- tree.plot %<+% clades +
  geom_point(aes(color = Group), size = 1.5) +
  scale_color_manual(values = tree_point_colors) +
  theme_tree2() + 
  theme_bw() +
  theme(legend.position = "right",
        panel.grid.minor = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank()) +
  scale_x_date(date_breaks = "month", date_labels = "%y/%m/%d")
tree.plot.meta

tree.plot.meta.heat <- gheatmap(tree.plot.meta, clades[, c("has_11782", "has_12331", "has_13914"),  drop = FALSE], width = 0.05, color = "white", offset = 0) +
  scale_fill_manual(values = mutation_colors, name = "") +
  theme(legend.position = "right")
tree.plot.meta.heat # PDF 8 by 12

# For divergence tree
tree.plot <- ggtree(tree, color = "black")

tree.plot.div <- tree.plot %<+% clades +
  geom_point(aes(color = Group), size = 1.5) +
  scale_color_manual(values = tree_point_colors) +
  geom_treescale(x = 0, y = 30) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank())
tree.plot.div

tree.plot.div.heat <- gheatmap(tree.plot.div, clades[, c("has_11782", "has_12331", "has_13914"),  drop = FALSE], width = 0.05, color = "white", offset = 0) +
  scale_fill_manual(values = mutation_colors, name = "")
tree.plot.div.heat # PDF 8 by 12

