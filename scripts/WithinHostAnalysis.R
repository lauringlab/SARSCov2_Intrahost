
### Author: Andrew Valesano
### Purpose: Analysis of within-host variation in clinical SARS-CoV-2 samples.

# ========================= Load modules and data ============================

library(tidyverse)
library(patchwork)

variants <- read.csv("data/raw/all.variants.csv", stringsAsFactors = FALSE) 
coverage <- read.csv("data/raw/coverage.csv", stringsAsFactors = FALSE) %>% mutate(ID = gsub("/", "_", ID))
metadata <- read.csv("data/metadata/sample_qpcr_dpso.csv", stringsAsFactors = FALSE)

# ========================== Filter by coverage depth and evenness =================================

consensus_depth_threshold = 10
consensus_length_threshold = 29000
variant_threshold = 200 
variant_completeness_threshold = 0.8

coverage %>%
  group_by(ID) %>%
  filter(cov >= consensus_depth_threshold) %>%
  filter(n() >= consensus_length_threshold) -> coverage_consensus

coverage %>%
  group_by(ID) %>%
  filter(cov >= variant_threshold) %>%
  filter(n() >= max(coverage$pos)*variant_completeness_threshold) -> coverage_variant

variants_analyze <- filter(variants, ID %in% coverage_consensus$ID & ID %in% coverage_variant$ID)

# ========================== Filter by copy number, based on mixing experiment =====================

metadata_quality_copynumber <- filter(metadata, N1_copynum_final > 1e3 & !is.na(DPSO) & DPSO > -10)

variants_analyze <- filter(variants_analyze, ID %in% metadata_quality_copynumber$ID) 

length(unique(variants_analyze$ID)) # Down to 178 samples

# ================================= Plot coverage ===============================

coverage_varquality <- filter(coverage, ID %in% variants_analyze$ID)

coverage_varquality %>%
  group_by(pos) %>%
  summarize(percentile_5th = quantile(cov, probs = c(0.05)),
            median = median(cov),
            percentile_95th = quantile(cov, probs = c(0.95))) -> coverage_summary

cov.plot <- ggplot() +
  geom_ribbon(data = coverage_summary, aes(x = pos, ymin = percentile_5th, ymax = percentile_95th), fill = "#989C97", alpha = 0.5) +
  geom_line(data = coverage_summary, aes(x = pos, y = median), color = "#575294", size = 0.3) +
  scale_y_log10() +
  theme_classic() +
  xlab("Genome Position") +
  ylab("Read Depth") +
  theme(text = element_text(size = 10)) # PDF, 4 by 6

coverage_varquality_masked <- filter(coverage_varquality, pos > 55 & pos < 29753)
coverage_varquality_masked %>% 
  group_by(ID) %>%
  summarize(mean = mean(cov)) -> coverage_varquality_masked_means

ggplot(coverage_varquality_masked_means, aes(mean)) +
  geom_histogram()


# =================================== Process and filter variants ======================================

variants_analyze %>% 
  filter(ALT_FREQ < 0.5) %>%
  select(-X) %>%
  mutate(REF_FWD = REF_DP - REF_RV, ALT_FWD = ALT_DP - ALT_RV) %>%
  mutate(mutation = paste0(REF, POS, ALT)) %>%
  filter(!str_detect(string = ALT, pattern = "\\+") & !str_detect(string = ALT, pattern = "\\-")) %>%
  filter(TOTAL_DP > 100 & PVAL < 0.0001 & ALT_QUAL > 35 & ALT_DP >= 10) -> variants_filtered

# strand bias filter
variants_filtered$strand.p.val <- apply(variants_filtered, 1, function(x) 
{
  data <- x[c("REF_FWD", "ALT_FWD", "REF_RV", "ALT_RV")]
  tbl <- matrix(as.numeric(data), ncol = 2, byrow = FALSE)
  fisher.test(tbl, alternative = "two.sided")$p.value
})
variants_filtered <- mutate(variants_filtered, strand.p.val.adj = p.adjust(strand.p.val, method = "bonferroni", n = nrow(variants_filtered)))
variants_filtered_strandbias <- filter(variants_filtered, strand.p.val.adj > 0.05)
variants_final <- filter(variants_filtered_strandbias, ALT_FREQ > 0.02)

# Get coding information.
variants_final <- mutate(variants_final, type = NA)
variants_final <- mutate(variants_final, type = ifelse(is.na(GFF_FEATURE), "Noncoding", type))
variants_final <- mutate(variants_final, type = ifelse(!is.na(GFF_FEATURE) & REF_AA == ALT_AA, "Synonymous", type))
variants_final <- mutate(variants_final, type = ifelse(!is.na(GFF_FEATURE) & REF_AA != ALT_AA, "Non-synonymous", type))
variants_final <- mutate(variants_final, type = ifelse(!is.na(GFF_FEATURE) & ALT_AA == "*", "Stop", type))

variants_final_processed <- filter(variants_final, type != "Stop" & POS != 11083) # Remove position 11083, error-prone

#write.csv(variants_final_processed, "data/processed/processed.variants.csv")
#variants_final_processed <- read.csv("data/processed/processed.variants.csv", stringsAsFactors = FALSE)

# ================================= Plot variant coding changes ================================

freq.pos.palette <- c("Synonymous" = "#575294", "Non-synonymous" = "coral2", "Noncoding" = "#989C97")
freq.pos <- ggplot(filter(variants_final_processed, type != "Stop"), aes(x = POS, y = ALT_FREQ, fill = type)) +
  geom_point(shape = 21, size = 2) +
  theme_bw() +
  xlab("Genome Position") +
  ylab("Frequency") +
  scale_fill_manual(name = "Type", values = freq.pos.palette)

variants_quality_compare_NS_S <- filter(variants_final_processed, type %in% c("Non-synonymous", "Synonymous"))
freq.histogram <- ggplot(variants_quality_compare_NS_S, aes(x = ALT_FREQ, fill = type)) + 
  geom_histogram(color = "white", binwidth = 0.05, position = "dodge", boundary = 0.02) +
  xlab("Frequency") + 
  ylab("Number of Minor iSNV") +
  scale_fill_manual(name = "" , values = freq.pos.palette) +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", text = element_text(size = 10))

variants_quality_compare_NS_S$GFF_FEATURE <- factor(variants_quality_compare_NS_S$GFF_FEATURE, levels = c("ORF1a", "ORF1b", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF8", "N"))

mutation.type.by.gene <- ggplot(variants_quality_compare_NS_S, aes(x = GFF_FEATURE, fill = type)) + 
  geom_bar(position = "dodge") +
  xlab("Gene in coding region") + 
  ylab("Number of Minor iSNV") +
  scale_fill_manual(name = "" , values = freq.pos.palette) +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 10)) 

coding.changes.plot <- freq.histogram | mutation.type.by.gene # 5 by 12

# =================================== Relationship of iSNV richness to viral load and DPSO ======================================

metadata_variants <- filter(metadata, ID %in% variants_analyze$ID) 
nrow(metadata_variants) # 178
#write.csv(metadata_variants, "data/metadata/metadata_called_variants.csv")

variants_final_processed %>%
  group_by(ID) %>%
  summarize(iSNV = n()) -> minor_counts

left_join(select(metadata_variants, ID, DPSO, N1_copynum_final, group), minor_counts, by = "ID") -> minor_counts_meta
minor_counts_meta$iSNV[is.na(minor_counts_meta$iSNV)] <- 0 # If none were found, make it zero instead of NA

richness.hist <- ggplot(minor_counts_meta, aes(iSNV)) + 
  geom_histogram(binwidth = 1, fill = "#575294", color = "white") +
  theme_classic() +
  xlab("Number of Minor iSNV") +
  ylab("Number of Specimens") +
  scale_x_continuous(breaks = c(seq(0, 50, 5))) +
  theme(text = element_text(size = 10)) # PDF 4 by 4

median(minor_counts_meta$iSNV)
quantile(minor_counts_meta$iSNV)

# iSNV by copy number

snv_plot_colors <- c("hospitalized" = "#00B2A9", "employee" = "#575294")
snv.by.copies <- ggplot(minor_counts_meta, aes(x = N1_copynum_final, y = iSNV, color = group)) + 
  geom_point(size = 1.5) + 
  ylab("Minor iSNV Per Sample") + 
  theme_bw() +
  xlab("Genome Copies") +
  scale_x_log10() +
  theme(text = element_text(size = 10), legend.position = "none") +
  scale_color_manual(values = snv_plot_colors) # PDF 4 by 5

# Look at high-richness samples

samples_high_richness <- filter(minor_counts_meta, iSNV > 15)

freq.pos.high.richness <- ggplot(filter(variants_final, ID %in% samples_high_richness$ID & type != "Stop"), aes(x = POS, y = ALT_FREQ, fill = type)) +
  geom_point(shape = 21, size = 2) +
  theme_bw() +
  xlab("Genome Position") +
  ylab("Frequency") +
  scale_fill_manual(name = "Type", values = freq.pos.palette) +
  facet_wrap(~ID) +
  ylim(c(0, 0.5))

variants_final_lowrichness <- filter(variants_final, !ID %in% samples_high_richness$ID)
write.csv(variants_final_lowrichness, "data/processed/processed.variants.lowrichness.csv")

variants_highrichness <- filter(variants_analyze, ID %in% samples_high_richness$ID) %>%
  filter(ALT_FREQ > 0.02) %>%
  filter(!str_detect(string = ALT, pattern = "\\+") & !str_detect(string = ALT, pattern = "\\-")) %>%
  filter(TOTAL_DP > 100 & PVAL < 0.0001 & ALT_QUAL > 35 & ALT_DP >= 10)

freq.pos.high.richness <- ggplot(variants_highrichness, aes(x = POS, y = ALT_FREQ)) +
  geom_point(shape = 21, size = 2) +
  theme_bw() +
  xlab("Genome Position") +
  ylab("Frequency") +
  scale_fill_manual(name = "Type", values = freq.pos.palette) +
  facet_wrap(~ID)

# iSNV by DPSO

max_isnv <- 15
snv.by.day <- ggplot() + 
  geom_point(data = filter(minor_counts_meta, iSNV <= max_isnv), aes(x = DPSO, y = iSNV, fill = group), size = 2.5, alpha = 0.9, shape = 21, position = position_jitter(width = 0.2, height = 0.1)) +
  ylab("Minor iSNV Per Sample") + 
  theme_classic() +
  xlab("Day Post Symptom Onset") +
  scale_fill_manual(name = "", values = c("#575294", "#00B2A9")) +
  theme(text = element_text(size = 15)) +
  scale_x_continuous(breaks = c(-5, 0, 5, 10, 15, 20, 25)) +
  scale_y_continuous(breaks = seq(0, 50, 10)) +
  theme(text = element_text(size = 10), legend.position = "none") # PDF 4 by 8 for full plot, 4 by 6 for max_isnv of 15

snv.model <- lm(iSNV ~ log(N1_copynum_final, 10) + DPSO, data = minor_counts_meta)
summary(snv.model)

# No difference in iSNV richness between the two groups. p = 0.29.
wilcox.test(iSNV ~ group, data = minor_counts_meta)

compare.groups.plot <- ggplot() +
  geom_jitter(data = minor_counts_meta, aes(x = as.factor(group), y = iSNV, color = group), alpha = 0.7, width = 0.25) +
  geom_boxplot(data = minor_counts_meta, aes(x = as.factor(group), y = iSNV), alpha = 0) +
  scale_color_manual(name = "", values = c("#575294", "#00B2A9")) +
  geom_boxplot(alpha = 0, outlier.shape = NULL) +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("") +
  ylab("Minor iSNV Per Sample") +
  theme(text = element_text(size = 10)) # PDF 4 by 4
