

### Author: Andrew Valesano
### Purpose: Analyze sensitivity and specificity of iSNV identification for SARS-CoV-2 (ARTIC on Illumina).

# =============================== Modules and data ========================

library(tidyverse)
library(cowplot)

variants <- read.csv("data/mixing/all.variants.final.csv", stringsAsFactors = FALSE) 
coverage <- read.csv("data/mixing/coverage.csv")
expectedMutations <- read.csv("data/mixing/ExpectedMutations.csv", stringsAsFactors = FALSE)
groups <- read.csv("data/mixing/MixGroups.csv", stringsAsFactors = FALSE) %>% rename(ID = samplename)
left_join(variants, groups, by = "ID") -> variants_meta

# ======================= Examine coverage by input copy number =====================

coverage %>%
  mutate(ID = gsub(pattern = "/", "_", ID)) -> coverage
left_join(coverage, groups, by = "ID") -> coverage_meta

coverage_meta %>%
  group_by(copies, pos) %>%
  summarize(mean_cov = mean(cov)) %>%
  mutate(copies = as.character(copies)) -> mean_cov_by_copies

copy_colors <- c("100" = "#D86018", 
                    "100000" = "#75988d", 
                    "1000" = "#00274C", 
                    "10000" = "#00B2A9")

cov.plot <- ggplot(mean_cov_by_copies, aes(x = pos, y = mean_cov, color = as.factor(copies))) +
  geom_line() +
  xlab("Genome Position") +
  ylab("Mean Coverage") +
  theme_bw() +
  scale_y_log10() +
  scale_color_manual(values = copy_colors, name = "Copies")

# ======================== Variant processing =============================

variants_meta %>% 
  select(-X) %>%
  mutate(REF_FWD = REF_DP - REF_RV, ALT_FWD = ALT_DP - ALT_RV) %>%
  mutate(mutation = paste0(REF, POS, ALT)) %>%
  filter(!str_detect(string = ALT, pattern = "\\+") & !str_detect(string = ALT, pattern = "\\-")) %>%
  mutate(expected = ifelse(mutation %in% expectedMutations$mutation, TRUE, FALSE)) %>%
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
variants_final <- filter(variants_filtered_strandbias, ALT_FREQ > 0.02) # 5% cutoff gets rid of most FPs in high-titer samples. Is there a smarter way?

# =============================== Variant plots ==========================

variants_meta %>% 
  select(-X) %>%
  mutate(REF_FWD = REF_DP - REF_RV, ALT_FWD = ALT_DP - ALT_RV) %>%
  mutate(mutation = paste0(REF, POS, ALT)) %>%
  filter(!str_detect(string = ALT, pattern = "\\+") & !str_detect(string = ALT, pattern = "\\-") & ALT_FREQ > 0.01) %>%
  mutate(expected = ifelse(mutation %in% expectedMutations$mutation, TRUE, FALSE)) -> variants_unfiltered

set2 <- brewer.pal(n = 8, name = 'Set2')
variant.plot <- ggplot(variants_final, aes(x = POS, y = ALT_FREQ, fill = expected)) +
  geom_point(alpha = 1, shape = 21, size = 2) +
  xlab("Genome Position") +
  ylab("Frequency") +
  theme_bw() +
  ylim(c(0, 0.4)) +
  facet_wrap(~copies) +
  scale_fill_manual(name = "", values = c(set2[2], set2[3]))

variant.plot.log <- ggplot(variants_final, aes(x = POS, y = ALT_FREQ, color = expected)) +
  geom_point() +
  xlab("Genome Position") +
  ylab("Frequency") +
  theme_bw() +
  facet_wrap(~copies) +
  scale_y_log10() +
  scale_color_manual(name = "", values = c(wes_D2[2], wes_D1[1]))

qual.by.pval <- ggplot(filter(variants_unfiltered, copies > 10), aes(x = as.factor(ALT_QUAL), y = PVAL, color = expected)) +
  geom_jitter(width = 0.2) +
  theme_bw() +
  scale_color_manual(name = "", values = c(wes_D2[2], wes_D1[1]))

# ================================= Plot for figure =============================

tf.plot.colors <- c("TRUE" = "#575294", "FALSE" =  "coral2")

variant.plot.e2 <- ggplot(filter(variants_final, copies == 100), aes(x = POS, y = ALT_FREQ, fill = expected)) +
  geom_point(alpha = 1, shape = 21, size = 2) +
  xlab("Genome Position") +
  ylab("") +
  theme_bw() +
  ylim(c(0, 0.4)) +
  scale_fill_manual(name = "", values = tf.plot.colors) +
  ggtitle("100 copies") +
  theme(legend.position = "none")

variant.plot.e3 <- ggplot(filter(variants_final, copies == 1000), aes(x = POS, y = ALT_FREQ, fill = expected)) +
  geom_point(alpha = 1, shape = 21, size = 2) +
  xlab("Genome Position") +
  ylab("Frequency") +
  theme_bw() +
  ylim(c(0, 0.4)) +
  scale_fill_manual(name = "", values = tf.plot.colors) +
  ggtitle("1000 copies") +
  theme(legend.position = "none")

variant.plot.e4 <- ggplot(filter(variants_final, copies == 10000), aes(x = POS, y = ALT_FREQ, fill = expected)) +
  geom_point(alpha = 1, shape = 21, size = 2) +
  xlab("") +
  ylab("") +
  theme_bw() +
  ylim(c(0, 0.4)) +
  scale_fill_manual(name = "", values = tf.plot.colors) +
  ggtitle("10000 copies") +
  theme(legend.position = "none")

variant.plot.e5 <- ggplot(filter(variants_final, copies == 100000), aes(x = POS, y = ALT_FREQ, fill = expected)) +
  geom_point(alpha = 1, shape = 21, size = 2) +
  xlab("") +
  ylab("Frequency") +
  theme_bw() +
  ylim(c(0, 0.4)) +
  scale_fill_manual(name = "", values = tf.plot.colors) +
  ggtitle("100000 copies") +
  theme(legend.position = "none")

specificity.plot <- plot_grid(variant.plot.e5, variant.plot.e4, variant.plot.e3, variant.plot.e2) # PDF 7 by 10

# =========================== Sensitivity and specificity quick look =======================

variants_filtered_strandbias %>%
  group_by(ID) %>%
  filter(expected == TRUE) %>%
  summarize(true_positives = n()) -> sensitivity

#write.csv(sensitivity, "data/processed/mixing_sensitivity.csv")

left_join(sensitivity, groups, by = "ID") %>%
  filter(copies > 100) -> sensitivity_groups

#arrange(sensitivity_groups, -copies, -frequency) %>% select(copies, frequency, true_positives) %>% View()

variants_final %>%
  filter(expected == FALSE) %>%
  group_by(mutation) %>%
  summarize(count = n()) -> FP_counts # why are some positions piling up?

variants_final %>%
  group_by(ID) %>%
  filter(expected == FALSE) %>%
  summarize(FP = n()) -> false_positives

left_join(false_positives, groups, by = "ID") -> false_positives_meta

fp.by.copies <- ggplot(false_positives_meta, aes(x = as.factor(copies), y = FP)) +
  geom_jitter(width = 0.1, color = "coral2") +
  geom_boxplot(alpha = 0) +
  theme_bw() +
  scale_y_log10() +
  xlab("Copies") +
  ylab("False Positives Per Sample") # PDF 4 by 4

# ============================ Observed by expected for true positives =========================

variants_final_TP <- filter(variants_filtered_strandbias, expected == TRUE) # no frequency filter applied here

observed.by.expected <- ggplot(variants_final_TP, aes(y = ALT_FREQ, x = frequency)) +
  geom_point() +
  theme_bw() +
  xlab("Expected Frequency") +
  ylab("Observed Frequency") +
  facet_wrap(~copies) +
  geom_abline(slope = 1, linetype = "dotted", color = "blue")

observed.by.expected.low.colors <- c("#D86018", "#702082", "#FFCB05", "#00B2A9", "#575294")
observed.by.expected.low <- ggplot(filter(variants_final_TP, frequency < 1), aes(y = ALT_FREQ, x = frequency, fill = as.factor(frequency))) +
  geom_point(size = 3, shape = 21) +
  theme_bw() +
  xlab("Expected Frequency") +
  ylab("Observed Frequency") +
  facet_wrap(~copies) +
  geom_abline(slope = 1, linetype = "dotted", color = "black") +
  ylim(c(0, 0.25)) +
  xlim(c(0, 0.25)) +
  scale_fill_manual(name = "Frequency", values = observed.by.expected.low.colors) +
  theme(strip.background = element_rect(colour = "white", fill = "white"), strip.text.x = element_text(face = "bold")) # PDF 6 by 8

# For figure:

observed.by.expected.low.colors <- c("0.005" = "#D86018", "0.01" =  "#702082", "0.02" = "#FFCB05", "0.05" = "#00B2A9", "0.1" = "#575294")
observed.by.expected.low.e5 <- ggplot(filter(variants_final_TP, frequency < 1 & copies == 100000), aes(y = ALT_FREQ, x = frequency, fill = as.factor(frequency))) +
  geom_point(size = 3, shape = 21) +
  theme_bw() +
  xlab("Expected Frequency") +
  ylab("Observed Frequency") +
  facet_wrap(~copies) +
  geom_abline(slope = 1, linetype = "dotted", color = "black") +
  ylim(c(0, 0.25)) +
  xlim(c(0, 0.25)) +
  scale_fill_manual(name = "Frequency", values = observed.by.expected.low.colors) +
  theme(strip.background = element_rect(colour = "white", fill = "white"), strip.text.x = element_text(face = "bold"), legend.position = "none") # PDF 6 by 8

observed.by.expected.low.e4 <- ggplot(filter(variants_final_TP, frequency < 1 & copies == 10000), aes(y = ALT_FREQ, x = frequency, fill = as.factor(frequency))) +
  geom_point(size = 3, shape = 21) +
  theme_bw() +
  xlab("Expected Frequency") +
  ylab("Observed Frequency") +
  facet_wrap(~copies) +
  geom_abline(slope = 1, linetype = "dotted", color = "black") +
  ylim(c(0, 0.25)) +
  xlim(c(0, 0.25)) +
  scale_fill_manual(name = "Frequency", values = observed.by.expected.low.colors) +
  theme(strip.background = element_rect(colour = "white", fill = "white"), strip.text.x = element_text(face = "bold"), legend.position = "none") # PDF 6 by 8

observed.by.expected.low.e3 <- ggplot(filter(variants_final_TP, frequency < 1 & copies == 1000), aes(y = ALT_FREQ, x = frequency, fill = as.factor(frequency))) +
  geom_point(size = 3, shape = 21) +
  theme_bw() +
  xlab("Expected Frequency") +
  ylab("Observed Frequency") +
  facet_wrap(~copies) +
  geom_abline(slope = 1, linetype = "dotted", color = "black") +
  ylim(c(0, 0.25)) +
  xlim(c(0, 0.25)) +
  scale_fill_manual(name = "Frequency", values = observed.by.expected.low.colors) +
  theme(strip.background = element_rect(colour = "white", fill = "white"), strip.text.x = element_text(face = "bold"), legend.position = "none") # PDF 6 by 8

observed.by.expected.low.e2 <- ggplot(filter(variants_final_TP, frequency < 1 & copies == 100), aes(y = ALT_FREQ, x = frequency, fill = as.factor(frequency))) +
  geom_point(size = 3, shape = 21) +
  theme_bw() +
  xlab("Expected Frequency") +
  ylab("Observed Frequency") +
  facet_wrap(~copies) +
  geom_abline(slope = 1, linetype = "dotted", color = "black") +
  ylim(c(0, 0.25)) +
  xlim(c(0, 0.25)) +
  scale_fill_manual(name = "Frequency", values = observed.by.expected.low.colors) +
  theme(strip.background = element_rect(colour = "white", fill = "white"), strip.text.x = element_text(face = "bold"), legend.position = "none") # PDF 6 by 8

sensitivity.plot <- plot_grid(observed.by.expected.low.e5, 
                              observed.by.expected.low.e4, 
                              observed.by.expected.low.e3, 
                              observed.by.expected.low.e2) # PDF 6 by 8



