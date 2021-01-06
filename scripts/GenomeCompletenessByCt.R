

### Author: Andrew Valesano
### Purpose: Compare genome completeness to input copy number.


# ================================ Modules and data ================================

library(tidyverse)

coverage <- read.csv("data/raw/coverage.csv", stringsAsFactors = FALSE)
metadata <- read.csv("data/metadata/sample_qpcr_dpso.csv", stringsAsFactors = FALSE) %>%
  filter(!is.na(DPSO) & DPSO > -10)

# ================================ Genome copies by DPSO =======================

ct_plot_colors <- c("hospitalized" = "#00B2A9", "employee" = "#575294") # Seahawks green: #69be28

copies.by.day.jitter <- ggplot(filter(metadata, N1_copynum_final > 0), aes(x = as.factor(DPSO), y = N1_copynum_final, color = group)) +
  geom_jitter(width = 0.15, size = 2, height = 0.05, alpha = 1) +
  theme_classic() +
  xlab("Day Post Symptom Onset") +
  ylab("Genome Copies (N1 qPCR)") +
  scale_color_manual(name = "", values = ct_plot_colors) +
  scale_y_continuous(breaks = c(0, 1e2, 1e4, 1e6, 1e8, 1e10), trans = "log10") +
  theme(text = element_text(size = 10)) # PDF, 4 by 6

copies.by.day.point <- ggplot() +
  geom_point(data = filter(metadata, N1_copynum_final > 0), aes(x = DPSO, y = N1_copynum_final, color = group)) +
  theme_classic() +
  xlab("Day Post Symptom Onset") +
  ylab("Genome Copies (N1 qPCR)") +
  #geom_smooth(data = filter(metadata, N1_copynum_final > 0), aes(x = DPSO, y = N1_copynum_final), method = "lm", color = "black") +
  scale_color_manual(name = "", values = ct_plot_colors) +
  scale_y_continuous(breaks = c(0, 1e2, 1e4, 1e6, 1e8, 1e10), trans = "log10") +
  theme(text = element_text(size = 10)) # PDF, 4 by 6

# Negatively correlated with DPSO
shedding.model <- lm(log(N1_copynum_final, 10) ~ DPSO, data = filter(metadata, N1_copynum_final > 0))
summary(shedding.model)

# =================================== Relationship of coverage to Ct, all samples ======================================

GetFraction <- function(cov_df, depth) 
{
  genome_len <- length(unique(cov_df$pos))
  
  cov_df %>%
    filter(cov >= depth) %>%
    group_by(ID) %>%
    summarize(fraction = n()/genome_len) -> cov_df_threshold
  
  return(cov_df_threshold)
}

# Expand to more depth cutoffs: 20x, 200x, 500x, 1000x
genome_end <- max(coverage$pos)
depth_thresholds <- c(10, 50, 100, 400, 500, 1000)

stepwise_data <- data.frame()
for(d in depth_thresholds)
{
  print(d)
  coverage %>%
    filter(pos > 100 & pos < (genome_end - 50)) %>%
    do(GetFraction(., depth = d)) -> fractions
  
  fractions <- mutate(fractions, threshold = d)
  stepwise_data <- rbind(stepwise_data, fractions)
}

stepwise_data$threshold <- factor(stepwise_data$threshold, levels = c(10, 50, 100, 400, 500, 1000))

stepwise_data_10x <- filter(stepwise_data, threshold == 10) %>% 
  mutate(ID = gsub(pattern = "/", "_", ID)) %>% 
  rename(fraction_10x = fraction) %>% 
  select(-threshold)

metadata_cov_data <- left_join(select(metadata, ID, N1_copynum_final, group), stepwise_data_10x, by = "ID")
metadata_cov_data <- filter(metadata_cov_data, !is.na(N1_copynum_final))

point_colors <- c("hospitalized" = "#00B2A9", "employee" = "#575294")
completeness.by.copies <- ggplot(metadata_cov_data, aes(x = N1_copynum_final, y = fraction_10x, color = group)) +
  geom_point() +
  xlab("Genome Copies") +
  ylab("Genome Completeness") +
  theme_classic() +
  scale_x_continuous(breaks = c(0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9), trans = "log10") +
  scale_color_manual(name = "", values = point_colors) +
  theme(legend.position = "right") +
  theme(text = element_text(size = 10)) # PDF 4 by 6


