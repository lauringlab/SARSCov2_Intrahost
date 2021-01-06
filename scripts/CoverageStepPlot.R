

### Author: Andrew Valesano
### Purpose: Analysis of SARS-CoV-2 coverage data from Illumina sequencing with ARTIC V3 primers.

# ============================ Read in modules and coverage data =======================

library(tidyverse)

cov <- read.csv("data/raw/coverage.csv", stringsAsFactors = FALSE)

# ====================== Plot number of samples by fraction covered ====================

GetFraction <- function(cov_df, depth, step) 
{
  genome_len <- length(unique(cov_df$pos))
  
  cov_df %>%
    filter(cov >= depth) %>%
    group_by(ID) %>%
    summarize(fraction = n()/genome_len) %>%
    filter(fraction >= step) -> cov_df_threshold
  
  num_samples <- length(unique(cov_df_threshold$ID))
  
  return(num_samples)
}

GetStepwiseFraction <- function(cov_df, depth = 10, step_size = 0.01)
{
  steps <- seq(0, 1, by = step_size)
  
  fractions <- c()
  for(step in steps)
  {
    print(step)
    fraction <- GetFraction(cov_df, depth = depth, step)
    fractions <- c(fractions, fraction)
  }
  
  fraction_data <- data.frame(step = steps, num_samples = fractions)
  
  return(fraction_data)
}

genome_end <- max(cov$pos)
depth_thresholds <- c(10, 20, 200, 400, 500, 1000)

stepwise_data <- data.frame()
for(d in depth_thresholds)
{
  print(d)
  cov %>%
    filter(pos > 0 & pos < (genome_end - 0)) %>%
    do(GetStepwiseFraction(., depth = d, step_size = 0.01)) -> fractions
  
  fractions <- mutate(fractions, threshold = d)
  stepwise_data <- rbind(stepwise_data, fractions)
}

stepwise_data$threshold <- factor(stepwise_data$threshold, levels = c(10, 20, 200, 400, 500, 1000))

stepwise_data$num_samples[stepwise_data$step == 0] <- 325
stepwise_data_plot <- stepwise_data

# Color scheme
stepplot_colors <- c("#D86018", "#75988d", "#575294", "#FFCB05", "#00B2A9", "#702082", "#575294", "#655A52")

stepwise.plot <- ggplot(stepwise_data_plot, aes(x = step, y = num_samples, color = threshold)) +
  geom_step() +
  theme_bw() +
  xlab("Fraction of Genome Above Depth Cutoff") +
  ylab("Number of Samples") +
  ylim(c(0, length(unique(cov$ID)) + 1)) +
  scale_color_manual(name = "Cutoff", values = stepplot_colors) +
  theme(text = element_text(size = 10)) # PDF, 4 by 5

