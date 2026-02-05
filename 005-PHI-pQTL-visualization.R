#!/usr/bin/env Rscript
# 005-PHI-pQTL-visualization.R
# Visualize Log_BETA distribution grouped by PHI
# Categorize PHI into 4 groups: 0-0.25, 0.25-0.5, 0.5-0.75, 0.75-1
# Create scientific-style violin plots and perform significance tests

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggsignif)

# ============================================================================
# 1. Configuration parameters
# ============================================================================

# Input file path (merged PHI and pQTL data)
input_file <- "pQTL-PHI/PHI_pQTL_merged_12_22.csv"

# Output directory
output_dir <- "pQTL-PHI"
figs_dir <- file.path(output_dir, "Figs")

# ============================================================================
# 2. Data loading and preprocessing
# ============================================================================

cat("=== PHI-pQTL Visualization Analysis ===\n\n")
cat("Loading data files...\n")

# Read merged data
if (!file.exists(input_file)) {
  stop(sprintf("Error: Input file not found: %s\nPlease run merge_phi_pqtl.py script first", input_file))
}

df <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)

cat(sprintf("Data loaded: %d rows, %d columns\n", nrow(df), ncol(df)))
cat("Column names:\n")
print(colnames(df))

# Check required columns
required_cols <- c("PHI", "Log_BETA")
missing_cols <- required_cols[!required_cols %in% colnames(df)]
if (length(missing_cols) > 0) {
  stop(sprintf("Error: Required columns not found: %s\n", paste(missing_cols, collapse = ", ")))
}

# Create output directory
if (!dir.exists(figs_dir)) {
  dir.create(figs_dir, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", figs_dir))
}

# ============================================================================
# 3. Create PHI groups
# ============================================================================

cat("\nCreating PHI groups...\n")

# Remove missing values
df_clean <- df %>%
  filter(!is.na(PHI) & !is.na(Log_BETA)) %>%
  mutate(
    PHI_bin = case_when(
      PHI >= 0 & PHI < 0.25 ~ "0-0.25",
      PHI >= 0.25 & PHI < 0.5 ~ "0.25-0.5",
      PHI >= 0.5 & PHI < 0.75 ~ "0.5-0.75",
      PHI >= 0.75 & PHI <= 1 ~ "0.75-1",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(PHI_bin))

# Set group order
df_clean$PHI_bin <- factor(df_clean$PHI_bin, 
                           levels = c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1"))

cat(sprintf("Cleaned data: %d rows\n", nrow(df_clean)))
cat("Sample counts per group:\n")
print(table(df_clean$PHI_bin))

# ============================================================================
# 4. Perform inter-group significance tests
# ============================================================================

cat("\nPerforming inter-group significance tests...\n")

# Extract Log_BETA values for each group
group_0_25 <- df_clean$Log_BETA[df_clean$PHI_bin == "0-0.25"]
group_25_50 <- df_clean$Log_BETA[df_clean$PHI_bin == "0.25-0.5"]
group_50_75 <- df_clean$Log_BETA[df_clean$PHI_bin == "0.5-0.75"]
group_75_100 <- df_clean$Log_BETA[df_clean$PHI_bin == "0.75-1"]

# Perform pairwise comparisons (using t-test)
comparisons <- list(
  c("0-0.25", "0.25-0.5"),
  c("0.25-0.5", "0.5-0.75"),
  c("0.5-0.75", "0.75-1"),
  c("0-0.25", "0.5-0.75"),
  c("0.25-0.5", "0.75-1"),
  c("0-0.25", "0.75-1")
)

# Store test results
test_results <- data.frame(
  Group1 = character(),
  Group2 = character(),
  P_value = numeric(),
  Statistic = numeric(),
  stringsAsFactors = FALSE
)

for (comp in comparisons) {
  group1_data <- df_clean$Log_BETA[df_clean$PHI_bin == comp[1]]
  group2_data <- df_clean$Log_BETA[df_clean$PHI_bin == comp[2]]
  
  if (length(group1_data) >= 2 && length(group2_data) >= 2) {
    # Perform t-test
    test_result <- t.test(group1_data, group2_data)
    p_val <- test_result$p.value
    stat_val <- test_result$statistic
    
    test_results <- rbind(test_results, data.frame(
      Group1 = comp[1],
      Group2 = comp[2],
      P_value = p_val,
      Statistic = as.numeric(stat_val),
      stringsAsFactors = FALSE
    ))
    
    # Format p-value display
    if (p_val < 0.001) {
      p_display <- "P < 0.001"
    } else if (p_val < 0.01) {
      p_display <- "P < 0.01"
    } else if (p_val < 0.05) {
      p_display <- "P < 0.05"
    } else {
      p_display <- sprintf("P = %.3f", p_val)
    }
    
    cat(sprintf("%s vs %s: %s (t=%.3f, n1=%d, n2=%d)\n", 
                comp[1], comp[2], p_display, stat_val, 
                length(group1_data), length(group2_data)))
  }
}

# Save test results
test_results_file <- file.path(output_dir, "PHI_bin_statistical_tests.csv")
write.csv(test_results, test_results_file, row.names = FALSE)
cat(sprintf("\nTest results saved to: %s\n", test_results_file))

# ============================================================================
# 5. Create violin plots
# ============================================================================

cat("\nGenerating violin plots...\n")

# Calculate median and mean for each group for annotation
summary_stats <- df_clean %>%
  group_by(PHI_bin) %>%
  summarise(
    n = n(),
    median = median(Log_BETA, na.rm = TRUE),
    mean = mean(Log_BETA, na.rm = TRUE),
    .groups = 'drop'
  )

cat("Summary statistics for each group:\n")
print(summary_stats)

# Define colors (using scientific-style color palette)
color_palette <- c(
  "0-0.25" = "#E8E8E8",      # Light gray
  "0.25-0.5" = "#B0D4E8",    # Light blue
  "0.5-0.75" = "#5BA3D0",    # Medium blue
  "0.75-1" = "#1E5F8F"        # Dark blue
)

# Create violin plot
p_violin <- ggplot(df_clean, aes(x = PHI_bin, y = Log_BETA, fill = PHI_bin)) +
  geom_violin(trim = FALSE, alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.8, 
               outlier.alpha = 0.5, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = color_palette) +
  labs(
    x = "PHI bin",
    y = "pQTL-log10(abs(BETA))"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    plot.margin = margin(15, 15, 15, 15)
  )

# Add significance annotations (comparisons between adjacent groups)
# Calculate y-axis maximum for annotation positioning
y_max <- max(df_clean$Log_BETA, na.rm = TRUE)
y_range <- max(df_clean$Log_BETA, na.rm = TRUE) - min(df_clean$Log_BETA, na.rm = TRUE)

p_violin <- p_violin +
  geom_signif(
    comparisons = list(
      c("0-0.25", "0.25-0.5"),
      c("0.25-0.5", "0.5-0.75"),
      c("0.5-0.75", "0.75-1")
    ),
    map_signif_level = function(p) {
      if (p < 0.001) return("***")
      if (p < 0.01) return("**")
      if (p < 0.05) return("*")
      return("ns")
    },
    textsize = 4,
    y_position = y_max + y_range * c(0.15, 0.25, 0.35),
    step_increase = 0
  )

# Save plot (reduced width)
output_file <- file.path(figs_dir, "violin_PHI_bins_Log_BETA.png")
ggsave(output_file, plot = p_violin, width = 6, height = 4, dpi = 300)
cat(sprintf("Saved: %s\n", output_file))

# ============================================================================
# 6. Generate version with more significance annotations (all inter-group comparisons)
# ============================================================================

cat("\nGenerating full significance annotation version...\n")

# Create version with all comparisons
p_violin_full <- ggplot(df_clean, aes(x = PHI_bin, y = Log_BETA, fill = PHI_bin)) +
  geom_violin(trim = FALSE, alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.8, 
               outlier.alpha = 0.5, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = color_palette) +
  labs(
    x = "PHI bin",
    y = "pQTL-log10(abs(BETA))"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    plot.margin = margin(15, 15, 15, 15)
  )

# Add all significance annotations
y_max_full <- max(df_clean$Log_BETA, na.rm = TRUE)
y_range_full <- max(df_clean$Log_BETA, na.rm = TRUE) - min(df_clean$Log_BETA, na.rm = TRUE)

p_violin_full <- p_violin_full +
  geom_signif(
    comparisons = comparisons,
    map_signif_level = function(p) {
      if (p < 0.001) return("***")
      if (p < 0.01) return("**")
      if (p < 0.05) return("*")
      return("ns")
    },
    textsize = 3.5,
    y_position = y_max_full + y_range_full * seq(0.15, 0.5, length.out = length(comparisons)),
    step_increase = 0,
    tip_length = 0.01
  )

# Save full version (reduced width)
output_file_full <- file.path(figs_dir, "violin_PHI_bins_Log_BETA_full_comparisons.png")
ggsave(output_file_full, plot = p_violin_full, width = 6, height = 5, dpi = 300)
cat(sprintf("Saved: %s\n", output_file_full))

# ============================================================================
# 7. Generate summary statistics table
# ============================================================================

cat("\nGenerating summary statistics...\n")

# Create detailed summary statistics
detailed_stats <- df_clean %>%
  group_by(PHI_bin) %>%
  summarise(
    n = n(),
    mean = mean(Log_BETA, na.rm = TRUE),
    median = median(Log_BETA, na.rm = TRUE),
    sd = sd(Log_BETA, na.rm = TRUE),
    se = sd / sqrt(n),
    q25 = quantile(Log_BETA, 0.25, na.rm = TRUE),
    q75 = quantile(Log_BETA, 0.75, na.rm = TRUE),
    min = min(Log_BETA, na.rm = TRUE),
    max = max(Log_BETA, na.rm = TRUE),
    .groups = 'drop'
  )

# Save summary statistics
stats_file <- file.path(output_dir, "PHI_bin_summary_statistics.csv")
write.csv(detailed_stats, stats_file, row.names = FALSE)
cat(sprintf("Summary statistics saved to: %s\n", stats_file))

cat("\nCompleted!\n")
