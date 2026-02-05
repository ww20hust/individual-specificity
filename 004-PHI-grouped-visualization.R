#!/usr/bin/env Rscript
# 004-PHI-grouped-visualization.R
# Visualize grouped correlation results from script 003
# Includes frequency distribution plots, violin plots, paired tests, and difference bar charts

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(cowplot)

# ============================================================================
# 1. Configuration parameters (modifiable)
# ============================================================================

# Set visit pair to analyze (must match VISIT1 and VISIT2 in script 003)
VISIT1 <- "02"  # First visit
VISIT2 <- "22"  # Second visit

# Data paths
OUTPUT_DIR <- "./results"
INPUT_SUBDIR <- "age-sex-seperate-result"
OUTPUT_SUBDIR <- "age-sex-seperate-result/Figs"

# ============================================================================
# 2. Data loading and preprocessing
# ============================================================================

cat("=== 004-PHI Grouped Correlation Visualization ===\n\n")
cat(sprintf("Analyzing visit pair: %s vs %s\n", VISIT1, VISIT2))
cat("Loading data files...\n")

# Build input file path
pair_name <- paste(VISIT1, VISIT2, sep = "_")
input_file <- file.path(OUTPUT_DIR, INPUT_SUBDIR, 
                        sprintf("003-PHI-grouped-correlation-%s-results.csv", pair_name))

# Read result file
if (!file.exists(input_file)) {
  stop(sprintf("Error: Input file not found: %s\nPlease run 003-PHI-grouped-correlation.R script first", input_file))
}

df <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)

cat(sprintf("Data loaded: %d rows, %d columns\n", nrow(df), ncol(df)))
cat("Column names:\n")
print(colnames(df))

# Function: Find column name (handle R's automatic X prefix addition)
find_column <- function(df, pattern) {
  # Try exact match first
  if (pattern %in% colnames(df)) {
    return(pattern)
  }
  # Try adding X prefix
  x_pattern <- paste0("X", pattern)
  if (x_pattern %in% colnames(df)) {
    return(x_pattern)
  }
  # Try case-insensitive match
  matches <- grep(paste0("^", pattern, "$"), colnames(df), ignore.case = TRUE, value = TRUE)
  if (length(matches) > 0) {
    return(matches[1])
  }
  return(NULL)
}

# Identify column names
sex_male_corr <- find_column(df, "Sex_Male_corr")
sex_female_corr <- find_column(df, "Sex_Female_corr")
age_40_65_corr <- find_column(df, "Age_40-65_corr")
age_65_90_corr <- find_column(df, "Age_65-90_corr")
entrez_col <- find_column(df, "EntrezGeneSymbol")

# Verify required columns exist
required_cols <- list(
  "Sex_Male_corr" = sex_male_corr,
  "Sex_Female_corr" = sex_female_corr,
  "Age_40-65_corr" = age_40_65_corr,
  "Age_65-90_corr" = age_65_90_corr,
  "EntrezGeneSymbol" = entrez_col
)

missing_cols <- names(required_cols)[sapply(required_cols, is.null)]
if (length(missing_cols) > 0) {
  stop(sprintf("Error: Required columns not found: %s\n", paste(missing_cols, collapse = ", ")))
}

cat("All required columns found\n")

# Create output directory
figs_dir <- file.path(OUTPUT_DIR, OUTPUT_SUBDIR)
if (!dir.exists(figs_dir)) {
  dir.create(figs_dir, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", figs_dir))
}

# ============================================================================
# 3. Frequency distribution plots
# ============================================================================

cat("\nGenerating frequency distribution plots...\n")

# Function to create distribution plot
create_distribution_plot <- function(values, title, xlabel, output_file) {
  # Remove missing values
  values_clean <- values[!is.na(values)]
  
  if (length(values_clean) == 0) {
    cat(sprintf("Warning: %s has no valid data, skipping\n", title))
    return(NULL)
  }
  
  plot_data <- data.frame(Value = values_clean)
  
  p <- ggplot(plot_data, aes(x = Value)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
    labs(x = xlabel, y = "Count", title = title) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 15),
      plot.title = element_text(size = 16, hjust = 0.5),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  ggsave(output_file, plot = p, width = 6, height = 4, dpi = 300)
  cat(sprintf("Saved: %s\n", output_file))
  return(p)
}

# Create distribution plots for each group
create_distribution_plot(
  df[[sex_male_corr]], 
  "Male Group PHI Distribution",
  "PHI",
  file.path(figs_dir, "distribution_sex_male.png")
)

create_distribution_plot(
  df[[sex_female_corr]], 
  "Female Group PHI Distribution",
  "PHI",
  file.path(figs_dir, "distribution_sex_female.png")
)

create_distribution_plot(
  df[[age_40_65_corr]], 
  "Age 40-65 Group PHI Distribution",
  "PHI",
  file.path(figs_dir, "distribution_age_40-65.png")
)

create_distribution_plot(
  df[[age_65_90_corr]], 
  "Age 65-90 Group PHI Distribution",
  "PHI",
  file.path(figs_dir, "distribution_age_65-90.png")
)

# ============================================================================
# 4. Sex group violin plots and paired tests
# ============================================================================

cat("\nGenerating sex group violin plots...\n")

# Prepare long format data
sex_data <- data.frame(
  Protein = 1:nrow(df),
  Male = df[[sex_male_corr]],
  Female = df[[sex_female_corr]],
  stringsAsFactors = FALSE
)

# Convert to long format
sex_data_long <- sex_data %>%
  tidyr::pivot_longer(cols = c(Male, Female), 
                      names_to = "Group", 
                      values_to = "Correlation") %>%
  filter(!is.na(Correlation))

# Perform Wilcoxon paired test
# Only use proteins with both Male and Female data
valid_pairs <- !is.na(sex_data$Male) & !is.na(sex_data$Female)
if (sum(valid_pairs) >= 2) {
  wilcox_result <- wilcox.test(sex_data$Male[valid_pairs], 
                                sex_data$Female[valid_pairs], 
                                paired = TRUE)
  p_value <- wilcox_result$p.value
  statistic_value <- wilcox_result$statistic
  
  # Format p-value display
  if (p_value < 0.001) {
    p_display <- "P < 0.001"
  } else if (p_value < 0.01) {
    p_display <- "P < 0.01"
  } else if (p_value < 0.05) {
    p_display <- "P < 0.05"
  } else {
    p_display <- sprintf("P = %.3f", p_value)
  }
  
  # Save test results
  sex_test_result <- data.frame(
    Comparison = "Sex (Male vs Female)",
    Test = "Wilcoxon signed-rank test",
    Statistic = as.numeric(statistic_value),
    P_value = p_value,
    N_pairs = sum(valid_pairs),
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("Wilcoxon paired test result: %s (n=%d pairs)\n", p_display, sum(valid_pairs)))
} else {
  p_display <- "Insufficient data"
  sex_test_result <- data.frame(
    Comparison = "Sex (Male vs Female)",
    Test = "Wilcoxon signed-rank test",
    Statistic = NA,
    P_value = NA,
    N_pairs = sum(valid_pairs),
    stringsAsFactors = FALSE
  )
  cat("Warning: Insufficient paired data for testing\n")
}

# Create violin plot
p_violin_sex <- ggplot(sex_data_long, aes(x = Group, y = Correlation, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
  scale_fill_manual(values = c("Male" = "#4A90E2", "Female" = "#E24A90")) +
  labs(x = "Group", y = "PHI", 
       title = "Sex Group Comparison") +
  annotate("text", 
           x = 1.5, y = Inf,
           label = p_display,
           hjust = 0.5, vjust = 1.2,
           size = 5) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(file.path(figs_dir, "violin_sex_comparison.png"), 
       plot = p_violin_sex, width = 5, height = 6, dpi = 300)
cat(sprintf("Saved: %s\n", file.path(figs_dir, "violin_sex_comparison.png")))

# ============================================================================
# 5. Age group violin plots and paired tests
# ============================================================================

cat("\nGenerating age group violin plots...\n")

# Prepare long format data
age_data <- data.frame(
  Protein = 1:nrow(df),
  Age_40_65 = df[[age_40_65_corr]],
  Age_65_90 = df[[age_65_90_corr]],
  stringsAsFactors = FALSE
)

# Convert to long format
age_data_long <- age_data %>%
  tidyr::pivot_longer(cols = c(Age_40_65, Age_65_90), 
                     names_to = "Group", 
                     values_to = "Correlation") %>%
  mutate(Group = ifelse(Group == "Age_40_65", "40-65", "65-90")) %>%
  filter(!is.na(Correlation))

# Perform Wilcoxon paired test
valid_pairs_age <- !is.na(age_data$Age_40_65) & !is.na(age_data$Age_65_90)
if (sum(valid_pairs_age) >= 2) {
  wilcox_result_age <- wilcox.test(age_data$Age_40_65[valid_pairs_age], 
                                    age_data$Age_65_90[valid_pairs_age], 
                                    paired = TRUE)
  p_value_age <- wilcox_result_age$p.value
  statistic_value_age <- wilcox_result_age$statistic
  
  # Format p-value display
  if (p_value_age < 0.001) {
    p_display_age <- "P < 0.001"
  } else if (p_value_age < 0.01) {
    p_display_age <- "P < 0.01"
  } else if (p_value_age < 0.05) {
    p_display_age <- "P < 0.05"
  } else {
    p_display_age <- sprintf("P = %.3f", p_value_age)
  }
  
  # Save test results
  age_test_result <- data.frame(
    Comparison = "Age (65-90 vs 40-65)",
    Test = "Wilcoxon signed-rank test",
    Statistic = as.numeric(statistic_value_age),
    P_value = p_value_age,
    N_pairs = sum(valid_pairs_age),
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("Wilcoxon paired test result: %s (n=%d pairs)\n", p_display_age, sum(valid_pairs_age)))
} else {
  p_display_age <- "Insufficient data"
  age_test_result <- data.frame(
    Comparison = "Age (65-90 vs 40-65)",
    Test = "Wilcoxon signed-rank test",
    Statistic = NA,
    P_value = NA,
    N_pairs = sum(valid_pairs_age),
    stringsAsFactors = FALSE
  )
  cat("Warning: Insufficient paired data for testing\n")
}

# Create violin plot
p_violin_age <- ggplot(age_data_long, aes(x = Group, y = Correlation, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
  scale_fill_manual(values = c("40-65" = "#5CB85C", "65-90" = "#D9534F")) +
  labs(x = "Age Group", y = "PHI", 
       title = "Age Group Comparison") +
  annotate("text", 
           x = 1.5, y = Inf,
           label = p_display_age,
           hjust = 0.5, vjust = 1.2,
           size = 5) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(file.path(figs_dir, "violin_age_comparison.png"), 
       plot = p_violin_age, width = 5, height = 6, dpi = 300)
cat(sprintf("Saved: %s\n", file.path(figs_dir, "violin_age_comparison.png")))

# Save paired test results to CSV
test_results <- rbind(sex_test_result, age_test_result)
test_results_file <- file.path(figs_dir, "paired_test_results.csv")
write.csv(test_results, test_results_file, row.names = FALSE)
cat(sprintf("Paired test results saved to: %s\n", test_results_file))

# ============================================================================
# 5.5. Sex group scatter plot (Male vs Female)
# ============================================================================

cat("\nGenerating sex group scatter plot...\n")

# Prepare data (only use proteins with both Male and Female data)
valid_pairs_scatter <- !is.na(sex_data$Male) & !is.na(sex_data$Female)
if (sum(valid_pairs_scatter) >= 2) {
  scatter_data_sex <- data.frame(
    Male = sex_data$Male[valid_pairs_scatter],
    Female = sex_data$Female[valid_pairs_scatter]
  )
  
  # Calculate Pearson correlation coefficient and p-value
  cor_result_sex <- cor.test(scatter_data_sex$Male, scatter_data_sex$Female, method = "pearson")
  corr_coef_sex <- cor_result_sex$estimate
  p_value_sex <- cor_result_sex$p.value
  
  # Format p-value display
  if (p_value_sex < 0.001) {
    p_display_sex <- "P < 0.001"
  } else if (p_value_sex < 0.01) {
    p_display_sex <- "P < 0.01"
  } else if (p_value_sex < 0.05) {
    p_display_sex <- "P < 0.05"
  } else {
    p_display_sex <- sprintf("P = %.3f", p_value_sex)
  }
  
  # Create scatter plot
  p_scatter_sex <- ggplot(scatter_data_sex, aes(x = Male, y = Female)) +
    geom_point(alpha = 0.7, size = 2, color = "#2E86AB") +
    geom_smooth(method = "lm", se = TRUE, color = "#A23B72", linetype = "dashed", linewidth = 1) +
    labs(x = "PHI (Male)", y = "PHI (Female)") +
    annotate("text", 
             x = -Inf, y = Inf,
             label = sprintf("R = %.3f\n%s", corr_coef_sex, p_display_sex),
             hjust = -0.1, vjust = 1.5,
             size = 5) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 15),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  ggsave(file.path(figs_dir, "scatter_sex_male_vs_female.png"), 
         plot = p_scatter_sex, width = 5, height = 5, dpi = 300)
  cat(sprintf("Saved: %s\n", file.path(figs_dir, "scatter_sex_male_vs_female.png")))
  cat(sprintf("  Correlation coefficient: %.4f, P-value: %.4f (n=%d)\n", corr_coef_sex, p_value_sex, nrow(scatter_data_sex)))
} else {
  cat("Warning: Insufficient paired data to generate sex scatter plot\n")
}

# ============================================================================
# 5.6. Age group scatter plot (40-65 vs 65-90)
# ============================================================================

cat("\nGenerating age group scatter plot...\n")

# Prepare data (only use proteins with both 40-65 and 65-90 data)
valid_pairs_age_scatter <- !is.na(age_data$Age_40_65) & !is.na(age_data$Age_65_90)
if (sum(valid_pairs_age_scatter) >= 2) {
  scatter_data_age <- data.frame(
    Age_40_65 = age_data$Age_40_65[valid_pairs_age_scatter],
    Age_65_90 = age_data$Age_65_90[valid_pairs_age_scatter]
  )
  
  # Calculate Pearson correlation coefficient and p-value
  cor_result_age <- cor.test(scatter_data_age$Age_40_65, scatter_data_age$Age_65_90, method = "pearson")
  corr_coef_age <- cor_result_age$estimate
  p_value_age_scatter <- cor_result_age$p.value
  
  # Format p-value display
  if (p_value_age_scatter < 0.001) {
    p_display_age_scatter <- "P < 0.001"
  } else if (p_value_age_scatter < 0.01) {
    p_display_age_scatter <- "P < 0.01"
  } else if (p_value_age_scatter < 0.05) {
    p_display_age_scatter <- "P < 0.05"
  } else {
    p_display_age_scatter <- sprintf("P = %.3f", p_value_age_scatter)
  }
  
  # Create scatter plot
  p_scatter_age <- ggplot(scatter_data_age, aes(x = Age_40_65, y = Age_65_90)) +
    geom_point(alpha = 0.7, size = 2, color = "#2E86AB") +
    geom_smooth(method = "lm", se = TRUE, color = "#A23B72", linetype = "dashed", linewidth = 1) +
    labs(x = "PHI (Age 40-65)", y = "PHI (Age 65-90)") +
    annotate("text", 
             x = -Inf, y = Inf,
             label = sprintf("R = %.3f\n%s", corr_coef_age, p_display_age_scatter),
             hjust = -0.1, vjust = 1.5,
             size = 5) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 15),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  ggsave(file.path(figs_dir, "scatter_age_40-65_vs_65-90.png"), 
         plot = p_scatter_age, width = 5, height = 5, dpi = 300)
  cat(sprintf("Saved: %s\n", file.path(figs_dir, "scatter_age_40-65_vs_65-90.png")))
  cat(sprintf("  Correlation coefficient: %.4f, P-value: %.4f (n=%d)\n", corr_coef_age, p_value_age_scatter, nrow(scatter_data_age)))
} else {
  cat("Warning: Insufficient paired data to generate age scatter plot\n")
}

# ============================================================================
# 6. Sex difference bar chart (Male - Female)
# ============================================================================

cat("\nGenerating sex difference bar chart...\n")

# Calculate difference
df$diff_sex <- df[[sex_male_corr]] - df[[sex_female_corr]]

# Add gene symbols
df$GeneLabel <- ifelse(df[[entrez_col]] != "Null" & !is.na(df[[entrez_col]]) & df[[entrez_col]] != "",
                       df[[entrez_col]], 
                       NA)

# Sort by difference
df_sorted_sex <- df %>%
  arrange(diff_sex) %>%
  mutate(Index = row_number())

# Find top10 and bottom10 (exclude NA)
valid_diff_sex <- !is.na(df_sorted_sex$diff_sex)
df_valid_sex <- df_sorted_sex[valid_diff_sex, ]

if (nrow(df_valid_sex) > 0) {
  # Top 20 (largest positive values)
  top20_sex <- df_valid_sex %>%
    filter(diff_sex > 0) %>%
    arrange(desc(diff_sex)) %>%
    head(20) %>%
    filter(!is.na(GeneLabel))
  
  # Bottom 20 (smallest negative values)
  bottom20_sex <- df_valid_sex %>%
    filter(diff_sex < 0) %>%
    arrange(diff_sex) %>%
    head(20) %>%
    filter(!is.na(GeneLabel))
  
  # Calculate y-axis range
  y_max <- max(df_valid_sex$diff_sex, na.rm = TRUE)
  y_min <- min(df_valid_sex$diff_sex, na.rm = TRUE)
  y_range <- y_max - y_min
  
  # Calculate x-axis range
  x_max <- max(df_valid_sex$Index, na.rm = TRUE)
  x_min <- min(df_valid_sex$Index, na.rm = TRUE)
  x_range <- x_max - x_min
  x_center <- x_min + x_range * 0.5
  
  # top20 and bottom20 have been filtered for subsequent separate bar chart plotting
  
  # Create main bar chart (without gene labels and connecting lines)
  p_bar_sex <- ggplot(df_valid_sex, aes(x = Index, y = diff_sex)) +
    geom_bar(stat = "identity", width = 0.1, fill = "gray70", color = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    labs(x = "Protein Index (sorted by difference)", 
         y = "Difference (Male - Female)",
         title = "Sex Group Difference (Male - Female)") +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, hjust = 0.5),
      plot.margin = margin(20, 20, 20, 20)
    )
  
  ggsave(file.path(figs_dir, "barplot_sex_difference.png"), 
         plot = p_bar_sex, width = 8, height = 8, dpi = 300)
  cat(sprintf("Saved: %s\n", file.path(figs_dir, "barplot_sex_difference.png")))
  
  # Create separate dot plot for top 20 genes with sex difference > 0
  if (nrow(top20_sex) > 0) {
    top20_sex_plot <- top20_sex %>%
      arrange(desc(diff_sex)) %>%
      mutate(
        GeneLabel = factor(GeneLabel, levels = rev(GeneLabel)),
        x_pos = 1  # All points have the same x coordinate
      )
    
    p_top20_sex <- ggplot(top20_sex_plot, aes(x = x_pos, y = GeneLabel, color = diff_sex)) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_gradient(low = "lightpink", high = "darkred", guide = "none") +
      scale_x_continuous(limits = c(0.99, 1.01), expand = c(0, 0)) +
      theme_bw() +
      theme(
        axis.text.y = element_text(size = 9),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        plot.title = element_blank(),
        plot.margin = margin(5, 2, 5, 5, unit = "mm"),  # Fixed right margin (small), left margin auto-expands
        panel.spacing = unit(0, "mm")
      ) +
      coord_fixed(ratio = 0.01)  # Fixed aspect ratio to make dot plot very narrow
    
    ggsave(file.path(figs_dir, "barplot_sex_top20.png"), 
           plot = p_top20_sex, width = NA, height = 8 * 0.3, dpi = 300, units = "in", 
           limitsize = FALSE)
    cat(sprintf("Saved: %s\n", file.path(figs_dir, "barplot_sex_top20.png")))
  }
  
  # Create separate dot plot for bottom 20 genes with sex difference < 0
  if (nrow(bottom20_sex) > 0) {
    bottom20_sex_plot <- bottom20_sex %>%
      arrange(diff_sex) %>%
      mutate(
        GeneLabel = factor(GeneLabel, levels = rev(GeneLabel)),
        x_pos = 1  # All points have the same x coordinate
      )
    
    p_bottom20_sex <- ggplot(bottom20_sex_plot, aes(x = x_pos, y = GeneLabel, color = diff_sex)) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_gradient(low = "darkblue", high = "lightblue", guide = "none") +
      scale_x_continuous(limits = c(0.99, 1.01), expand = c(0, 0)) +
      theme_bw() +
      theme(
        axis.text.y = element_text(size = 9),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        plot.title = element_blank(),
        plot.margin = margin(5, 2, 5, 5, unit = "mm"),  # Fixed right margin (small), left margin auto-expands
        panel.spacing = unit(0, "mm")
      ) +
      coord_fixed(ratio = 0.01)  # Fixed aspect ratio to make dot plot very narrow
    
    ggsave(file.path(figs_dir, "barplot_sex_bottom20.png"), 
           plot = p_bottom20_sex, width = NA, height = 8 * 0.3, dpi = 300, units = "in", 
           limitsize = FALSE)
    cat(sprintf("Saved: %s\n", file.path(figs_dir, "barplot_sex_bottom20.png")))
  }
  
  cat(sprintf("  Top 20 (largest positive): %d genes\n", nrow(top20_sex)))
  cat(sprintf("  Bottom 20 (smallest negative): %d genes\n", nrow(bottom20_sex)))
} else {
  cat("Warning: No valid sex difference data\n")
}

# ============================================================================
# 7. Age difference bar chart (65-90 - 40-65)
# ============================================================================

cat("\nGenerating age difference bar chart...\n")

# Calculate difference
df$diff_age <- df[[age_65_90_corr]] - df[[age_40_65_corr]]

# Sort by difference
df_sorted_age <- df %>%
  arrange(diff_age) %>%
  mutate(Index = row_number())

# Find top10 and bottom10 (exclude NA)
valid_diff_age <- !is.na(df_sorted_age$diff_age)
df_valid_age <- df_sorted_age[valid_diff_age, ]

if (nrow(df_valid_age) > 0) {
  # Top 20 (largest positive values)
  top20_age <- df_valid_age %>%
    filter(diff_age > 0) %>%
    arrange(desc(diff_age)) %>%
    head(20) %>%
    filter(!is.na(GeneLabel))
  
  # Bottom 20 (smallest negative values)
  bottom20_age <- df_valid_age %>%
    filter(diff_age < 0) %>%
    arrange(diff_age) %>%
    head(20) %>%
    filter(!is.na(GeneLabel))
  
  # Calculate y-axis range
  y_max_age <- max(df_valid_age$diff_age, na.rm = TRUE)
  y_min_age <- min(df_valid_age$diff_age, na.rm = TRUE)
  y_range_age <- y_max_age - y_min_age
  
  # Calculate x-axis range
  x_max_age <- max(df_valid_age$Index, na.rm = TRUE)
  x_min_age <- min(df_valid_age$Index, na.rm = TRUE)
  x_range_age <- x_max_age - x_min_age
  x_center_age <- x_min_age + x_range_age * 0.5
  
  # top20 and bottom20 have been filtered for subsequent separate bar chart plotting
  
  # Create main bar chart (without gene labels and connecting lines)
  p_bar_age <- ggplot(df_valid_age, aes(x = Index, y = diff_age)) +
    geom_bar(stat = "identity", width = 0.1, fill = "gray70", color = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    labs(x = "Protein Index (sorted by difference)", 
         y = "Difference (65-90 - 40-65)",
         title = "Age Group Difference (65-90 - 40-65)") +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, hjust = 0.5),
      plot.margin = margin(20, 20, 20, 20)
    )
  
  ggsave(file.path(figs_dir, "barplot_age_difference.png"), 
         plot = p_bar_age, width = 8, height = 8, dpi = 300)
  cat(sprintf("Saved: %s\n", file.path(figs_dir, "barplot_age_difference.png")))
  
  # Create separate dot plot for top 20 genes with age difference > 0
  if (nrow(top20_age) > 0) {
    top20_age_plot <- top20_age %>%
      arrange(desc(diff_age)) %>%
      mutate(
        GeneLabel = factor(GeneLabel, levels = rev(GeneLabel)),
        x_pos = 1  # All points have the same x coordinate
      )
    
    p_top20_age <- ggplot(top20_age_plot, aes(x = x_pos, y = GeneLabel, color = diff_age)) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_gradient(low = "lightpink", high = "darkred", guide = "none") +
      scale_x_continuous(limits = c(0.99, 1.01), expand = c(0, 0)) +
      theme_bw() +
      theme(
        axis.text.y = element_text(size = 9),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        plot.title = element_blank(),
        plot.margin = margin(5, 2, 5, 5, unit = "mm"),  # Fixed right margin (small), left margin auto-expands
        panel.spacing = unit(0, "mm")
      ) +
      coord_fixed(ratio = 0.01)  # Fixed aspect ratio to make dot plot very narrow
    
    ggsave(file.path(figs_dir, "barplot_age_top20.png"), 
           plot = p_top20_age, width = NA, height = 8 * 0.3, dpi = 300, units = "in", 
           limitsize = FALSE)
    cat(sprintf("Saved: %s\n", file.path(figs_dir, "barplot_age_top20.png")))
  }
  
  # Create separate dot plot for bottom 20 genes with age difference < 0
  if (nrow(bottom20_age) > 0) {
    bottom20_age_plot <- bottom20_age %>%
      arrange(diff_age) %>%
      mutate(
        GeneLabel = factor(GeneLabel, levels = rev(GeneLabel)),
        x_pos = 1  # All points have the same x coordinate
      )
    
    p_bottom20_age <- ggplot(bottom20_age_plot, aes(x = x_pos, y = GeneLabel, color = diff_age)) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_gradient(low = "darkblue", high = "lightblue", guide = "none") +
      scale_x_continuous(limits = c(0.99, 1.01), expand = c(0, 0)) +
      theme_bw() +
      theme(
        axis.text.y = element_text(size = 9),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        plot.title = element_blank(),
        plot.margin = margin(5, 2, 5, 5, unit = "mm"),  # Fixed right margin (small), left margin auto-expands
        panel.spacing = unit(0, "mm")
      ) +
      coord_fixed(ratio = 0.01)  # Fixed aspect ratio to make dot plot very narrow
    
    ggsave(file.path(figs_dir, "barplot_age_bottom20.png"), 
           plot = p_bottom20_age, width = NA, height = 8 * 0.3, dpi = 300, units = "in", 
           limitsize = FALSE)
    cat(sprintf("Saved: %s\n", file.path(figs_dir, "barplot_age_bottom20.png")))
  }
  
  cat(sprintf("  Top 20 (largest positive): %d genes\n", nrow(top20_age)))
  cat(sprintf("  Bottom 20 (smallest negative): %d genes\n", nrow(bottom20_age)))
} else {
  cat("Warning: No valid age difference data\n")
}

# ============================================================================
# 8. Output top 50 proteins for four difference types to table
# ============================================================================

cat("\nGenerating top 50 protein table...\n")

# Prepare result list
top50_results <- data.frame(
  DifferenceType = character(),
  ProteinList = character(),
  stringsAsFactors = FALSE
)

# 1. Top 50 with sex difference > 0
if (exists("df_valid_sex") && nrow(df_valid_sex) > 0) {
  top50_sex_pos <- df_valid_sex %>%
    filter(diff_sex > 0, !is.na(GeneLabel)) %>%
    arrange(desc(diff_sex)) %>%
    head(50)
  
  if (nrow(top50_sex_pos) > 0) {
    protein_list <- paste(top50_sex_pos$GeneLabel, collapse = ", ")
    top50_results <- rbind(top50_results, data.frame(
      DifferenceType = "Sex Difference > 0 (Male - Female)",
      ProteinList = protein_list,
      stringsAsFactors = FALSE
    ))
    cat(sprintf("  Top 50 with sex difference > 0: %d proteins\n", nrow(top50_sex_pos)))
  }
}

# 2. Top 50 with sex difference < 0
if (exists("df_valid_sex") && nrow(df_valid_sex) > 0) {
  top50_sex_neg <- df_valid_sex %>%
    filter(diff_sex < 0, !is.na(GeneLabel)) %>%
    arrange(diff_sex) %>%
    head(50)
  
  if (nrow(top50_sex_neg) > 0) {
    protein_list <- paste(top50_sex_neg$GeneLabel, collapse = ", ")
    top50_results <- rbind(top50_results, data.frame(
      DifferenceType = "Sex Difference < 0 (Male - Female)",
      ProteinList = protein_list,
      stringsAsFactors = FALSE
    ))
    cat(sprintf("  Top 50 with sex difference < 0: %d proteins\n", nrow(top50_sex_neg)))
  }
}

# 3. Top 50 with age difference > 0
if (exists("df_valid_age") && nrow(df_valid_age) > 0) {
  top50_age_pos <- df_valid_age %>%
    filter(diff_age > 0, !is.na(GeneLabel)) %>%
    arrange(desc(diff_age)) %>%
    head(50)
  
  if (nrow(top50_age_pos) > 0) {
    protein_list <- paste(top50_age_pos$GeneLabel, collapse = ", ")
    top50_results <- rbind(top50_results, data.frame(
      DifferenceType = "Age Difference > 0 (65-90 - 40-65)",
      ProteinList = protein_list,
      stringsAsFactors = FALSE
    ))
    cat(sprintf("  Top 50 with age difference > 0: %d proteins\n", nrow(top50_age_pos)))
  }
}

# 4. Top 50 with age difference < 0
if (exists("df_valid_age") && nrow(df_valid_age) > 0) {
  top50_age_neg <- df_valid_age %>%
    filter(diff_age < 0, !is.na(GeneLabel)) %>%
    arrange(diff_age) %>%
    head(50)
  
  if (nrow(top50_age_neg) > 0) {
    protein_list <- paste(top50_age_neg$GeneLabel, collapse = ", ")
    top50_results <- rbind(top50_results, data.frame(
      DifferenceType = "Age Difference < 0 (65-90 - 40-65)",
      ProteinList = protein_list,
      stringsAsFactors = FALSE
    ))
    cat(sprintf("  Top 50 with age difference < 0: %d proteins\n", nrow(top50_age_neg)))
  }
}

# Save results
if (nrow(top50_results) > 0) {
  top50_file <- file.path(figs_dir, "top50_proteins_by_difference.csv")
  write.csv(top50_results, top50_file, row.names = FALSE, quote = TRUE)
  cat(sprintf("Top 50 protein table saved to: %s\n", top50_file))
} else {
  cat("Warning: No valid top 50 protein data found\n")
}

cat("\nAll visualizations completed!\n")
cat(sprintf("All figures saved to: %s\n", figs_dir))
