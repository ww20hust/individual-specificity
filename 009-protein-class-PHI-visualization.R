#!/usr/bin/env Rscript
# 009-protein-class-PHI-visualization.R
# Analyze protein_class distribution and PHI 12-22 relationship
# Generate pie chart, violin plot, and stacked bar chart

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggsignif)
library(cowplot)

# ============================================================================
# 1. Configuration
# ============================================================================

# Conservation data file path
conservation_file <- "/home/ww/Project/AnzhenData/figure3-gtex或者tony不同类别蛋白的保守度/conservation_final.csv"

# PHI results file path
phi_file <- "./results/001-PHI-different-visit-results.csv"

# Output directory
output_dir <- "./results/Organ-Secreted-Analysis"
figs_dir <- file.path(output_dir, "Figs")

# ============================================================================
# 2. Helper function to find column names
# ============================================================================

find_column <- function(df, pattern) {
  if (pattern %in% colnames(df)) {
    return(pattern)
  }
  x_pattern <- paste0("X", pattern)
  if (x_pattern %in% colnames(df)) {
    return(x_pattern)
  }
  matches <- grep(paste0("^", pattern, "$"), colnames(df), ignore.case = TRUE, value = TRUE)
  if (length(matches) > 0) {
    return(matches[1])
  }
  return(NULL)
}

# ============================================================================
# 3. Read conservation data
# ============================================================================

cat("=== Protein Class PHI Visualization ===\n\n")
cat("Reading conservation data...\n")

if (!file.exists(conservation_file)) {
  stop(sprintf("Error: File not found: %s\n", conservation_file))
}

df_conservation <- read.csv(conservation_file, stringsAsFactors = FALSE, check.names = FALSE)

cat(sprintf("Conservation data loaded: %d rows, %d columns\n", nrow(df_conservation), ncol(df_conservation)))
cat("Column names:\n")
print(colnames(df_conservation))

# Check required columns
if (!"Protein_seq" %in% colnames(df_conservation)) {
  stop("Error: Protein_seq column not found\n")
}
if (!"protein_class" %in% colnames(df_conservation)) {
  stop("Error: protein_class column not found\n")
}

# Check all unique protein_class values
cat("All unique protein_class values:\n")
print(unique(df_conservation$protein_class))
cat("Protein class distribution (including NA):\n")
print(table(df_conservation$protein_class, useNA = "ifany"))

# Extract required columns and clean protein_class
# Replace NA with "extracellular\nmembrane" (with line break)
df_conservation_clean <- df_conservation %>%
  select(Protein_seq, protein_class) %>%
  filter(!is.na(Protein_seq) & Protein_seq != "") %>%
  mutate(
    protein_class = trimws(protein_class),
    protein_class = ifelse(is.na(protein_class) | protein_class == "", "extracellular\nmembrane", protein_class)
  )

cat(sprintf("Cleaned conservation data: %d rows\n", nrow(df_conservation_clean)))
cat("Protein class distribution after cleaning:\n")
print(table(df_conservation_clean$protein_class, useNA = "ifany"))

# ============================================================================
# 4. Read PHI data
# ============================================================================

cat("\nReading PHI data...\n")

if (!file.exists(phi_file)) {
  stop(sprintf("Error: File not found: %s\n", phi_file))
}

df_phi <- read.csv(phi_file, stringsAsFactors = FALSE, check.names = FALSE)

cat(sprintf("PHI data loaded: %d rows, %d columns\n", nrow(df_phi), ncol(df_phi)))

# Find required columns
aptname_col <- find_column(df_phi, "AptName")
phi_col <- find_column(df_phi, "12_22_corr")

if (is.null(aptname_col)) {
  stop("Error: AptName column not found\n")
}
if (is.null(phi_col)) {
  stop("Error: 12_22_corr column not found\n")
}

cat(sprintf("Using columns: AptName=%s, PHI=%s\n", aptname_col, phi_col))

# Extract required columns
df_phi_clean <- df_phi %>%
  select(AptName = !!aptname_col, PHI = !!phi_col) %>%
  filter(!is.na(AptName) & AptName != "",
         !is.na(PHI))

cat(sprintf("Cleaned PHI data: %d rows\n", nrow(df_phi_clean)))

# ============================================================================
# 5. Merge data
# ============================================================================

cat("\nMerging data...\n")

df_merged <- df_conservation_clean %>%
  inner_join(df_phi_clean, by = c("Protein_seq" = "AptName"))

cat(sprintf("Merged data: %d rows\n", nrow(df_merged)))
cat("Protein class distribution in merged data (before cleaning):\n")
print(table(df_merged$protein_class, useNA = "ifany"))

# Clean protein_class in merged data
# Replace NA with "extracellular\nmembrane" (with line break)
df_merged$protein_class <- trimws(df_merged$protein_class)
df_merged$protein_class <- ifelse(is.na(df_merged$protein_class) | df_merged$protein_class == "", 
                                   "extracellular\nmembrane", 
                                   df_merged$protein_class)

cat("Protein class distribution in merged data (after cleaning):\n")
print(table(df_merged$protein_class, useNA = "ifany"))

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", output_dir))
}
if (!dir.exists(figs_dir)) {
  dir.create(figs_dir, recursive = TRUE)
  cat(sprintf("Created figures directory: %s\n", figs_dir))
}

# Save merged data
merged_file <- file.path(output_dir, "protein_class_PHI_merged.csv")
write.csv(df_merged, merged_file, row.names = FALSE)
cat(sprintf("Merged data saved to: %s\n", merged_file))


# ============================================================================
# 6. Pie chart: protein class distribution（修正版）
# ============================================================================

cat("\nGenerating pie chart...\n")

pie_data <- df_conservation_clean %>%
  count(protein_class) %>%
  mutate(percentage = n / sum(n) * 100)

pie_data$protein_class <- factor(
  pie_data$protein_class,
  levels = c("intracellular", "extracellular\nmembrane", "secreted")
)

color_palette <- c(
  "intracellular" = "#B4C7E7",
  "extracellular\nmembrane" = "#FFE699",
  "secreted" = "#C5E0B4"
)

p_pie <- ggplot(pie_data, aes(x = "", y = percentage, fill = protein_class)) +
  geom_bar(
    stat = "identity",
    width = 1,
    color = "white"
  ) +
  coord_polar(theta = "y") +
  geom_text(
    aes(label = sprintf("%.1f%%", percentage)),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    fontface = "bold"
  ) +
  scale_fill_manual(values = color_palette) +
  theme_void() +
  theme(legend.position = "none")

ggsave(
  file.path(figs_dir, "protein_class_pie_chart.png"),
  p_pie,
  width = 2.5,
  height = 2.5,
  dpi = 300
)

# legend
p_legend <- ggplot(pie_data, aes(x = protein_class, fill = protein_class)) +
  geom_bar() +
  scale_fill_manual(values = color_palette, name = "Protein class") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

legend_only <- cowplot::get_legend(p_legend)

ggsave(
  file.path(figs_dir, "protein_class_legend.png"),
  legend_only,
  width = 1.5,
  height = 1.5,
  dpi = 300
)

# ============================================================================
# 7. Violin plot: PHI distribution by protein class
# ============================================================================

cat("\nGenerating violin plot...\n")

# Clean and set factor levels for protein_class (intracellular left, secreted right)
df_merged$protein_class <- trimws(df_merged$protein_class)
df_merged$protein_class <- ifelse(is.na(df_merged$protein_class) | df_merged$protein_class == "", 
                                   "extracellular\nmembrane", 
                                   df_merged$protein_class)
df_merged$protein_class <- factor(df_merged$protein_class,
                                   levels = c("intracellular", "extracellular\nmembrane", "secreted"))

cat("Protein class distribution in merged data (after factor setting):\n")
print(table(df_merged$protein_class, useNA = "ifany"))

# Scientific color palette (light colors)
color_palette <- c(
  "intracellular" = "#B4C7E7",                    # Light blue
  "extracellular\nmembrane" = "#FFE699",         # Light yellow/orange
  "secreted" = "#C5E0B4"                          # Light green
)

# Statistical tests between groups
cat("\nPerforming statistical tests between protein classes...\n")

# Extract PHI values for each group
group_intra <- df_merged$PHI[df_merged$protein_class == "intracellular"]
group_extra <- df_merged$PHI[df_merged$protein_class == "extracellular\nmembrane"]
group_secreted <- df_merged$PHI[df_merged$protein_class == "secreted"]

# Perform pairwise comparisons
comparisons <- list(
  c("intracellular", "extracellular\nmembrane"),
  c("extracellular\nmembrane", "secreted"),
  c("intracellular", "secreted")
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
  group1_data <- df_merged$PHI[df_merged$protein_class == comp[1]]
  group2_data <- df_merged$PHI[df_merged$protein_class == comp[2]]
  
  if (length(group1_data) >= 2 && length(group2_data) >= 2) {
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
test_results_file <- file.path(output_dir, "protein_class_statistical_tests.csv")
write.csv(test_results, test_results_file, row.names = FALSE)
cat(sprintf("Test results saved to: %s\n", test_results_file))

# Calculate y-axis range for significance annotations
y_max <- max(df_merged$PHI, na.rm = TRUE)
y_range <- max(df_merged$PHI, na.rm = TRUE) - min(df_merged$PHI, na.rm = TRUE)

# Create violin plot
p_violin <- ggplot(df_merged, aes(x = protein_class, y = PHI, fill = protein_class)) +
  geom_violin(trim = FALSE, alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.8,
               outlier.alpha = 0.5, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = color_palette) +
  labs(
    x = "",
    y = "PHI"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    plot.margin = margin(15, 15, 15, 15)
  ) +
  geom_signif(
    comparisons = comparisons,
    map_signif_level = function(p) {
      if (p < 0.001) return("***")
      if (p < 0.01) return("**")
      if (p < 0.05) return("*")
      return("ns")
    },
    textsize = 3.5,
    y_position = y_max + y_range * c(0.15, 0.25, 0.35),
    step_increase = 0
  )

# Save violin plot
violin_file <- file.path(figs_dir, "protein_class_phi_violin.png")
ggsave(violin_file, plot = p_violin, width = 5, height = 4, dpi = 300)
cat(sprintf("Violin plot saved to: %s\n", violin_file))

# ============================================================================
# 8. Stacked bar chart: protein class proportion in PHI bins（修正版）
# ============================================================================

cat("\nGenerating stacked bar chart...\n")

df_merged_binned <- df_merged %>%
  mutate(
    PHI_bin = case_when(
      PHI >= 0 & PHI < 0.25 ~ "0-0.25",
      PHI >= 0.25 & PHI < 0.5 ~ "0.25-0.5",
      PHI >= 0.5 & PHI < 0.75 ~ "0.5-0.75",
      PHI >= 0.75 & PHI <= 1 ~ "0.75-1"
    )
  ) %>%
  filter(!is.na(PHI_bin))

df_merged_binned$PHI_bin <- factor(
  df_merged_binned$PHI_bin,
  levels = c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1")
)

df_merged_binned$protein_class <- factor(
  df_merged_binned$protein_class,
  levels = c("intracellular", "extracellular\nmembrane", "secreted")
)

stacked_data <- df_merged_binned %>%
  count(PHI_bin, protein_class) %>%
  group_by(PHI_bin) %>%
  mutate(proportion = n / sum(n) * 100) %>%
  ungroup()

p_stacked <- ggplot(
  stacked_data,
  aes(x = PHI_bin, y = proportion, fill = protein_class)
) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(
    aes(label = sprintf("%.1f%%", proportion)),
    position = position_stack(vjust = 0.5),
    size = 3,
    fontface = "bold"
  ) +
  scale_fill_manual(values = color_palette) +
  labs(
    x = "PHI bin",
    y = "Proportion (%)"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(figs_dir, "protein_class_phi_stacked_bar.png"),
  p_stacked,
  width = 5,
  height = 4,
  dpi = 300
)

cat("\nCompleted!\n")


# ============================================================================
# 9. Summary statistics
# ============================================================================

cat("\n=== Summary Statistics ===\n")
cat(sprintf("Total proteins in conservation data: %d\n", nrow(df_conservation_clean)))
cat(sprintf("Total proteins in PHI data: %d\n", nrow(df_phi_clean)))
cat(sprintf("Matched proteins: %d\n", nrow(df_merged)))

cat("\nProtein class distribution:\n")
print(table(df_merged$protein_class))

cat("\nPHI statistics by protein class:\n")
phi_stats <- df_merged %>%
  group_by(protein_class) %>%
  summarise(
    n = n(),
    mean = mean(PHI, na.rm = TRUE),
    median = median(PHI, na.rm = TRUE),
    sd = sd(PHI, na.rm = TRUE),
    min = min(PHI, na.rm = TRUE),
    max = max(PHI, na.rm = TRUE),
    .groups = "drop"
  )
print(phi_stats)

cat("\nPHI bin distribution:\n")
print(table(df_merged_binned$PHI_bin))

cat("\nProtein class distribution in each PHI bin:\n")
print(table(df_merged_binned$PHI_bin, df_merged_binned$protein_class))

cat("\nCompleted!\n")

# ============================================================================
# 10. Tissue enriched analysis by PHI groups
# ============================================================================

cat("\n=== Tissue Enriched Analysis ===\n")

# Read tissue enriched data
tissue_file <- "/home/ww/Project/AnzhenData/figure3-gtex或者tony不同类别蛋白的保守度/conservation_Tissue_enriched.csv"

cat("Reading tissue enriched data...\n")

if (!file.exists(tissue_file)) {
  stop(sprintf("Error: File not found: %s\n", tissue_file))
}

df_tissue <- read.csv(tissue_file, stringsAsFactors = FALSE, check.names = FALSE)

cat(sprintf("Tissue enriched data loaded: %d rows, %d columns\n", nrow(df_tissue), ncol(df_tissue)))
cat("Column names:\n")
print(colnames(df_tissue))

# Check required columns
if (!"Protein_seq" %in% colnames(df_tissue)) {
  stop("Error: Protein_seq column not found\n")
}
if (!"tony_organ" %in% colnames(df_tissue)) {
  stop("Error: tony_organ column not found\n")
}

# Extract required columns
df_tissue_clean <- df_tissue %>%
  select(Protein_seq, tony_organ) %>%
  filter(!is.na(Protein_seq) & Protein_seq != "",
         !is.na(tony_organ) & tony_organ != "")

cat(sprintf("Cleaned tissue data: %d rows\n", nrow(df_tissue_clean)))
cat("Organ distribution:\n")
print(table(df_tissue_clean$tony_organ))

# First, create PHI groups for all PHI data
cat("\nCreating PHI groups for all proteins...\n")
df_phi_grouped <- df_phi_clean %>%
  mutate(
    PHI_group = ifelse(PHI > 0.5, "PHI > 0.5", "PHI <= 0.5")
  )

cat("PHI group distribution (all proteins):\n")
print(table(df_phi_grouped$PHI_group))

# Merge with tissue data (left join to keep all PHI proteins)
cat("\nMerging tissue data with PHI data...\n")
df_tissue_merged <- df_phi_grouped %>%
  left_join(df_tissue_clean, by = c("AptName" = "Protein_seq"))

cat(sprintf("Merged tissue-PHI data: %d rows\n", nrow(df_tissue_merged)))
cat("Proteins with organ information: %d\n", sum(!is.na(df_tissue_merged$tony_organ)))
cat("Proteins without organ information: %d\n", sum(is.na(df_tissue_merged$tony_organ)))

# Count proteins by organ and PHI group (calculate percentage)
organ_counts <- df_tissue_merged %>%
  filter(!is.na(tony_organ)) %>%
  group_by(PHI_group, tony_organ) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(PHI_group) %>%
  mutate(
    total_in_group = sum(count),
    percentage = count / total_in_group * 100
  ) %>%
  ungroup() %>%
  select(PHI_group, tony_organ, count, percentage) %>%
  group_by(PHI_group) %>%
  arrange(desc(percentage)) %>%
  ungroup()

cat("\nOrgan counts by PHI group:\n")
print(organ_counts)

# Create bar charts for each PHI group
for (phi_group in c("PHI > 0.5", "PHI <= 0.5")) {
  cat(sprintf("\nGenerating bar chart for %s...\n", phi_group))
  
  # Filter data for this PHI group
  group_data <- organ_counts %>%
    filter(PHI_group == phi_group) %>%
    arrange(desc(count))
  
  # Set factor levels to maintain order (top to bottom)
  group_data$tony_organ <- factor(group_data$tony_organ, 
                                   levels = rev(group_data$tony_organ))
  
  # Create bar chart with percentage
  p_bar <- ggplot(group_data, aes(x = tony_organ, y = percentage)) +
    geom_bar(stat = "identity", fill = "#4472C4", width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", percentage)),
              hjust = -0.1,
              size = 3,
              fontface = "bold") +
    coord_flip() +
    labs(
      x = "Organ",
      y = "Percentage (%)"
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 11, color = "black"),
      axis.title = element_text(size = 13, color = "black"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(15, 15, 15, 15)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  
  # Save bar chart
  bar_file <- file.path(figs_dir, sprintf("organ_enriched_%s.png", 
                                           gsub(" ", "_", phi_group)))
  ggsave(bar_file, plot = p_bar, width = 6, height = max(4, nrow(group_data) * 0.3), 
         dpi = 300, limitsize = FALSE)
  cat(sprintf("Bar chart saved to: %s\n", bar_file))
}

# Save organ counts data
organ_counts_file <- file.path(output_dir, "organ_enriched_counts_by_PHI.csv")
write.csv(organ_counts, organ_counts_file, row.names = FALSE)
cat(sprintf("\nOrgan counts data saved to: %s\n", organ_counts_file))

cat("\nTissue enriched analysis completed!\n")

# ============================================================================
# 11. Organ-specific protein class analysis (PHI > 0.75)
# ============================================================================

cat("\n=== Organ-Specific Protein Class Analysis (PHI > 0.75) ===\n")

# Define organs to analyze
target_organs <- c("Brain", "Liver", "Immune")

# Get protein_class information from merged data
# We need to merge with conservation data to get protein_class
cat("Preparing data for organ-specific analysis...\n")

# Merge PHI data with conservation data to get protein_class
df_phi_with_class <- df_phi_clean %>%
  left_join(df_conservation_clean, by = c("AptName" = "Protein_seq")) %>%
  filter(!is.na(protein_class))

# Clean protein_class (replace NA with "extracellular\nmembrane")
df_phi_with_class$protein_class <- trimws(df_phi_with_class$protein_class)
df_phi_with_class$protein_class <- ifelse(
  is.na(df_phi_with_class$protein_class) | df_phi_with_class$protein_class == "",
  "extracellular\nmembrane",
  df_phi_with_class$protein_class
)

# Set factor levels
df_phi_with_class$protein_class <- factor(
  df_phi_with_class$protein_class,
  levels = c("intracellular", "extracellular\nmembrane", "secreted")
)

# Merge with tissue data
df_organ_analysis <- df_phi_with_class %>%
  left_join(df_tissue_clean, by = c("AptName" = "Protein_seq")) %>%
  filter(
    PHI > 0.75,
    !is.na(tony_organ),
    tony_organ %in% target_organs
  )

cat(sprintf("Proteins with PHI > 0.75 in target organs: %d\n", nrow(df_organ_analysis)))

# Scientific color palette (same as before)
color_palette <- c(
  "intracellular" = "#B4C7E7",
  "extracellular\nmembrane" = "#FFE699",
  "secreted" = "#C5E0B4"
)

# Analyze each organ
for (organ in target_organs) {
  cat(sprintf("\nAnalyzing %s...\n", organ))
  
  # Filter data for this organ
  organ_data <- df_organ_analysis %>%
    filter(tony_organ == organ)
  
  if (nrow(organ_data) == 0) {
    cat(sprintf("No data found for %s, skipping...\n", organ))
    next
  }
  
  cat(sprintf("Total proteins in %s (PHI > 0.75): %d\n", organ, nrow(organ_data)))
  
  # Count protein class distribution
  class_counts <- organ_data %>%
    count(protein_class) %>%
    mutate(
      percentage = n / sum(n) * 100
    ) %>%
    arrange(match(protein_class, c("intracellular", "extracellular\nmembrane", "secreted")))
  
  cat("Protein class distribution:\n")
  print(class_counts)
  
  # Set factor levels to maintain order
  class_counts$protein_class <- factor(
    class_counts$protein_class,
    levels = c("intracellular", "extracellular\nmembrane", "secreted")
  )
  
  # Create bar chart
  p_bar <- ggplot(class_counts, aes(x = protein_class, y = percentage, fill = protein_class)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", percentage)),
              vjust = -0.5,
              size = 3.5,
              fontface = "bold") +
    scale_fill_manual(values = color_palette) +
    labs(
      x = "Protein class",
      y = "Percentage (%)"
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 11, color = "black"),
      axis.text.x = element_text(size = 11, color = "black"),
      axis.title = element_text(size = 13, color = "black"),
      legend.position = "none",
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(15, 15, 15, 15)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  
  # Save bar chart
  bar_file <- file.path(figs_dir, sprintf("organ_%s_PHI_high_protein_class.png", organ))
  ggsave(bar_file, plot = p_bar, width = 4, height = 4, dpi = 300)
  cat(sprintf("Bar chart saved to: %s\n", bar_file))
  
  # Save data
  data_file <- file.path(output_dir, sprintf("organ_%s_PHI_high_protein_class.csv", organ))
  write.csv(class_counts, data_file, row.names = FALSE)
  cat(sprintf("Data saved to: %s\n", data_file))
}

cat("\nOrgan-specific protein class analysis completed!\n")


