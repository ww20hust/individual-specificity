#!/usr/bin/env Rscript
# 002-PHI-visualization.R
# Visualization and statistical analysis of PHI correlation results

# Load required libraries
library(ggplot2)
library(dplyr)
library(gridExtra)

# ============================================================================
# 1. Load data
# ============================================================================

cat("Loading results from 001-PHI-different-visit.R...\n")

# Read the results file
results_path <- "./results/001-PHI-different-visit-results.csv"
df <- read.csv(results_path, stringsAsFactors = FALSE, check.names = FALSE)

cat(sprintf("Data loaded: %d rows, %d columns\n", nrow(df), ncol(df)))
cat("Column names:\n")
print(colnames(df))

# Function to find column name (handle R's automatic name modification)
find_column <- function(df, pattern) {
  # Try exact match first
  if (pattern %in% colnames(df)) {
    return(pattern)
  }
  # Try with X prefix (R adds X to column names starting with numbers)
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

# Define visit pairs
visit_pairs <- c("02_07", "02_12", "02_22", "07_12", "07_22", "12_22")

# Create output directory for figures
output_base <- "./results"
figs_dir <- file.path(output_base, "Figs")
if (!dir.exists(figs_dir)) {
  dir.create(figs_dir, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", figs_dir))
}

# ============================================================================
# 2. Generate distribution plots for each visit pair
# ============================================================================

cat("\nGenerating distribution plots for each visit pair...\n")

for (pair in visit_pairs) {
  corr_col <- paste0(pair, "_corr")
  actual_col <- find_column(df, corr_col)
  
  if (is.null(actual_col)) {
    cat(sprintf("Warning: Column %s not found, skipping...\n", corr_col))
    next
  }
  
  # Extract correlation values (remove NA)
  corr_values <- df[[actual_col]]
  corr_values <- corr_values[!is.na(corr_values)]
  
  if (length(corr_values) == 0) {
    cat(sprintf("Warning: No valid data for %s, skipping...\n", pair))
    next
  }
  
  # Create data frame for plotting
  plot_data <- data.frame(Correlation = corr_values)
  
  # Create histogram
  p <- ggplot(plot_data, aes(x = Correlation)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
    labs(x = "PHI", y = "Count") +
    theme_bw() +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 15),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  # Save plot (height smaller than width)
  output_file <- file.path(figs_dir, paste0("protein-r-distribution_", pair, ".png"))
  ggsave(output_file, plot = p, width = 4, height = 3, dpi = 300)
  cat(sprintf("Saved: %s\n", output_file))
}

# ============================================================================
# 3. Generate summary statistics table
# ============================================================================

cat("\nGenerating summary statistics table...\n")

summary_stats <- data.frame(
  VisitPair = character(),
  Pvalue_lt_005_count = integer(),
  Pvalue_lt_005_proportion = numeric(),
  Median_corr = numeric(),
  Mean_corr = numeric(),
  Variance_corr = numeric(),
  Min_corr = numeric(),
  Max_corr = numeric(),
  stringsAsFactors = FALSE
)

for (pair in visit_pairs) {
  corr_col <- paste0(pair, "_corr")
  pval_col <- paste0(pair, "_pval")
  
  actual_corr_col <- find_column(df, corr_col)
  actual_pval_col <- find_column(df, pval_col)
  
  if (is.null(actual_corr_col) || is.null(actual_pval_col)) {
    next
  }
  
  # Get valid data (both correlation and p-value available)
  valid_idx <- !is.na(df[[actual_corr_col]]) & !is.na(df[[actual_pval_col]])
  corr_values <- df[[actual_corr_col]][valid_idx]
  pval_values <- df[[actual_pval_col]][valid_idx]
  
  if (length(corr_values) == 0) {
    next
  }
  
  # Calculate statistics
  pval_lt_005_count <- sum(pval_values < 0.05, na.rm = TRUE)
  pval_lt_005_prop <- pval_lt_005_count / length(pval_values)
  median_corr <- median(corr_values, na.rm = TRUE)
  mean_corr <- mean(corr_values, na.rm = TRUE)
  var_corr <- var(corr_values, na.rm = TRUE)
  min_corr <- min(corr_values, na.rm = TRUE)
  max_corr <- max(corr_values, na.rm = TRUE)
  
  summary_stats <- rbind(summary_stats, data.frame(
    VisitPair = pair,
    Pvalue_lt_005_count = pval_lt_005_count,
    Pvalue_lt_005_proportion = round(pval_lt_005_prop, 4),
    Median_corr = round(median_corr, 4),
    Mean_corr = round(mean_corr, 4),
    Variance_corr = round(var_corr, 4),
    Min_corr = round(min_corr, 4),
    Max_corr = round(max_corr, 4),
    stringsAsFactors = FALSE
  ))
}

# Save summary statistics
summary_file <- file.path(figs_dir, "summary_statistics.csv")
write.csv(summary_stats, summary_file, row.names = FALSE)
cat(sprintf("Summary statistics saved to: %s\n", summary_file))

# ============================================================================
# 4. Generate scatter plots for specific pairs
# ============================================================================

cat("\nGenerating scatter plots...\n")

# Function to create scatter plot
create_scatter_plot <- function(df, pair1, pair2, output_file, find_col_func) {
  corr_col1 <- paste0(pair1, "_corr")
  corr_col2 <- paste0(pair2, "_corr")
  
  actual_col1 <- find_col_func(df, corr_col1)
  actual_col2 <- find_col_func(df, corr_col2)
  
  if (is.null(actual_col1) || is.null(actual_col2)) {
    cat(sprintf("Warning: Columns not found for %s vs %s\n", pair1, pair2))
    cat(sprintf("  Looking for: %s, %s\n", corr_col1, corr_col2))
    cat(sprintf("  Found: %s, %s\n", 
                ifelse(is.null(actual_col1), "NULL", actual_col1),
                ifelse(is.null(actual_col2), "NULL", actual_col2)))
    return(NULL)
  }
  
  # Get valid data
  valid_idx <- !is.na(df[[actual_col1]]) & !is.na(df[[actual_col2]])
  x_values <- df[[actual_col1]][valid_idx]
  y_values <- df[[actual_col2]][valid_idx]
  
  if (length(x_values) < 2) {
    cat(sprintf("Warning: Insufficient data for %s vs %s\n", pair1, pair2))
    return(NULL)
  }
  
  # Calculate correlation and p-value
  cor_result <- cor.test(x_values, y_values, method = "pearson")
  corr_coef <- cor_result$estimate
  p_value <- cor_result$p.value
  
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
  
  # Create data frame for plotting
  plot_data <- data.frame(
    X = x_values,
    Y = y_values
  )
  
  # Create scatter plot with scientific style colors
  p <- ggplot(plot_data, aes(x = X, y = Y)) +
    geom_point(alpha = 0.7, size = 2, color = "#2E86AB") +
    geom_smooth(method = "lm", se = TRUE, color = "#A23B72", linetype = "dashed", linewidth = 1) +
    labs(x = paste0("PHI: ", pair1), y = paste0("PHI: ", pair2)) +
    annotate("text", 
             x = -Inf, y = Inf,
             label = sprintf("R = %.3f\n%s", corr_coef, p_display),
             hjust = -0.1, vjust = 1.5,
             size = 5) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 15),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  # Save plot
  ggsave(output_file, plot = p, width = 4, height = 4, dpi = 300)
  cat(sprintf("Saved: %s\n", output_file))
  
  return(p)
}

# Plot 1: 02-07 vs 02-12
create_scatter_plot(df, "02_07", "02_12", 
                    file.path(figs_dir, "scatter_02_07_vs_02_12.png"),
                    find_column)

# Plot 2: 02-07 vs 02-22
create_scatter_plot(df, "02_07", "02_22", 
                    file.path(figs_dir, "scatter_02_07_vs_02_22.png"),
                    find_column)

# Plot 3: 07-12 vs 07-22
create_scatter_plot(df, "07_12", "07_22", 
                    file.path(figs_dir, "scatter_07_12_vs_07_22.png"),
                    find_column)

# ============================================================================
# 5. Generate circular bar chart for top 30 EntrezGeneSymbol based on 02-22
# ============================================================================

cat("\nGenerating circular bar chart for top 30 EntrezGeneSymbol...\n")

# Get 02-22 correlation data
corr_col_0222 <- "02_22_corr"
entrez_col <- "EntrezGeneSymbol"

actual_corr_col_0222 <- find_column(df, corr_col_0222)
actual_entrez_col <- find_column(df, entrez_col)

if (!is.null(actual_corr_col_0222) && !is.null(actual_entrez_col)) {
  # Filter out Null EntrezGeneSymbol and NA correlations
  valid_data <- df[df[[actual_entrez_col]] != "Null" & 
                   !is.na(df[[actual_corr_col_0222]]) & 
                   !is.na(df[[actual_entrez_col]]) & 
                   df[[actual_entrez_col]] != "", ]
  
  if (nrow(valid_data) > 0) {
    # For each EntrezGeneSymbol, keep the maximum correlation
    top_data <- valid_data %>%
      group_by(.data[[actual_entrez_col]]) %>%
      summarise(MaxCorr = max(.data[[actual_corr_col_0222]], na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(MaxCorr)) %>%
      head(30)
    
    if (nrow(top_data) > 0) {
      # Prepare data for circular bar chart
      max_corr_val <- max(top_data$MaxCorr, na.rm = TRUE)
      top_data <- top_data %>%
        mutate(
          id = row_number(),
          angle = 90 - 360 * (id - 0.5) / nrow(top_data),
          hjust = ifelse(angle < -90, 1, 0),
          angle_adj = ifelse(angle < -90, angle + 180, angle),
          label_y = MaxCorr + max_corr_val * 0.12
        )
      
      # Calculate ylim
      ylim_max <- max_corr_val * 1.3
      
      # Create circular bar chart with gradient colors
      p <- ggplot(top_data, aes(x = factor(id), y = MaxCorr, fill = MaxCorr)) +
        geom_bar(stat = "identity", width = 0.8) +
        scale_fill_gradient(low = "#E3F2FD", high = "#1976D2", 
                           name = "PHI",
                           guide = guide_colorbar(
                             title.position = "top",
                             title.hjust = 0.5,
                             barwidth = 0.5,
                             barheight = 8,
                             direction = "vertical"
                           )) +
        coord_polar(start = 0) +
        ylim(-max_corr_val * 0.3, ylim_max) +
        geom_text(aes(x = id, y = label_y,
                      label = .data[[actual_entrez_col]],
                      angle = angle_adj, hjust = hjust),
                  size = 4.5) +
        annotate("text", x = 0, y = -0.23, label = "PHI", 
                 size = 8) +
        theme_void() +
        theme(
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          plot.margin = margin(5, 5, 5, 5),
          legend.position = "right",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10)
        )
      
      # Save plot
      output_file <- file.path(figs_dir, "circular_bar_top30_entrez_02_22.png")
      ggsave(output_file, plot = p, width = 5, height = 5, dpi = 300)
      cat(sprintf("Saved: %s\n", output_file))
      
      # Also save the top 30 data
      top30_file <- file.path(figs_dir, "top30_entrez_02_22.csv")
      # Rename column back to original name for output
      output_data <- top_data
      colnames(output_data)[colnames(output_data) == actual_entrez_col] <- entrez_col
      write.csv(output_data[, c(entrez_col, "MaxCorr")], 
                top30_file, row.names = FALSE)
      cat(sprintf("Top 30 data saved to: %s\n", top30_file))
    } else {
      cat("Warning: No valid data for circular bar chart\n")
    }
  } else {
    cat("Warning: No valid EntrezGeneSymbol data found\n")
  }
} else {
  cat("Warning: Required columns not found for circular bar chart\n")
}

cat("\nAll visualizations completed!\n")
cat(sprintf("All figures saved to: %s\n", figs_dir))
