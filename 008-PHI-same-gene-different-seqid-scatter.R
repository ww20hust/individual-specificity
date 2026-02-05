#!/usr/bin/env Rscript
# 008-PHI-same-gene-different-seqid-scatter.R
# Identify genes with the same name but different seqids
# Generate pairwise PHI combinations for all seqids
# Create scatter plots for all genes and all combinations (using 12-22 PHI)

# Load required libraries
library(ggplot2)
library(dplyr)

# ============================================================================
# 1. Configuration parameters
# ============================================================================

# PHI results file path
phi_file <- "./results/001-PHI-different-visit-results.csv"

# Output directory
output_base <- "./results"
output_dir <- file.path(output_base, "same-gene-different-seqid-PHI")

# ============================================================================
# 2. Read PHI results data
# ============================================================================

cat("=== PHI Scatter Plot Analysis for Same Gene Different Seqid ===\n\n")
cat("Reading PHI results data...\n")

if (!file.exists(phi_file)) {
  stop(sprintf("Error: File not found: %s\n", phi_file))
}

df_phi <- read.csv(phi_file, stringsAsFactors = FALSE, check.names = FALSE)

cat(sprintf("PHI data loaded: %d rows, %d columns\n", nrow(df_phi), ncol(df_phi)))
cat("Column names:\n")
print(colnames(df_phi))

# Function: Find column name (handle R's automatic X prefix addition)
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

# Find required columns
entrez_col <- find_column(df_phi, "EntrezGeneSymbol")
aptname_col <- find_column(df_phi, "AptName")
phi_col <- find_column(df_phi, "12_22_corr")

if (is.null(entrez_col)) {
  stop("Error: EntrezGeneSymbol column not found\n")
}
if (is.null(aptname_col)) {
  stop("Error: AptName column not found\n")
}
if (is.null(phi_col)) {
  stop("Error: 12_22_corr column not found\n")
}

cat(sprintf("Using columns: EntrezGeneSymbol=%s, AptName=%s, PHI=%s\n", entrez_col, aptname_col, phi_col))

# ============================================================================
# 3. Data preprocessing
# ============================================================================

cat("\nPreprocessing data...\n")

# Extract required columns and filter
df_clean <- df_phi %>%
  select(
    Gene = !!entrez_col,
    seqid = !!aptname_col,
    PHI = !!phi_col
  ) %>%
  filter(
    !is.na(Gene) & Gene != "" & Gene != "Null",
    !is.na(seqid) & seqid != "",
    !is.na(PHI)
  )

cat(sprintf("Cleaned data: %d rows\n", nrow(df_clean)))
cat(sprintf("Unique gene count: %d\n", length(unique(df_clean$Gene))))
cat(sprintf("Unique seqid count: %d\n", length(unique(df_clean$seqid))))

# ============================================================================
# 4. Identify genes with same name but different seqids
# ============================================================================

cat("\nIdentifying genes with same name but different seqids...\n")

# Find how many different seqids each gene has
gene_seqid_count <- df_clean %>%
  group_by(Gene) %>%
  summarise(
    n_seqid = n_distinct(seqid),
    seqids = paste(unique(seqid), collapse = "; "),
    .groups = "drop"
  ) %>%
  filter(n_seqid > 1)  # Keep only genes with multiple seqids

cat(sprintf("Found %d genes with multiple different seqids\n", nrow(gene_seqid_count)))

if (nrow(gene_seqid_count) == 0) {
  stop("Error: No genes with multiple seqids found\n")
}

# Display information for first few genes
cat("\nSeqid counts for first 10 genes:\n")
print(head(gene_seqid_count, 10))

# ============================================================================
# 5. Generate pairwise combinations
# ============================================================================

cat("\nGenerating pairwise combinations of seqids...\n")

# Generate all pairwise combinations for each gene
pairs_list <- list()

for (i in 1:nrow(gene_seqid_count)) {
  gene <- gene_seqid_count$Gene[i]
  
  # Get all seqids and their PHI values for this gene
  gene_data <- df_clean %>%
    filter(Gene == gene) %>%
    select(seqid, PHI)
  
  seqids <- gene_data$seqid
  n_seqid <- length(seqids)
  
  if (n_seqid < 2) {
    next
  }
  
  # Generate all pairwise combinations
  # Use combn function to generate combination indices
  comb_indices <- combn(1:n_seqid, 2)
  
  # Create data for each pair
  for (j in 1:ncol(comb_indices)) {
    idx1 <- comb_indices[1, j]
    idx2 <- comb_indices[2, j]
    
    pairs_list[[length(pairs_list) + 1]] <- data.frame(
      Gene = gene,
      seqid1 = seqids[idx1],
      seqid2 = seqids[idx2],
      PHI1 = gene_data$PHI[idx1],
      PHI2 = gene_data$PHI[idx2],
      stringsAsFactors = FALSE
    )
  }
  
  if (i %% 10 == 0) {
    cat(sprintf("  Processed %d/%d genes...\n", i, nrow(gene_seqid_count)))
  }
}

# Combine all combinations
df_pairs <- do.call(rbind, pairs_list)

cat(sprintf("Generation completed: %d combination pairs\n", nrow(df_pairs)))

# ============================================================================
# 6. Create output directory
# ============================================================================

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", output_dir))
}

# ============================================================================
# 7. Save gene and seqid combination information
# ============================================================================

cat("\nSaving gene and seqid combination information...\n")

# Save combination pair information (Gene, seqid1, seqid2)
genes_combinations_file <- file.path(output_dir, "genes_with_combinations.csv")
write.csv(
  df_pairs %>% select(Gene, seqid1, seqid2),
  genes_combinations_file,
  row.names = FALSE
)
cat(sprintf("Gene combination information saved to: %s\n", genes_combinations_file))

# ============================================================================
# 8. Save all combination data
# ============================================================================

cat("\nSaving all combination data...\n")

pairs_data_file <- file.path(output_dir, "PHI_pairs_data.csv")
write.csv(df_pairs, pairs_data_file, row.names = FALSE)
cat(sprintf("Combination data saved to: %s\n", pairs_data_file))

# ============================================================================
# 9. Create scatter plot
# ============================================================================

cat("\nCreating scatter plot...\n")

# Calculate correlation coefficient and p-value
cor_test <- cor.test(df_pairs$PHI1, df_pairs$PHI2, method = "pearson")
cor_coef <- cor_test$estimate
cor_pval <- cor_test$p.value

cat(sprintf("Correlation coefficient: %.4f\n", cor_coef))
cat(sprintf("P-value: %.2e\n", cor_pval))

# Format correlation coefficient and p-value text
cor_text <- sprintf("r = %.2f", cor_coef)
if (cor_pval < 0.001) {
  pval_text <- "P < 0.001"
} else {
  pval_text <- sprintf("P = %.3f", cor_pval)
}

# Create scatter plot with optimized color scheme
p <- ggplot(df_pairs, aes(x = PHI1, y = PHI2)) +
  geom_point(alpha = 0.5, size = 1.3, color = "#4472C4") +
  geom_smooth(method = "lm", se = TRUE, color = "#C55A11", linetype = "solid", linewidth = 1.2, fill = "#FFC000", alpha = 0.25) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray60", linewidth = 0.6) +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = paste(cor_text, pval_text, sep = "\n"),
    hjust = 1.1, vjust = -0.5,
    size = 5.5,
    color = "black"
  ) +
  labs(
    x = "PHI1",
    y = "PHI2"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

# Save plot
plot_file <- file.path(output_dir, "PHI_scatter_plot.png")
ggsave(
  plot_file,
  plot = p,
  width = 5,
  height = 5,
  dpi = 300
)
cat(sprintf("Scatter plot saved to: %s\n", plot_file))

# ============================================================================
# 10. Output statistics
# ============================================================================

cat("\n=== Statistics ===\n")
cat(sprintf("Number of genes found: %d\n", nrow(gene_seqid_count)))
cat(sprintf("Total combination pairs: %d\n", nrow(df_pairs)))
cat(sprintf("Average pairs per gene: %.2f\n", nrow(df_pairs) / nrow(gene_seqid_count)))

cat("\nSeqid counts and pair counts for each gene:\n")
gene_stats <- df_pairs %>%
  group_by(Gene) %>%
  summarise(
    n_seqid = n_distinct(c(seqid1, seqid2)),
    n_pairs = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(n_pairs))

print(gene_stats)

cat("\nPHI value statistics:\n")
cat(sprintf("  PHI1 range: %.4f ~ %.4f\n", min(df_pairs$PHI1, na.rm = TRUE), max(df_pairs$PHI1, na.rm = TRUE)))
cat(sprintf("  PHI2 range: %.4f ~ %.4f\n", min(df_pairs$PHI2, na.rm = TRUE), max(df_pairs$PHI2, na.rm = TRUE)))
cat(sprintf("  PHI1 mean: %.4f\n", mean(df_pairs$PHI1, na.rm = TRUE)))
cat(sprintf("  PHI2 mean: %.4f\n", mean(df_pairs$PHI2, na.rm = TRUE)))

cat("\nCorrelation analysis:\n")
cat(sprintf("  Correlation coefficient (r): %.4f\n", cor_coef))
cat(sprintf("  P-value: %.2e\n", cor_pval))

cat("\nCompleted!\n")
