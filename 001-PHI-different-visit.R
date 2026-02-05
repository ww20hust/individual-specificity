#!/usr/bin/env Rscript
# 001-PHI-different-visit.R
# Calculate Pearson correlation coefficients for each protein across different visit pairs
# Calculate correlations for 6 visit pairs: 02-07, 02-12, 02-22, 07-12, 07-22, 12-22

# Load required libraries
library(dplyr)

# ============================================================================
# 1. Data loading and preprocessing
# ============================================================================

cat("Loading data files...\n")

# Read main data file
data_path <- "/home/ww/Project/AnzhenData/Figure2/age_sex_normalized_filtered.csv"
df <- read.csv(data_path, stringsAsFactors = FALSE)

cat(sprintf("Data loaded: %d rows, %d columns\n", nrow(df), ncol(df)))

# Read annotation information file
anno_path <- "/home/ww/Project/AnzhenData/annoinfo_df_SOMA.csv"
anno_df <- read.csv(anno_path, stringsAsFactors = FALSE)

cat(sprintf("Annotation information loaded: %d rows\n", nrow(anno_df)))

# Identify protein columns (columns starting with seq.)
protein_cols <- grep("^seq\\.", colnames(df), value = TRUE)
cat(sprintf("Found %d protein columns\n", length(protein_cols)))

# Dataset mapping: Dataset 1=2002(02), 2=2007(07), 3=2012(12), 4=2022(22)
dataset_mapping <- c("1" = "02", "2" = "07", "3" = "12", "4" = "22")

# ============================================================================
# 2. Extract AptName and match annotation information
# ============================================================================

cat("\nExtracting AptName and matching annotation information...\n")

# Extract AptName from protein column names (part before underscore)
# Example: seq.10000.28_crybb2 -> seq.10000.28
extract_aptname <- function(col_name) {
  # Extract part before underscore
  parts <- strsplit(col_name, "_")[[1]]
  return(parts[1])
}

aptnames <- sapply(protein_cols, extract_aptname)
names(aptnames) <- protein_cols

# Create mapping from AptName to annotation information
# Ensure AptName column exists in annotation file
if (!"AptName" %in% colnames(anno_df)) {
  stop("AptName column not found in annotation file")
}

# Create mapping dictionaries
aptname_to_uniprot <- setNames(anno_df$UniProt, anno_df$AptName)
aptname_to_entrez <- setNames(anno_df$EntrezGeneSymbol, anno_df$AptName)

# Get annotation information for each protein
protein_annotation <- data.frame(
  Protein = protein_cols,
  AptName = aptnames,
  stringsAsFactors = FALSE
)

# Match UniProt and EntrezGeneSymbol
protein_annotation$UniProt <- sapply(protein_annotation$AptName, function(apt) {
  if (apt %in% names(aptname_to_uniprot)) {
    result <- aptname_to_uniprot[[apt]]
    if (is.na(result) || result == "" || is.null(result) || length(result) == 0) {
      return("Null")
    }
    return(as.character(result))
  } else {
    return("Null")
  }
})

protein_annotation$EntrezGeneSymbol <- sapply(protein_annotation$AptName, function(apt) {
  if (apt %in% names(aptname_to_entrez)) {
    result <- aptname_to_entrez[[apt]]
    if (is.na(result) || result == "" || is.null(result) || length(result) == 0) {
      return("Null")
    }
    return(as.character(result))
  } else {
    return("Null")
  }
})

cat(sprintf("Annotation information matching completed: %d proteins\n", nrow(protein_annotation)))

# ============================================================================
# 3. Define visit pairs
# ============================================================================

cat("\nDefining visit pairs...\n")

# 6 visit pairs
visit_pairs <- list(
  c("02", "07"),  # Dataset 1 vs 2
  c("02", "12"),  # Dataset 1 vs 3
  c("02", "22"),  # Dataset 1 vs 4
  c("07", "12"),  # Dataset 2 vs 3
  c("07", "22"),  # Dataset 2 vs 4
  c("12", "22")   # Dataset 3 vs 4
)

# Corresponding Dataset IDs
visit_to_dataset <- list(
  "02" = 1,
  "07" = 2,
  "12" = 3,
  "22" = 4
)

cat(sprintf("Defined %d visit pairs\n", length(visit_pairs)))

# ============================================================================
# 4. Calculate Pearson correlation coefficients
# ============================================================================

cat("\nStarting to calculate Pearson correlation coefficients...\n")

# Minimum sample size requirement
min_samples <- 10

# Calculate correlation coefficient for a protein across a visit pair
calculate_correlation <- function(protein_col, visit1, visit2, data_df) {
  # Get corresponding Dataset IDs
  dataset1 <- visit_to_dataset[[visit1]]
  dataset2 <- visit_to_dataset[[visit2]]
  
  # Filter data for two visits, ensure ID is character type for matching
  data1 <- data_df[data_df$Dataset == dataset1, c("ID", protein_col)]
  data2 <- data_df[data_df$Dataset == dataset2, c("ID", protein_col)]
  
  # Convert ID to character type to ensure matching
  data1$ID <- as.character(data1$ID)
  data2$ID <- as.character(data2$ID)
  
  # Find common IDs
  common_ids <- intersect(data1$ID, data2$ID)
  
  if (length(common_ids) < min_samples) {
    return(list(corr = NA, pval = NA))
  }
  
  # Create data frame for better processing
  merged_data <- merge(
    data1[data1$ID %in% common_ids, c("ID", protein_col)],
    data2[data2$ID %in% common_ids, c("ID", protein_col)],
    by = "ID",
    suffixes = c("_1", "_2")
  )
  
  # Extract values
  values1 <- merged_data[[paste0(protein_col, "_1")]]
  values2 <- merged_data[[paste0(protein_col, "_2")]]
  
  # Remove missing values
  valid_idx <- !is.na(values1) & !is.na(values2)
  values1_clean <- values1[valid_idx]
  values2_clean <- values2[valid_idx]
  
  if (length(values1_clean) < min_samples) {
    return(list(corr = NA, pval = NA))
  }
  
  # Calculate Pearson correlation coefficient
  tryCatch({
    cor_result <- cor.test(values1_clean, values2_clean, method = "pearson")
    return(list(corr = as.numeric(cor_result$estimate), pval = cor_result$p.value))
  }, error = function(e) {
    return(list(corr = NA, pval = NA))
  })
}

# Calculate correlation coefficients for all visit pairs and all proteins
results_list <- list()

for (pair_idx in seq_along(visit_pairs)) {
  visit1 <- visit_pairs[[pair_idx]][1]
  visit2 <- visit_pairs[[pair_idx]][2]
  pair_name <- paste(visit1, visit2, sep = "_")
  
  cat(sprintf("Calculating correlation coefficients for %s pair...\n", pair_name))
  
  corr_col <- paste0(pair_name, "_corr")
  pval_col <- paste0(pair_name, "_pval")
  
  # Calculate correlation coefficient for each protein
  correlations <- numeric(length(protein_cols))
  pvalues <- numeric(length(protein_cols))
  
  for (i in seq_along(protein_cols)) {
    protein <- protein_cols[i]
    result <- calculate_correlation(protein, visit1, visit2, df)
    correlations[i] <- result$corr
    pvalues[i] <- result$pval
    
    if (i %% 100 == 0) {
      cat(sprintf("  Processing progress: %d/%d\n", i, length(protein_cols)))
    }
  }
  
  results_list[[corr_col]] <- correlations
  results_list[[pval_col]] <- pvalues
  
  cat(sprintf("%s pair calculation completed\n", pair_name))
}

# ============================================================================
# 5. Build result table
# ============================================================================

cat("\nBuilding result table...\n")

# Create result data frame
result_df <- data.frame(
  AptName = protein_annotation$AptName,
  UniProt = protein_annotation$UniProt,
  EntrezGeneSymbol = protein_annotation$EntrezGeneSymbol,
  stringsAsFactors = FALSE
)

# Add correlation coefficient and p-value columns
for (col_name in names(results_list)) {
  result_df[[col_name]] <- results_list[[col_name]]
}

cat(sprintf("Result table construction completed: %d rows, %d columns\n", nrow(result_df), ncol(result_df)))

# ============================================================================
# 6. Save results
# ============================================================================

cat("\nSaving results...\n")

# Output directory
output_dir <- "./results"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", output_dir))
}

# Save result CSV file
output_file <- file.path(output_dir, "001-PHI-different-visit-results.csv")
write.csv(result_df, output_file, row.names = FALSE, quote = FALSE)

cat(sprintf("Results saved to: %s\n", output_file))
cat("\nCompleted!\n")
