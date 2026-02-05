#!/usr/bin/env Rscript
# 003-PHI-grouped-correlation.R
# Calculate Pearson correlation coefficients for each protein across specified visit pairs
# Calculate separately by sex groups and age groups

# Load required libraries
library(dplyr)

# ============================================================================
# 1. Configuration parameters (modifiable)
# ============================================================================

# Set the two visits to analyze (modifiable)
VISIT1 <- "02"  # First visit, options: "02", "07", "12", "22"
VISIT2 <- "22"  # Second visit, options: "02", "07", "12", "22"

# Age group settings (modifiable)
# Define age group boundaries (based on age at first visit)
AGE_GROUPS <- list(
  "40-65" = c(40, 65),
  "65-90" = c(65, 90)
)

# Minimum sample size requirement
MIN_SAMPLES <- 10

# Data paths
DATA_PATH <- "/home/ww/Project/AnzhenData/Figure2/age_sex_normalized_filtered.csv"
ANNO_PATH <- "/home/ww/Project/AnzhenData/annoinfo_df_SOMA.csv"

# Output directory
OUTPUT_DIR <- "./results"

# ============================================================================
# 2. Data loading and preprocessing
# ============================================================================

cat("=== 003-PHI Grouped Correlation Analysis ===\n\n")
cat(sprintf("Analyzing visit pair: %s vs %s\n", VISIT1, VISIT2))
cat("Loading data files...\n")

# Read main data file
df <- read.csv(DATA_PATH, stringsAsFactors = FALSE)

cat(sprintf("Data loaded: %d rows, %d columns\n", nrow(df), ncol(df)))

# Read annotation information file
anno_df <- read.csv(ANNO_PATH, stringsAsFactors = FALSE)

cat(sprintf("Annotation information loaded: %d rows\n", nrow(anno_df)))

# Identify protein columns (columns starting with seq.)
protein_cols <- grep("^seq\\.", colnames(df), value = TRUE)
cat(sprintf("Found %d protein columns\n", length(protein_cols)))

# Dataset mapping: Dataset 1=2002(02), 2=2007(07), 3=2012(12), 4=2022(22)
visit_to_dataset <- list(
  "02" = 1,
  "07" = 2,
  "12" = 3,
  "22" = 4
)

# Validate visit settings
if (!VISIT1 %in% names(visit_to_dataset) || !VISIT2 %in% names(visit_to_dataset)) {
  stop(sprintf("Error: visit must be one of '02', '07', '12', or '22'\n"))
}

if (VISIT1 == VISIT2) {
  stop("Error: two visits cannot be the same\n")
}

dataset1 <- visit_to_dataset[[VISIT1]]
dataset2 <- visit_to_dataset[[VISIT2]]

cat(sprintf("Dataset mapping: %s -> %d, %s -> %d\n", VISIT1, dataset1, VISIT2, dataset2))

# ============================================================================
# 3. Extract AptName and match annotation information
# ============================================================================

cat("\nExtracting AptName and matching annotation information...\n")

# Extract AptName from protein column names (part before underscore)
extract_aptname <- function(col_name) {
  parts <- strsplit(col_name, "_")[[1]]
  return(parts[1])
}

aptnames <- sapply(protein_cols, extract_aptname)
names(aptnames) <- protein_cols

# Create mapping from AptName to annotation information
if (!"AptName" %in% colnames(anno_df)) {
  stop("AptName column not found in annotation file")
}

aptname_to_uniprot <- setNames(anno_df$UniProt, anno_df$AptName)
aptname_to_entrez <- setNames(anno_df$EntrezGeneSymbol, anno_df$AptName)

# 为每个蛋白质获取注释信息
protein_annotation <- data.frame(
  Protein = protein_cols,
  AptName = aptnames,
  stringsAsFactors = FALSE
)

# 匹配UniProt和EntrezGeneSymbol
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

cat(sprintf("注释信息匹配完成: %d 个蛋白质\n", nrow(protein_annotation)))

# ============================================================================
# 4. 计算相关系数的辅助函数
# ============================================================================

# 计算一对visit中某个蛋白质的相关系数
calculate_correlation <- function(protein_col, visit1, visit2, data_df, 
                                  dataset1_id, dataset2_id, 
                                  filter_ids = NULL) {
  # 筛选两个visit的数据
  data1 <- data_df[data_df$Dataset == dataset1_id, c("ID", protein_col, "age", "sexNew")]
  data2 <- data_df[data_df$Dataset == dataset2_id, c("ID", protein_col, "age", "sexNew")]
  
  # 将ID转换为字符类型以确保匹配
  data1$ID <- as.character(data1$ID)
  data2$ID <- as.character(data2$ID)
  
  # 如果提供了ID过滤列表，只保留这些ID
  if (!is.null(filter_ids)) {
    filter_ids <- as.character(filter_ids)
    data1 <- data1[data1$ID %in% filter_ids, ]
    data2 <- data2[data2$ID %in% filter_ids, ]
  }
  
  # 找到共同的ID
  common_ids <- intersect(data1$ID, data2$ID)
  
  if (length(common_ids) < MIN_SAMPLES) {
    return(list(corr = NA, pval = NA, n = length(common_ids)))
  }
  
  # 创建数据框以便更好地处理
  merged_data <- merge(
    data1[data1$ID %in% common_ids, c("ID", protein_col)],
    data2[data2$ID %in% common_ids, c("ID", protein_col)],
    by = "ID",
    suffixes = c("_1", "_2")
  )
  
  # 提取值
  values1 <- merged_data[[paste0(protein_col, "_1")]]
  values2 <- merged_data[[paste0(protein_col, "_2")]]
  
  # 移除缺失值
  valid_idx <- !is.na(values1) & !is.na(values2)
  values1_clean <- values1[valid_idx]
  values2_clean <- values2[valid_idx]
  
  if (length(values1_clean) < MIN_SAMPLES) {
    return(list(corr = NA, pval = NA, n = length(values1_clean)))
  }
  
  # 计算Pearson相关系数
  tryCatch({
    cor_result <- cor.test(values1_clean, values2_clean, method = "pearson")
    return(list(
      corr = as.numeric(cor_result$estimate), 
      pval = cor_result$p.value,
      n = length(values1_clean)
    ))
  }, error = function(e) {
    return(list(corr = NA, pval = NA, n = length(values1_clean)))
  })
}

# ============================================================================
# 5. 按性别分组计算相关系数
# ============================================================================

cat("\n开始按性别分组计算Pearson相关系数...\n")

# 性别分组：1=男性, 2=女性
sex_groups <- list(
  "Male" = 1,
  "Female" = 2
)

sex_results_list <- list()

# 首先获取每个人在第一个visit时的性别信息
cat("  获取第一个visit的性别信息...\n")
visit1_data <- df[df$Dataset == dataset1, c("ID", "sexNew")]
visit1_data$ID <- as.character(visit1_data$ID)

for (sex_name in names(sex_groups)) {
  sex_value <- sex_groups[[sex_name]]
  cat(sprintf("  处理性别组: %s (sexNew=%d)...\n", sex_name, sex_value))
  
  # 获取该性别组在第一个visit时的ID列表
  sex_ids <- visit1_data$ID[visit1_data$sexNew == sex_value]
  sex_ids <- unique(sex_ids)
  
  cat(sprintf("    该性别组在第一个visit的人数: %d\n", length(sex_ids)))
  
  # 为每个蛋白质计算相关系数
  correlations <- numeric(length(protein_cols))
  pvalues <- numeric(length(protein_cols))
  sample_sizes <- numeric(length(protein_cols))
  
  for (i in seq_along(protein_cols)) {
    protein <- protein_cols[i]
    result <- calculate_correlation(protein, VISIT1, VISIT2, df, 
                                    dataset1, dataset2, 
                                    filter_ids = sex_ids)
    correlations[i] <- result$corr
    pvalues[i] <- result$pval
    sample_sizes[i] <- result$n
    
    if (i %% 100 == 0) {
      cat(sprintf("    处理进度: %d/%d\n", i, length(protein_cols)))
    }
  }
  
  # 保存结果
  corr_col <- paste0("Sex_", sex_name, "_corr")
  pval_col <- paste0("Sex_", sex_name, "_pval")
  n_col <- paste0("Sex_", sex_name, "_n")
  
  sex_results_list[[corr_col]] <- correlations
  sex_results_list[[pval_col]] <- pvalues
  sex_results_list[[n_col]] <- sample_sizes
  
  cat(sprintf("  %s 组计算完成\n", sex_name))
}

# ============================================================================
# 6. 按年龄分组计算相关系数（根据第一个visit时的年龄分组）
# ============================================================================

cat("\n开始按年龄分组计算Pearson相关系数...\n")
cat(sprintf("注意: 年龄分组基于每个人在第一个visit (%s) 时的年龄\n", VISIT1))

# 首先获取每个人在第一个visit时的年龄信息
cat("  获取第一个visit的年龄信息...\n")
visit1_age_data <- df[df$Dataset == dataset1, c("ID", "age")]
visit1_age_data$ID <- as.character(visit1_age_data$ID)
# 去除重复，保留每个ID在第一个visit时的年龄
visit1_age_data <- visit1_age_data[!duplicated(visit1_age_data$ID), ]

age_results_list <- list()

for (age_group_name in names(AGE_GROUPS)) {
  age_range <- AGE_GROUPS[[age_group_name]]
  age_min <- age_range[1]
  age_max <- age_range[2]
  
  cat(sprintf("  处理年龄组: %s (%d-%d岁，基于第一个visit的年龄)...\n", age_group_name, age_min, age_max))
  
  # 根据第一个visit时的年龄筛选ID
  age_group_ids <- visit1_age_data$ID[
    visit1_age_data$age >= age_min & visit1_age_data$age <= age_max
  ]
  age_group_ids <- unique(age_group_ids)
  
  cat(sprintf("    该年龄组在第一个visit的人数: %d\n", length(age_group_ids)))
  
  # 为每个蛋白质计算相关系数
  correlations <- numeric(length(protein_cols))
  pvalues <- numeric(length(protein_cols))
  sample_sizes <- numeric(length(protein_cols))
  
  for (i in seq_along(protein_cols)) {
    protein <- protein_cols[i]
    result <- calculate_correlation(protein, VISIT1, VISIT2, df, 
                                    dataset1, dataset2, 
                                    filter_ids = age_group_ids)
    correlations[i] <- result$corr
    pvalues[i] <- result$pval
    sample_sizes[i] <- result$n
    
    if (i %% 100 == 0) {
      cat(sprintf("    处理进度: %d/%d\n", i, length(protein_cols)))
    }
  }
  
  # 保存结果
  corr_col <- paste0("Age_", age_group_name, "_corr")
  pval_col <- paste0("Age_", age_group_name, "_pval")
  n_col <- paste0("Age_", age_group_name, "_n")
  
  age_results_list[[corr_col]] <- correlations
  age_results_list[[pval_col]] <- pvalues
  age_results_list[[n_col]] <- sample_sizes
  
  cat(sprintf("  %s 组计算完成\n", age_group_name))
}

# ============================================================================
# 7. 构建结果表格
# ============================================================================

cat("\n构建结果表格...\n")

# 创建结果数据框
result_df <- data.frame(
  AptName = protein_annotation$AptName,
  UniProt = protein_annotation$UniProt,
  EntrezGeneSymbol = protein_annotation$EntrezGeneSymbol,
  stringsAsFactors = FALSE
)

# 添加性别分组结果
for (col_name in names(sex_results_list)) {
  result_df[[col_name]] <- sex_results_list[[col_name]]
}

# 添加年龄分组结果
for (col_name in names(age_results_list)) {
  result_df[[col_name]] <- age_results_list[[col_name]]
}

cat(sprintf("结果表格构建完成: %d 行, %d 列\n", nrow(result_df), ncol(result_df)))

# ============================================================================
# 8. 保存结果
# ============================================================================

cat("\n保存结果...\n")

# 创建输出目录（如果不存在）
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
  cat(sprintf("创建输出目录: %s\n", OUTPUT_DIR))
}

# 创建子文件夹
output_subdir <- file.path(OUTPUT_DIR, "age-sex-seperate-result")
if (!dir.exists(output_subdir)) {
  dir.create(output_subdir, recursive = TRUE)
  cat(sprintf("创建子文件夹: %s\n", output_subdir))
}

# 保存结果CSV文件
pair_name <- paste(VISIT1, VISIT2, sep = "_")
output_file <- file.path(output_subdir, sprintf("003-PHI-grouped-correlation-%s-results.csv", pair_name))
write.csv(result_df, output_file, row.names = FALSE, quote = FALSE)

cat(sprintf("结果已保存到: %s\n", output_file))

# ============================================================================
# 9. 生成汇总统计
# ============================================================================

cat("\n生成汇总统计...\n")

# 性别分组统计
cat("\n性别分组统计:\n")
for (sex_name in names(sex_groups)) {
  corr_col <- paste0("Sex_", sex_name, "_corr")
  pval_col <- paste0("Sex_", sex_name, "_pval")
  n_col <- paste0("Sex_", sex_name, "_n")
  
  valid_idx <- !is.na(result_df[[corr_col]]) & !is.na(result_df[[pval_col]])
  if (sum(valid_idx) > 0) {
    corr_vals <- result_df[[corr_col]][valid_idx]
    pval_vals <- result_df[[pval_col]][valid_idx]
    n_vals <- result_df[[n_col]][valid_idx]
    
    cat(sprintf("  %s:\n", sex_name))
    cat(sprintf("    有效蛋白数: %d\n", sum(valid_idx)))
    cat(sprintf("    平均样本数: %.1f (范围: %d-%d)\n", 
                mean(n_vals), min(n_vals), max(n_vals)))
    cat(sprintf("    平均相关系数: %.4f\n", mean(corr_vals)))
    cat(sprintf("    中位数相关系数: %.4f\n", median(corr_vals)))
    cat(sprintf("    P<0.05的蛋白数: %d (%.2f%%)\n", 
                sum(pval_vals < 0.05), 100 * sum(pval_vals < 0.05) / length(pval_vals)))
  }
}

# 年龄分组统计
cat("\n年龄分组统计:\n")
for (age_group_name in names(AGE_GROUPS)) {
  corr_col <- paste0("Age_", age_group_name, "_corr")
  pval_col <- paste0("Age_", age_group_name, "_pval")
  n_col <- paste0("Age_", age_group_name, "_n")
  
  valid_idx <- !is.na(result_df[[corr_col]]) & !is.na(result_df[[pval_col]])
  if (sum(valid_idx) > 0) {
    corr_vals <- result_df[[corr_col]][valid_idx]
    pval_vals <- result_df[[pval_col]][valid_idx]
    n_vals <- result_df[[n_col]][valid_idx]
    
    cat(sprintf("  %s:\n", age_group_name))
    cat(sprintf("    有效蛋白数: %d\n", sum(valid_idx)))
    cat(sprintf("    平均样本数: %.1f (范围: %d-%d)\n", 
                mean(n_vals), min(n_vals), max(n_vals)))
    cat(sprintf("    平均相关系数: %.4f\n", mean(corr_vals)))
    cat(sprintf("    中位数相关系数: %.4f\n", median(corr_vals)))
    cat(sprintf("    P<0.05的蛋白数: %d (%.2f%%)\n", 
                sum(pval_vals < 0.05), 100 * sum(pval_vals < 0.05) / length(pval_vals)))
  }
}

cat("\n完成！\n")
