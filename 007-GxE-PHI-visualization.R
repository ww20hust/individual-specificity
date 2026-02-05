#!/usr/bin/env Rscript
# 007-GxE-PHI-visualization.R
# Read GxE data and match with PHI data
# Visualize Exposures (E) distribution grouped by PHI
# Categorize PHI into 4 groups: 0-0.25, 0.25-0.5, 0.5-0.75, 0.75-1
# Create scientific-style violin plots and perform significance tests

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggsignif)

# ============================================================================
# 1. Configuration parameters
# ============================================================================

# GxE data file path
gxe_file <- "/home/ww/Project/AnzhenData/figure6-genetic-effect/final result/GxE-SEQID.csv"

# PHI results file path
phi_file <- "./results/001-PHI-different-visit-results.csv"

# Output directory
output_base <- "./results"
output_dir <- file.path(output_base, "GxE-PHI")
figs_dir <- file.path(output_dir, "Figs")

# ============================================================================
# 2. Read GxE data
# ============================================================================

cat("=== GxE-PHI Visualization Analysis ===\n\n")
cat("Reading GxE data...\n")

# Read GxE data file
if (!file.exists(gxe_file)) {
  stop(sprintf("Error: File not found: %s\n", gxe_file))
}

df_gxe <- read.csv(gxe_file, stringsAsFactors = FALSE, check.names = FALSE)

cat(sprintf("GxE data loaded: %d rows, %d columns\n", nrow(df_gxe), ncol(df_gxe)))
cat("Column names:\n")
print(colnames(df_gxe))

# Check required columns
if (!"Protein_seq" %in% colnames(df_gxe)) {
  stop("Error: Protein_seq column not found\n")
}

# Find Exposures (E) column (may contain spaces or special characters)
exposure_col <- NULL
possible_names <- c("Exposures (E)", "Exposures.E.", "Exposures_E", "Exposures")
for (name in possible_names) {
  if (name %in% colnames(df_gxe)) {
    exposure_col <- name
    break
  }
}

if (is.null(exposure_col)) {
  # Try fuzzy matching
  matches <- grep("Exposures", colnames(df_gxe), ignore.case = TRUE, value = TRUE)
  if (length(matches) > 0) {
    exposure_col <- matches[1]
    cat(sprintf("Using column '%s' as Exposures (E)\n", exposure_col))
  } else {
    stop("Error: Exposures (E) column not found\n")
  }
}

cat(sprintf("Using columns: Protein_seq=%s, Exposures (E)=%s\n", "Protein_seq", exposure_col))

# ============================================================================
# 3. Read PHI results data
# ============================================================================

cat("\nReading PHI results data...\n")

if (!file.exists(phi_file)) {
  stop(sprintf("Error: File not found: %s\n", phi_file))
}

df_phi <- read.csv(phi_file, stringsAsFactors = FALSE, check.names = FALSE)

cat(sprintf("PHI data loaded: %d rows, %d columns\n", nrow(df_phi), ncol(df_phi)))

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
aptname_col <- find_column(df_phi, "AptName")
phi_col <- find_column(df_phi, "12_22_corr")

if (is.null(aptname_col)) {
  stop("Error: AptName column not found\n")
}
if (is.null(phi_col)) {
  stop("Error: 12_22_corr column not found\n")
}

cat(sprintf("Using columns: AptName=%s, PHI=%s\n", aptname_col, phi_col))

# ============================================================================
# 4. Data matching
# ============================================================================

cat("\nMatching data...\n")

# Prepare GxE data (keep only required columns)
df_gxe_subset <- df_gxe %>%
  select(Protein_seq, Exposure = !!exposure_col) %>%
  filter(!is.na(Protein_seq) & !is.na(Exposure))

# Prepare PHI data (keep only required columns)
df_phi_subset <- df_phi %>%
  select(AptName = !!aptname_col, PHI = !!phi_col) %>%
  filter(!is.na(AptName) & !is.na(PHI))

# Match: by seqid (Protein_seq = AptName)
df_merged <- df_phi_subset %>%
  inner_join(df_gxe_subset, by = c("AptName" = "Protein_seq"))

cat(sprintf("Matching completed: %d rows\n", nrow(df_merged)))
cat(sprintf("Matched seqid count: %d\n", length(unique(df_merged$AptName))))

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", output_dir))
}
if (!dir.exists(figs_dir)) {
  dir.create(figs_dir, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", figs_dir))
}

# ============================================================================
# 5. Z-score normalization of Exposures (E)
# ============================================================================

cat("\nPerforming Z-score normalization of Exposures (E)...\n")

# Calculate mean and standard deviation
exposure_mean <- mean(df_merged$Exposure, na.rm = TRUE)
exposure_sd <- sd(df_merged$Exposure, na.rm = TRUE)

cat(sprintf("Exposures (E) original statistics: mean=%.4f, sd=%.4f\n", exposure_mean, exposure_sd))

# Perform Z-score normalization
df_merged <- df_merged %>%
  mutate(
    Exposure_zscore = (Exposure - exposure_mean) / exposure_sd
  )

cat(sprintf("Z-score normalization completed\n"))
cat(sprintf("Normalized statistics: mean=%.4f, sd=%.4f\n", 
            mean(df_merged$Exposure_zscore, na.rm = TRUE), 
            sd(df_merged$Exposure_zscore, na.rm = TRUE)))

# Save matched data (including original and Z-score values)
merged_file <- file.path(output_dir, "PHI_GxE_merged.csv")
write.csv(df_merged, merged_file, row.names = FALSE)
cat(sprintf("Matched data saved to: %s\n", merged_file))

# ============================================================================
# 6. Create PHI groups
# ============================================================================

cat("\nCreating PHI groups...\n")

# Remove missing values and create groups (using Z-score normalized values)
df_clean <- df_merged %>%
  filter(!is.na(PHI) & !is.na(Exposure_zscore)) %>%
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

# 设置分组顺序
df_clean$PHI_bin <- factor(df_clean$PHI_bin, 
                           levels = c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1"))

cat(sprintf("清理后数据: %d 行\n", nrow(df_clean)))
cat("各分组样本数:\n")
print(table(df_clean$PHI_bin))

# ============================================================================
# 6. 进行组间显著性检验
# ============================================================================

cat("\n正在进行组间显著性检验...\n")

# 提取各组的Exposure_zscore值
group_0_25 <- df_clean$Exposure_zscore[df_clean$PHI_bin == "0-0.25"]
group_25_50 <- df_clean$Exposure_zscore[df_clean$PHI_bin == "0.25-0.5"]
group_50_75 <- df_clean$Exposure_zscore[df_clean$PHI_bin == "0.5-0.75"]
group_75_100 <- df_clean$Exposure_zscore[df_clean$PHI_bin == "0.75-1"]

# 进行两两比较（使用t检验）
comparisons <- list(
  c("0-0.25", "0.25-0.5"),
  c("0.25-0.5", "0.5-0.75"),
  c("0.5-0.75", "0.75-1"),
  c("0-0.25", "0.5-0.75"),
  c("0.25-0.5", "0.75-1"),
  c("0-0.25", "0.75-1")
)

# 存储检验结果
test_results <- data.frame(
  Group1 = character(),
  Group2 = character(),
  P_value = numeric(),
  Statistic = numeric(),
  stringsAsFactors = FALSE
)

for (comp in comparisons) {
  group1_data <- df_clean$Exposure_zscore[df_clean$PHI_bin == comp[1]]
  group2_data <- df_clean$Exposure_zscore[df_clean$PHI_bin == comp[2]]
  
  if (length(group1_data) >= 2 && length(group2_data) >= 2) {
    # 进行t检验
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
    
    # 格式化p值显示
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

# 保存检验结果
test_results_file <- file.path(output_dir, "PHI_bin_statistical_tests_GxE.csv")
write.csv(test_results, test_results_file, row.names = FALSE)
cat(sprintf("\n检验结果已保存到: %s\n", test_results_file))

# ============================================================================
# 7. 绘制小提琴图
# ============================================================================

cat("\n正在生成小提琴图...\n")

# 计算各组的中位数和均值用于标注
summary_stats <- df_clean %>%
  group_by(PHI_bin) %>%
  summarise(
    n = n(),
    median = median(Exposure_zscore, na.rm = TRUE),
    mean = mean(Exposure_zscore, na.rm = TRUE),
    .groups = 'drop'
  )

cat("各组统计摘要:\n")
print(summary_stats)

# 定义颜色（使用科研风格的配色）
color_palette <- c(
  "0-0.25" = "#E8E8E8",      # 浅灰色
  "0.25-0.5" = "#B0D4E8",    # 浅蓝色
  "0.5-0.75" = "#5BA3D0",    # 中蓝色
  "0.75-1" = "#1E5F8F"        # 深蓝色
)

# 创建小提琴图（使用Z-score标准化后的值）
p_violin <- ggplot(df_clean, aes(x = PHI_bin, y = Exposure_zscore, fill = PHI_bin)) +
  geom_violin(trim = FALSE, alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.8, 
               outlier.alpha = 0.5, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = color_palette) +
  labs(
    x = "PHI bin",
    y = "Exposures (E) Z-score"
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

# 添加显著性标注（相邻组之间的比较）
y_max <- max(df_clean$Exposure_zscore, na.rm = TRUE)
y_range <- max(df_clean$Exposure_zscore, na.rm = TRUE) - min(df_clean$Exposure_zscore, na.rm = TRUE)

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

# 保存图片（减小宽度）
output_file <- file.path(figs_dir, "violin_PHI_bins_GxE_Exposures.png")
ggsave(output_file, plot = p_violin, width = 4.5, height = 4, dpi = 300)
cat(sprintf("已保存: %s\n", output_file))

# ============================================================================
# 8. 生成带更多显著性标注的版本（所有组间比较）
# ============================================================================

cat("\n正在生成完整显著性标注版本...\n")

# 创建包含所有比较的版本（使用Z-score标准化后的值）
p_violin_full <- ggplot(df_clean, aes(x = PHI_bin, y = Exposure_zscore, fill = PHI_bin)) +
  geom_violin(trim = FALSE, alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.8, 
               outlier.alpha = 0.5, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = color_palette) +
  labs(
    x = "PHI bin",
    y = "Exposures (E) Z-score"
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

# 添加所有显著性标注
y_max_full <- max(df_clean$Exposure_zscore, na.rm = TRUE)
y_range_full <- max(df_clean$Exposure_zscore, na.rm = TRUE) - min(df_clean$Exposure_zscore, na.rm = TRUE)

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

# 保存完整版本（减小宽度）
output_file_full <- file.path(figs_dir, "violin_PHI_bins_GxE_Exposures_full_comparisons.png")
ggsave(output_file_full, plot = p_violin_full, width = 4.5, height = 4, dpi = 300)
cat(sprintf("已保存: %s\n", output_file_full))

# ============================================================================
# 9. 生成统计摘要表格
# ============================================================================

cat("\n正在生成统计摘要...\n")

# 创建详细的统计摘要（使用Z-score标准化后的值）
detailed_stats <- df_clean %>%
  group_by(PHI_bin) %>%
  summarise(
    n = n(),
    mean = mean(Exposure_zscore, na.rm = TRUE),
    median = median(Exposure_zscore, na.rm = TRUE),
    sd = sd(Exposure_zscore, na.rm = TRUE),
    se = sd / sqrt(n),
    q25 = quantile(Exposure_zscore, 0.25, na.rm = TRUE),
    q75 = quantile(Exposure_zscore, 0.75, na.rm = TRUE),
    min = min(Exposure_zscore, na.rm = TRUE),
    max = max(Exposure_zscore, na.rm = TRUE),
    .groups = 'drop'
  )

# 保存统计摘要
stats_file <- file.path(output_dir, "PHI_bin_summary_statistics_GxE.csv")
write.csv(detailed_stats, stats_file, row.names = FALSE)
cat(sprintf("统计摘要已保存到: %s\n", stats_file))

cat("\n完成！\n")
