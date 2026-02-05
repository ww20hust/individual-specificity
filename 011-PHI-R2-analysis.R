#!/usr/bin/env Rscript
# 011-PHI-R2-analysis.R
# 绘制PHI与linear_r2的散点图，筛选PHI>0.7且gap_p<0.05的蛋白，绘制综合图

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)

options(warn = -1)

# =========================
# 路径设置
# =========================

# 输入文件路径
survival_results_path <- "/home/ww/Project/Protein-Composition/result/004-longevity-target/protein_gap_survival_results.csv"
phi_file <- "./results/001-PHI-different-visit-results.csv"

# 输出目录
out_root <- "./results/Person-Setpoint"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# 注释文件路径（用于获取gene name）
anno_path <- "/home/ww/Project/AnzhenData/annoinfo_df_SOMA.csv"

cat("========================================\n")
cat("PHI与R2散点图及综合图绘制\n")
cat("========================================\n\n")

# =========================
# 函数：查找列名（处理R可能添加的X前缀）
# =========================
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

# =========================
# 读取数据
# =========================
cat("正在读取数据...\n")

# 读取生存分析结果
survival_df <- read_csv(survival_results_path, show_col_types = FALSE)
cat(sprintf("  读取 survival results: %d 行 x %d 列\n", nrow(survival_df), ncol(survival_df)))

# 检查列名
cat("\nSurvival results 列名:\n")
print(names(survival_df))

# 读取PHI数据
if (!file.exists(phi_file)) {
  stop(sprintf("错误: 找不到PHI数据文件: %s\n", phi_file))
}

df_phi <- read_csv(phi_file, show_col_types = FALSE)
cat(sprintf("  读取 PHI data: %d 行 x %d 列\n", nrow(df_phi), ncol(df_phi)))

cat("\nPHI data 列名:\n")
print(names(df_phi))

# 读取注释文件（用于获取gene name）
if (file.exists(anno_path)) {
  anno_df <- read_csv(anno_path, show_col_types = FALSE)
  cat(sprintf("  读取 annotation data: %d 行 x %d 列\n", nrow(anno_df), ncol(anno_df)))
  
  # 创建AptName到EntrezGeneSymbol的映射
  if ("AptName" %in% names(anno_df) && "EntrezGeneSymbol" %in% names(anno_df)) {
    aptname_to_gene <- setNames(anno_df$EntrezGeneSymbol, anno_df$AptName)
    cat("  已创建AptName到EntrezGeneSymbol的映射\n")
  } else {
    aptname_to_gene <- NULL
    cat("  警告: 注释文件中未找到AptName或EntrezGeneSymbol列\n")
  }
} else {
  aptname_to_gene <- NULL
  cat("  警告: 找不到注释文件，将无法显示gene name\n")
}

# 查找必要的列
aptname_col <- find_column(df_phi, "AptName")
phi_col <- find_column(df_phi, "07_12_corr")

if (is.null(aptname_col)) {
  stop("错误: 找不到AptName列\n")
}
if (is.null(phi_col)) {
  stop("错误: 找不到07_12_corr列\n")
}

cat(sprintf("  使用列: AptName=%s, PHI=%s\n", aptname_col, phi_col))

# =========================
# 数据匹配
# =========================
cat("\n正在匹配数据...\n")

# 从PHI数据中提取AptName和07_12_corr
phi_data <- df_phi %>%
  select(AptName = !!aptname_col, PHI = !!phi_col) %>%
  filter(!is.na(AptName) & !is.na(PHI))

cat(sprintf("  PHI数据有效行数: %d\n", nrow(phi_data)))

# 合并数据：seqid直接作为AptName匹配
merged_data <- survival_df %>%
  left_join(phi_data, by = c("seqid" = "AptName"))

# 添加gene_symbol信息
if (!is.null(aptname_to_gene)) {
  merged_data <- merged_data %>%
    mutate(
      gene_symbol = ifelse(seqid %in% names(aptname_to_gene), 
                          aptname_to_gene[seqid], 
                          NA_character_),
      # 处理缺失或无效的gene_symbol
      gene_symbol = ifelse(is.na(gene_symbol) | gene_symbol == "" | gene_symbol == "Null", 
                          seqid, 
                          gene_symbol)
    )
} else {
  merged_data <- merged_data %>%
    mutate(gene_symbol = seqid)
}

cat(sprintf("  合并后数据: %d 行\n", nrow(merged_data)))
cat(sprintf("  有PHI值的蛋白数: %d\n", sum(!is.na(merged_data$PHI))))
cat(sprintf("  有linear_r2值的蛋白数: %d\n", sum(!is.na(merged_data$linear_r2))))

# =========================
# 绘制PHI与linear_r2散点图
# =========================
cat("\n正在绘制PHI与linear_r2散点图...\n")

# 准备散点图数据（去除缺失值）
scatter_data <- merged_data %>%
  filter(!is.na(PHI) & !is.na(linear_r2))

cat(sprintf("  有效数据点: %d\n", nrow(scatter_data)))

if (nrow(scatter_data) > 0) {
  # 计算相关系数和p-value
  cor_test_result <- cor.test(scatter_data$PHI, scatter_data$linear_r2, method = "pearson")
  cor_coef <- cor_test_result$estimate
  cor_pvalue <- cor_test_result$p.value
  
  cat(sprintf("  相关系数: %.4f\n", cor_coef))
  cat(sprintf("  p-value: %.4e\n", cor_pvalue))
  
  # 创建散点图
  scatter_plot <- ggplot(scatter_data, aes(x = PHI, y = linear_r2)) +
    geom_point(
      color = "#E74C3C",
      alpha = 0.7,
      size = 2,
      shape = 16
    ) +
    # 添加回归线
    geom_smooth(
      method = "lm",
      se = TRUE,
      color = "#2C3E50",
      linetype = "dashed",
      linewidth = 0.8,
      alpha = 0.3
    ) +
    # 添加相关系数和p-value标注（右下角，不加粗）
    annotate(
      "text",
      x = Inf,
      y = -Inf,
      label = sprintf("r = %.3f\np = %.2e", cor_coef, cor_pvalue),
      hjust = 1.1,
      vjust = -0.5,
      size = 4,
      color = "#2C3E50"
    ) +
    # 坐标轴标签
    labs(
      x = "PHI",
      y = "Linear R² (prediction accuracy)",
      title = "PHI vs Prediction Accuracy (R²)"
    ) +
    # 主题设置
    theme_minimal(base_size = 14) +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.6),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 14, color = "#2C3E50"),
      axis.text = element_text(size = 12, color = "#2C3E50"),
      plot.title = element_text(size = 16, color = "#2C3E50", hjust = 0.5, margin = margin(b = 15)),
      plot.margin = margin(20, 20, 20, 20),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  # 保存散点图
  scatter_filename <- "PHI_vs_linear_r2_scatter.png"
  scatter_path <- file.path(out_root, scatter_filename)
  
  ggsave(
    scatter_path,
    scatter_plot,
    width = 4,
    height = 4,
    dpi = 300,
    bg = "white"
  )
  
  cat(sprintf("✅ 散点图已保存: %s\n", scatter_path))
} else {
  cat("  警告：没有有效数据用于绘制散点图\n")
}

# =========================
# 筛选蛋白：PHI > 0.7 且 gap_p < 0.05
# =========================
cat("\n正在筛选蛋白...\n")

filtered_data <- merged_data %>%
  filter(
    !is.na(PHI) & PHI > 0.7,
    !is.na(gap_p) & gap_p < 0.05
  ) %>%
  arrange(gap_p) %>%
  mutate(
    protein_order = row_number(),
    seqid_factor = factor(seqid, levels = rev(unique(seqid))),
    # 创建protein_label：gene_symbol (seqid)
    protein_label = ifelse(gene_symbol == seqid, 
                          seqid, 
                          paste0(gene_symbol, " (", seqid, ")"))
  )

cat(sprintf("  筛选条件: PHI > 0.7 且 gap_p < 0.05\n"))
cat(sprintf("  筛选出的蛋白数: %d\n", nrow(filtered_data)))

if (nrow(filtered_data) == 0) {
  cat("  警告：没有蛋白满足筛选条件\n")
  stop("没有蛋白满足筛选条件，无法继续绘制综合图\n")
}

# =========================
# 统计信息：PHI > 0.7
# =========================
cat("\n正在统计PHI > 0.7的蛋白...\n")

phi_07_data <- merged_data %>%
  filter(!is.na(PHI) & PHI > 0.7)

phi_07_count <- nrow(phi_07_data)
cat(sprintf("  PHI > 0.7的蛋白数量: %d\n", phi_07_count))

if (phi_07_count > 0) {
  # 计算PHI > 0.7的蛋白中，最小的预测准确性
  phi_07_with_r2 <- phi_07_data %>%
    filter(!is.na(linear_r2))
  
  if (nrow(phi_07_with_r2) > 0) {
    min_r2 <- min(phi_07_with_r2$linear_r2, na.rm = TRUE)
    cat(sprintf("  PHI > 0.7的蛋白中，最小的预测准确性 (linear_r2): %.4f\n", min_r2))
  } else {
    cat("  警告：PHI > 0.7的蛋白中没有linear_r2数据\n")
    min_r2 <- NA
  }
} else {
  cat("  警告：没有蛋白满足PHI > 0.7\n")
  min_r2 <- NA
}

# =========================
# 保存筛选后的数据到CSV
# =========================
cat("\n正在保存筛选后的数据...\n")

output_csv_path <- file.path(out_root, "filtered_proteins_phi_r2_survival.csv")
write_csv(filtered_data, output_csv_path)
cat(sprintf("✅ 筛选后的数据已保存: %s\n", output_csv_path))

# =========================
# 保存显著蛋白的gene name列表到txt
# =========================
cat("\n正在保存显著蛋白的gene name列表...\n")

# 提取显著的gene name列表（去除重复和NA）
significant_genes <- filtered_data %>%
  filter(!is.na(gene_symbol) & gene_symbol != "" & gene_symbol != "Null") %>%
  distinct(gene_symbol) %>%
  arrange(gene_symbol) %>%
  pull(gene_symbol)

# 如果gene_symbol就是seqid，则使用seqid
if (length(significant_genes) == 0 || all(significant_genes %in% filtered_data$seqid)) {
  significant_genes <- filtered_data %>%
    distinct(seqid) %>%
    arrange(seqid) %>%
    pull(seqid)
}

gene_list_path <- file.path(out_root, "significant_proteins_gene_list.txt")

# 写入txt文件
writeLines(significant_genes, gene_list_path)
cat(sprintf("✅ 显著蛋白的gene name列表已保存: %s\n", gene_list_path))
cat(sprintf("  共 %d 个显著的蛋白\n", length(significant_genes)))

# =========================
# 绘制综合图：HR森林图 + R2线图
# =========================
cat("\n正在绘制综合图...\n")

# 准备综合数据
comprehensive_data <- filtered_data %>%
  filter(
    !is.na(gap_hr) & !is.na(gap_hr_ci_lower) & !is.na(gap_hr_ci_upper),
    !is.na(base_model_r2) & !is.na(full_model_r2)
  )

cat(sprintf("  综合数据蛋白数量: %d\n", nrow(comprehensive_data)))

if (nrow(comprehensive_data) > 0) {
  # 准备数据用于绘图
  plot_comprehensive <- comprehensive_data %>%
    mutate(
      gap_p_label = sprintf("%.4f", gap_p),
      lrt_p_label = ifelse(is.na(lrt_p), "NA", 
                          ifelse(lrt_p < 0.001, "<0.001", sprintf("%.4f", lrt_p)))
    )
  
  # 1. HR森林图
  hr_plot <- ggplot(plot_comprehensive, aes(x = seqid_factor, y = gap_hr)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "#95A5A6", linewidth = 0.5) +
    geom_errorbar(
      aes(ymin = gap_hr_ci_lower, ymax = gap_hr_ci_upper),
      width = 0.15,
      color = "#2C3E50",
      linewidth = 0.4
    ) +
    geom_point(color = "#E74C3C", size = 1.5, shape = 16) +
    geom_text(
      aes(x = seqid_factor, y = Inf, label = gap_p_label),
      inherit.aes = FALSE,
      data = plot_comprehensive,
      hjust = -0.1,
      vjust = 0.5,
      size = 2.5,
      color = "#2C3E50"
    ) +
    scale_x_discrete(labels = setNames(plot_comprehensive$protein_label, plot_comprehensive$seqid_factor)) +
    coord_flip(clip = "off") +
    labs(x = "", y = "HR (p)") +
    theme_minimal(base_size = 10) +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.4),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.1),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 8, color = "#2C3E50"),
      axis.text.x = element_text(size = 8, color = "#2C3E50"),
      axis.title = element_text(size = 9, color = "#2C3E50"),
      plot.margin = margin(5, 50, 5, 5),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  # 2. R2线图
  # 转换为长格式用于绘制连线
  r2_data_long <- plot_comprehensive %>%
    select(seqid, seqid_factor, base_model_r2, full_model_r2, lrt_p_label) %>%
    pivot_longer(
      cols = c(base_model_r2, full_model_r2),
      names_to = "model_type",
      values_to = "r2_value"
    ) %>%
    mutate(
      model_label = case_when(
        model_type == "base_model_r2" ~ "Base Model",
        model_type == "full_model_r2" ~ "Full Model",
        TRUE ~ NA_character_
      )
    )
  
  r2_plot <- ggplot(r2_data_long, aes(x = seqid_factor, y = r2_value, group = model_type, color = model_label)) +
    geom_line(linewidth = 0.8, alpha = 0.8) +
    geom_point(size = 1.5, alpha = 0.8) +
    geom_text(
      aes(x = seqid_factor, y = Inf, label = lrt_p_label),
      inherit.aes = FALSE,
      data = plot_comprehensive,
      hjust = -0.1,
      vjust = 0.5,
      size = 2.5,
      color = "#2C3E50"
    ) +
    scale_color_manual(
      values = c(
        "Base Model" = "#3498DB",
        "Full Model" = "#E74C3C"
      ),
      name = NULL
    ) +
    scale_x_discrete(labels = setNames(plot_comprehensive$protein_label, plot_comprehensive$seqid_factor)) +
    coord_flip(clip = "off") +
    labs(x = "", y = "R² (LR p)") +
    theme_minimal(base_size = 10) +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.4),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.1),
      panel.grid.minor = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 8, color = "#2C3E50"),
      axis.title = element_text(size = 9, color = "#2C3E50"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size = 7),
      legend.margin = margin(t = 2, b = 2),
      plot.margin = margin(5, 50, 5, 5),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  # 组合图形
  comprehensive_plot <- hr_plot + r2_plot +
    plot_layout(ncol = 2, widths = c(1.2, 1))
  
  # 保存综合图
  comprehensive_plot_filename <- "comprehensive_plot_hr_r2.png"
  comprehensive_plot_path <- file.path(out_root, comprehensive_plot_filename)
  
  ggsave(
    comprehensive_plot_path,
    comprehensive_plot,
    width = 8,
    height = max(8, nrow(plot_comprehensive) * 0.15),
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
  
  cat(sprintf("✅ 综合图已保存: %s\n", comprehensive_plot_path))
} else {
  cat("  警告：没有有效的综合数据用于绘制综合图\n")
}

# =========================
# 输出最终统计信息
# =========================
cat("\n========================================\n")
cat("统计信息汇总\n")
cat("========================================\n")
cat(sprintf("PHI > 0.7的蛋白数量: %d\n", phi_07_count))
if (!is.na(min_r2)) {
  cat(sprintf("PHI > 0.7的蛋白中，最小的预测准确性 (linear_r2): %.4f\n", min_r2))
}
cat(sprintf("筛选出的显著蛋白数量 (PHI > 0.7 且 gap_p < 0.05): %d\n", nrow(filtered_data)))
cat(sprintf("显著蛋白的gene name列表已保存到: %s\n", gene_list_path))
cat(sprintf("  共 %d 个显著的蛋白\n", length(significant_genes)))

cat("\n========================================\n")
cat("完成！\n")
cat("========================================\n")
