#!/usr/bin/env Rscript
# 010-protein-sankey-retention.R
# Generate Sankey plots for specified proteins and calculate retention ratios between adjacent visits

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial)

options(warn = -1)

# ============================================================================
# 1. Configuration
# ============================================================================

# Data paths
DATA_PATH <- "./results/age_sex_normalized_batch_removed.csv"
ANNO_PATH <- "/home/ww/Project/AnzhenData/annoinfo_df_SOMA.csv"

# Output directory
OUTPUT_DIR <- "./results/Sankey-Plot"

# Color scheme
group_colors <- c("#FFE5B4", "#B8E6D1", "#D4C5E8", "#E8E8E8")
group_labels <- c("0–25%", "25–50%", "50–75%", "75–100%")

# ============================================================================
# 2. Data loading and preprocessing
# ============================================================================

load_data <- function(data_path, anno_path) {
  cat("Loading data files...\n")
  
  # Read main data file
  df <- read.csv(data_path, stringsAsFactors = FALSE)
  cat(sprintf("Data loaded: %d rows, %d columns\n", nrow(df), ncol(df)))
  
  # Filter individuals who participated in all four visits (Dataset 1, 2, 3, 4)
  df <- df %>%
    filter(Dataset %in% c(1, 2, 3, 4)) %>%
    group_by(ID) %>%
    filter(n_distinct(Dataset) == 4) %>%
    ungroup()
  
  cat(sprintf("Filtered to individuals with all four visits: %d rows\n", nrow(df)))
  
  # Read annotation file
  anno_df <- read.csv(anno_path, stringsAsFactors = FALSE)
  cat(sprintf("Annotation loaded: %d rows\n", nrow(anno_df)))
  
  # Create mapping from AptName to gene symbol
  if (!"AptName" %in% colnames(anno_df)) {
    stop("AptName column not found in annotation file")
  }
  
  # Extract AptName from seqid (part before underscore)
  extract_aptname <- function(col_name) {
    parts <- strsplit(col_name, "_")[[1]]
    return(parts[1])
  }
  
  # Create mapping from seqid to AptName
  seqid_cols <- grep("^seq\\.", colnames(df), value = TRUE)
  aptnames <- sapply(seqid_cols, extract_aptname)
  names(aptnames) <- seqid_cols
  
  # Create mapping from AptName to gene symbol
  anno_map_symbol <- setNames(anno_df$EntrezGeneSymbol, anno_df$AptName)
  
  return(list(df = df, anno_map_symbol = anno_map_symbol, aptnames = aptnames))
}

# ============================================================================
# 3. 计算保留比例
# ============================================================================

calculate_retention_ratio <- function(df, seqid, visit1, visit2) {
  # 提取数据
  df_sub <- df %>%
    filter(Dataset %in% c(visit1, visit2)) %>%
    select(ID, Dataset, all_of(seqid)) %>%
    na.omit()
  
  if (nrow(df_sub) < 10) {
    return(NULL)
  }
  
  # 在每个visit内分组
  df_grouped <- df_sub %>%
    group_by(Dataset) %>%
    mutate(
      group = cut(
        .data[[seqid]],
        breaks = quantile(.data[[seqid]], probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE),
        include.lowest = TRUE,
        labels = group_labels
      )
    ) %>%
    ungroup()
  
  # 转换为宽格式
  df_wide <- df_grouped %>%
    select(ID, Dataset, group) %>%
    pivot_wider(names_from = Dataset, values_from = group, names_prefix = "V")
  
  visit1_col <- paste0("V", visit1)
  visit2_col <- paste0("V", visit2)
  
  if (!visit1_col %in% colnames(df_wide) || !visit2_col %in% colnames(df_wide)) {
    return(NULL)
  }
  
  # 只考虑在两个visit都有数据的人
  df_wide_complete <- df_wide %>%
    filter(!is.na(.data[[visit1_col]]) & !is.na(.data[[visit2_col]]))
  
  if (nrow(df_wide_complete) < 10) {
    return(NULL)
  }
  
  # 计算每个组在visit1的人数（只考虑在两个visit都有数据的人）
  group_counts_v1 <- df_wide_complete %>%
    count(.data[[visit1_col]]) %>%
    rename(group = .data[[visit1_col]], total = n)
  
  # 计算保留在同一组的人数
  retention_counts <- df_wide_complete %>%
    group_by(.data[[visit1_col]], .data[[visit2_col]]) %>%
    summarise(count = n(), .groups = "drop") %>%
    filter(.data[[visit1_col]] == .data[[visit2_col]]) %>%
    rename(group = .data[[visit1_col]], retained = count)
  
  # 合并计算保留比例
  retention_ratios <- group_counts_v1 %>%
    left_join(retention_counts, by = "group") %>%
    mutate(
      retained = ifelse(is.na(retained), 0, retained),
      ratio = retained / total,
      ratio_pct = round(ratio * 100, 1)
    ) %>%
    select(group, total, retained, ratio, ratio_pct)
  
  return(retention_ratios)
}

# ============================================================================
# 4. 绘制Sankey图并标注保留比例
# ============================================================================

plot_protein_sankey_with_retention <- function(seqid, df, anno_map_symbol, aptnames, output_dir) {
  tryCatch({
    cat(sprintf("处理蛋白: %s\n", seqid))
    
    # 检查seqid是否存在
    if (!seqid %in% colnames(df)) {
      cat(sprintf("  警告: %s 不在数据中，跳过\n", seqid))
      return(NULL)
    }
    
    # 准备数据
    df_sub <- df %>%
      select(ID, Dataset, all_of(seqid)) %>%
      na.omit()
    
    if (nrow(df_sub) < 40) {
      cat(sprintf("  警告: %s 数据量不足，跳过\n", seqid))
      return(NULL)
    }
    
    # 在每个visit内分组
    sankey_df <- df_sub %>%
      group_by(Dataset) %>%
      mutate(
        group = cut(
          .data[[seqid]],
          breaks = quantile(.data[[seqid]], probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE),
          include.lowest = TRUE,
          labels = group_labels
        )
      ) %>%
      ungroup() %>%
      mutate(
        visit = factor(Dataset, levels = c(1, 2, 3, 4)),
        alluvium = ID
      )
    
    # 计算相邻visit之间的保留比例
    retention_12 <- calculate_retention_ratio(df, seqid, 1, 2)
    retention_23 <- calculate_retention_ratio(df, seqid, 2, 3)
    retention_34 <- calculate_retention_ratio(df, seqid, 3, 4)
    
    # 创建Sankey图
    p_sankey <- ggplot(
      sankey_df,
      aes(
        x = visit,
        stratum = group,
        alluvium = alluvium,
        fill = group
      )
    ) +
      geom_flow(alpha = 0.5) +
      geom_stratum() +
      scale_fill_manual(values = group_colors, name = "Grouped in each visit") +
      scale_x_discrete(
        labels = c("1\n(2002)", "2\n(2007)", "3\n(2012)", "4\n(2022)")
      ) +
      labs(
        x = "Visit",
        y = "Count"
      ) +
      theme_bw() +
      theme(
        axis.text = element_text(size = 14, color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 16, color = "black"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.position = "bottom"
      ) +
      guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))
    
    # 获取基因符号
    aptname <- aptnames[seqid]
    gene_symbol <- if (!is.na(aptname) && aptname %in% names(anno_map_symbol)) {
      symbol <- anno_map_symbol[[aptname]]
      if (is.na(symbol) || symbol == "" || symbol == "Null") {
        seqid
      } else {
        symbol
      }
    } else {
      seqid
    }
    
    # 添加标题
    p_sankey <- p_sankey +
      labs(title = paste0(gene_symbol, " (", seqid, ")")) +
      theme(plot.title = element_text(hjust = 0.5, size = 16))
    
    # 计算每个visit中每个组的计数和y轴位置
    # 用于确定标注的位置
    group_counts <- sankey_df %>%
      group_by(visit, group) %>%
      summarise(count = n(), .groups = "drop") %>%
      arrange(visit, group)
    
    # 计算每个组在每个visit中的累积位置（y坐标）
    group_positions <- group_counts %>%
      group_by(visit) %>%
      arrange(visit, match(group, group_labels)) %>%
      mutate(
        y_bottom = lag(cumsum(count), default = 0),
        y_center = y_bottom + count / 2,
        y_top = y_bottom + count
      ) %>%
      ungroup()
    
    # 准备标注数据框
    annotation_df <- data.frame(
      x = numeric(),
      y = numeric(),
      label = character(),
      stringsAsFactors = FALSE
    )
    
    # 为每个相邻visit对添加标注
    visit_pairs <- list(
      list(visit1 = 1, visit2 = 2, x_pos = 1.5, retention = retention_12),
      list(visit1 = 2, visit2 = 3, x_pos = 2.5, retention = retention_23),
      list(visit1 = 3, visit2 = 4, x_pos = 3.5, retention = retention_34)
    )
    
    for (pair in visit_pairs) {
      if (!is.null(pair$retention) && nrow(pair$retention) > 0) {
        # 获取第一个visit中每个组的位置
        visit1_positions <- group_positions %>%
          filter(visit == pair$visit1) %>%
          select(group, y_center)
        
        # 合并保留比例数据
        annotation_data <- pair$retention %>%
          left_join(visit1_positions, by = "group") %>%
          filter(!is.na(y_center)) %>%
          mutate(
            x = pair$x_pos,
            y = y_center,
            label = sprintf("%.1f%%", ratio_pct)
          ) %>%
          select(x, y, label)
        
        annotation_df <- rbind(annotation_df, annotation_data)
      }
    }
    
    # 在图上添加标注
    if (nrow(annotation_df) > 0) {
      p_sankey <- p_sankey +
        geom_text(
          data = annotation_df,
          aes(x = x, y = y, label = label),
          inherit.aes = FALSE,
          size = 5,
          color = "black"
        )
    }
    
    # 保存图片
    filename <- paste0(seqid, "_", gene_symbol, ".png")
    filepath <- file.path(output_dir, filename)
    
    ggsave(
      filepath,
      plot = p_sankey,
      width = 7,
      height = 6,
      dpi = 150,
      limitsize = FALSE
    )
    
    cat(sprintf("  已保存: %s\n", filename))
    return(TRUE)
    
  }, error = function(e) {
    cat(sprintf("  错误: %s - %s\n", seqid, e$message))
    return(NULL)
  })
}

# ============================================================================
# 5. 主函数
# ============================================================================

main <- function(seqid = NULL) {
  cat("\n========================================\n")
  cat("蛋白Sankey图绘制（带保留比例标注）\n")
  cat("========================================\n\n")
  
  # 创建输出目录
  if (!dir.exists(OUTPUT_DIR)) {
    dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
    cat(sprintf("创建输出目录: %s\n", OUTPUT_DIR))
  }
  
  # 加载数据
  data_result <- load_data(DATA_PATH, ANNO_PATH)
  df <- data_result$df
  anno_map_symbol <- data_result$anno_map_symbol
  aptnames <- data_result$aptnames
  
  # 获取所有seqid列
  seqid_cols <- grep("^seq\\.", colnames(df), value = TRUE)
  cat(sprintf("找到 %d 个蛋白列\n", length(seqid_cols)))
  
  # 如果指定了seqid，只处理该蛋白；否则处理所有蛋白
  if (!is.null(seqid)) {
    if (seqid %in% seqid_cols) {
      seqid_cols <- seqid
    } else {
      stop(sprintf("错误: 未找到蛋白 %s\n", seqid))
    }
  }
  
  cat(sprintf("将处理 %d 个蛋白\n\n", length(seqid_cols)))
  
  # 处理每个蛋白
  for (seqid in seqid_cols) {
    plot_protein_sankey_with_retention(seqid, df, anno_map_symbol, aptnames, OUTPUT_DIR)
  }
  
  cat("\n========================================\n")
  cat("所有图片生成完成！\n")
  cat("========================================\n")
}

# ============================================================================
# 6. 执行主函数
# ============================================================================

# 如果作为脚本运行，可以指定seqid参数
# 例如：Rscript 010-protein-sankey-retention.R seq.10000.28
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  main(seqid = args[1])
} else {
  # 交互式运行：可以在这里指定seqid，或设置为NULL处理所有蛋白
  main(seqid = "seq.35381.3")  # 处理单个蛋白
#   main(seqid = NULL)  # 处理所有蛋白
}
