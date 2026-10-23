######################
#--------- Figure 5b----------
#######################

library(ggplot2)
library(dplyr)
library(ggpubr)
library(data.table)

file_stats <- fread('F5b_file_stats.csv')

# 从已有数据中提取LnF和LpF的diff值
if (exists("file_stats")) {
  # 确保只包含LnF和LpF
  diff_data <- file_stats %>%
    filter(tissue %in% c("LnF", "LpF")) %>%
    select(file, tissue, diff) %>%
    mutate(tissue = factor(tissue, levels = c("LnF", "LpF")))
  
  # 计算Wilcoxon秩和检验
  LnF_diff <- diff_data %>% filter(tissue == "LnF") %>% pull(diff)
  LpF_diff <- diff_data %>% filter(tissue == "LpF") %>% pull(diff)
  
  if (length(LnF_diff) > 0 && length(LpF_diff) > 0) {
    wilcox_test <- wilcox.test(diff ~ tissue, data = diff_data)
    
    # 计算基本统计量
    stats_summary <- diff_data %>%
      group_by(tissue) %>%
      summarise(
        n = n(),
        mean = mean(diff, na.rm = TRUE),
        median = median(diff, na.rm = TRUE),
        sd = sd(diff, na.rm = TRUE),
        se = sd / sqrt(n),
        .groups = "drop"
      )
    
    # 定义颜色
    tissue_colors <- c("LnF" = "#87CEEB", "LpF" = "#CD5C5C")
    
    # 创建箱线图+点图
    p_box <- ggplot(diff_data, aes(x = tissue, y = diff, fill = tissue)) +
      # 箱线图
      geom_boxplot(
        width = 0.6,
        alpha = 0.7,
        outlier.shape = NA,  # 不显示默认的离群点
        linewidth = 0.8
      ) +
      # 点图（带有抖动）
      geom_jitter(
        aes(color = tissue),
        width = 0.2,
        height = 0,
        size = 2.5,
        alpha = 0.7,
        shape = 19
      ) +
      # 添加平均值点（黑色菱形）
      geom_point(
        data = stats_summary,
        aes(x = tissue, y = mean),
        shape = 23,  # 菱形
        size = 4,
        fill = "black",
        color = "black",
        stroke = 1.2
      ) +
      # 添加平均值标签
      geom_text(
        data = stats_summary,
        aes(x = tissue, y = mean, label = sprintf("Mean: %.3f", mean)),
        vjust = 1.8,
        size = 3.2,
        fontface = "bold",
        color = "darkblue"
      ) +
      # 颜色设置
      scale_fill_manual(values = tissue_colors, name = "Tissue Type") +
      scale_color_manual(values = tissue_colors, guide = "none") +
      # 坐标轴和标签
      labs(
        title = "Comparison of Difference Values: LnF vs LpF",
        subtitle = sprintf("Wilcoxon rank-sum test: p = %.4f", wilcox_test$p.value),
        x = "Tissue Type",
        y = "Difference (MIMER mean - Random mean)",
        caption = sprintf("Sample size: LnF = %d, LpF = %d", 
                         stats_summary$n[stats_summary$tissue == "LnF"],
                         stats_summary$n[stats_summary$tissue == "LpF"])
      ) +
      # 调整y轴范围，为标签留出空间
      ylim(
        min(diff_data$diff, na.rm = TRUE) * 1.1,
        max(diff_data$diff, na.rm = TRUE) * 1.1
      ) +
      # 主题设置
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(size = 13, face = "bold", margin = margin(r = 10)),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
        plot.caption = element_text(hjust = 0.5, size = 10, color = "gray50"),
        plot.margin = margin(15, 15, 15, 15)
      )
    
    # 添加显著性标记
    p_box <- p_box + 
      stat_compare_means(
        method = "wilcox.test",
        comparisons = list(c("LnF", "LpF")),
        label = "p.signif",
        symnum.args = list(
          cutpoints = c(0, 0.001, 0.01, 0.05, 1),
          symbols = c("***", "**", "*", "NS")
        ),
        vjust = 0.2,
        size = 5
      )
    
    # 打印图表
    print(p_box)
    
    # 打印统计摘要
    cat("\n=== Statistical Summary ===\n")
    print(stats_summary)
    cat("\n=== Wilcoxon Rank-Sum Test Results ===\n")
    print(wilcox_test)
    
  } else {
    cat("Error: At least one group has zero samples\n")
  }
} else {
  cat("Error: file_stats data frame not found. Please run the previous code first.\n")
}
