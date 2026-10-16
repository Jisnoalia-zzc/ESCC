# --------Figure 3a
library(dplyr)
library(tidyr)
library(data.table)
library(tibble)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(pheatmap)
library(viridis)

# 1.4. 聚类
tls_zscore <- fread('C:/Users/Ron/Desktop/Figs/Fig3/F3_a_tls_zscore.csv') %>% tibble::column_to_rownames(var='V1') %>% as.matrix()
gene_cluster <- hclust(dist(tls_zscore, method = "euclidean"), method = "ward.D2")
tls_cluster <- hclust(dist(t(tls_zscore), method = "euclidean"), method = "ward.D2")

annotation_df <- fread('C:/Users/Ron/Desktop/Figs/Fig3/F3_a_annotation_df.csv') %>% tibble::column_to_rownames(var='V1')

# 对TLS进行ward.D2聚类并切分为3类
tls_cluster_3 <- cutree(tls_cluster, k = 3)

annotation_df$TLS_cluster <- factor(tls_cluster_3[rownames(annotation_df)],
                                    levels = 3:1, #顺序按照热图从左到右为G1-G3，所以level设定反过来，3:1,保持与ZXH的标注一致
                                    labels = paste0("TLS_G", 1:3))

# 创建调色板
blue_palette <- colorRampPalette(c("grey90", "#0066CC"))(100)
red_pallette <- colorRampPalette(c("grey90", "darkred"))(100)
lightblue_palette <- colorRampPalette(c("grey90", "#33A1DE"))(100)
lightred_palette <- colorRampPalette(c("grey90", "#FF6666"))(100)
darkgreen_palette <- colorRampPalette(c("grey90", "#006d2c"))(100)
lightgreen_palette <- colorRampPalette(c("grey90", "#66C2A5"))(100)  # 新增浅绿色

tls_cluster_colors <- setNames(brewer.pal(3, "Set1"),
                               levels(annotation_df$TLS_cluster))

annotation_colors <- list(
  patient = setNames(c(brewer.pal(8, "Set1"),brewer.pal(8, "Set3"))[seq_len(nlevels(factor(annotation_df$patient)))],
                     levels(factor(annotation_df$patient))),
  tissue = setNames(brewer.pal(8, "Set2")[seq_len(nlevels(factor(annotation_df$tissue)))],
                    levels(factor(annotation_df$tissue))),
  TLS_9_avg = viridis(100),
  TLS_12_avg = viridis(100),
  TLS_50_avg = viridis(100),
  
  # TLS聚类类别颜色设置
  TLS_cluster = tls_cluster_colors,
  
  # 6种BCR相关注释
  IGL.unique_seqs_avg = blue_palette,
  IGK.unique_seqs_avg = blue_palette,
  IGH.unique_seqs_avg = blue_palette,
  IGH.avg_mu_freq_avg = red_pallette,
  IGK.avg_mu_freq_avg = red_pallette,
  IGL.avg_mu_freq_avg = red_pallette,
  
  # 周边区域BCR指标
  peripheral_IGL.unique_seqs_avg = lightblue_palette,
  peripheral_IGK.unique_seqs_avg = lightblue_palette,
  peripheral_IGH.unique_seqs_avg = lightblue_palette,
  peripheral_IGH.avg_mu_freq_avg = lightred_palette,
  peripheral_IGK.avg_mu_freq_avg = lightred_palette,
  peripheral_IGL.avg_mu_freq_avg = lightred_palette,
  
  # ===== 自定义通路颜色设置 =====
  # 核心区域使用深绿渐变（与ALLOGRAFT一致）
  TLS_ALLOGRAFT_REJECTION_avg = darkgreen_palette,
  TLS_exhaustion_avg = darkgreen_palette,
  TLS_Cytolytics_effector_avg = darkgreen_palette,
  TLS_Phagocytosis_avg = darkgreen_palette,
  TLS_Anti_inflammatory_avg = darkgreen_palette,
  TLS_Angiogenesis_avg = darkgreen_palette,
  TLS_Inflammatory_avg = darkgreen_palette,
  TLS_Tumorsuppression_avg = darkgreen_palette,
  
  # 周边区域使用浅绿渐变
  peripheral_ALLOGRAFT_REJECTION_avg = lightgreen_palette,
  peripheral_exhaustion_avg = lightgreen_palette,
  peripheral_Cytolytics_effector_avg = lightgreen_palette,
  peripheral_Phagocytosis_avg = lightgreen_palette,
  peripheral_Anti_inflammatory_avg = lightgreen_palette,
  peripheral_Angiogenesis_avg = lightgreen_palette,
  peripheral_Inflammatory_avg = lightgreen_palette,
  peripheral_Tumorsuppression_avg = lightgreen_palette
)


#------------------4. 绘制热图
p <- pheatmap(
  mat = tls_zscore,
  cluster_rows = gene_cluster,
  cluster_cols = tls_cluster,
  show_rownames = TRUE,
  color = colorRampPalette(c("navy", "white", "firebrick"))(50),
  main = "TLS Subtypes by HVG Expression",
  annotation_col = annotation_df,
  annotation_colors = annotation_colors,
  fontsize_row = 8,
  fontsize_col = 8,
  annotation_legend = TRUE,
  gaps_col = c(7,49), # 三组区分的位置
  cutree_col = 3,
  silent = TRUE,
)

print(p)

#### Figure 3B
library(ggplot2)
library(dplyr)
library(concaveman)
library(ggpubr)

# 
plot_data <- read.csv("C:/Users/Ron/Desktop/Figs/Fig3/F3_b_all_spatial_data.csv")

# 
slices <- unique(plot_data$slice_name)
features <- unique(plot_data$feature)

# 
plot_spatial_simple <- function(data, target_slice, target_feature, 
                               color_palette = colorRampPalette(c("#094baf", '#ffff31', "#d40a0a"))(100)) {
 
  df <- data %>%
    filter(slice_name == target_slice, feature == target_feature)
  
  if (nrow(df) == 0) {
    stop("未找到数据: ", target_slice, " - ", target_feature)
  }
  
  df <- df %>%
    mutate(
      display_value = ifelse(TLS_cluster > 0 | peripheral_label > 0, expression_value, NA),
      point_type = case_when(
        TLS_cluster > 0 & !is_edge ~ "Core",
        TLS_cluster > 0 & is_edge ~ "Edge",
        peripheral_label > 0 ~ "Periphery",
        TRUE ~ "Other"
      ),
      imagecol = imagecol,
      imagerow = -imagerow
    )
  
  p <- ggplot(df, aes(x = imagecol, y = imagerow)) +
    geom_point(
      data = filter(df, point_type == "Other"),
      color = "gray90", size = 1, alpha = 0.6
    ) +
    geom_point(
      data = filter(df, point_type %in% c("Periphery", "Edge", "Core")),
      aes(color = display_value),
      size = 1, 
      alpha = 0.6
    ) +
    scale_color_gradientn(
      colors = color_palette,
      na.value = "transparent",
      name = target_feature
    ) +
    theme_void() +
    coord_fixed() +
    labs(title = paste0(target_slice, " : ", target_feature)) +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  if (any(df$TLS_cluster > 0, na.rm = TRUE)) {
    main_clusters <- df %>%
      filter(TLS_cluster > 0) %>%
      group_by(TLS_cluster) %>%
      summarise(size = n(), .groups = "drop") %>%
      filter(size >= 5) %>%
      arrange(desc(size))
    
    if (nrow(main_clusters) > 0) {
      for (cl in main_clusters$TLS_cluster) {
        cl_points <- df %>% filter(TLS_cluster == cl)
        hull <- concaveman::concaveman(as.matrix(cl_points[, c("imagecol", "imagerow")]))
        p <- p + 
          geom_polygon(
            data = as.data.frame(hull), 
            aes(x = V1, y = V2), 
            fill = NA, color = "black", 
            linetype = "solid", linewidth = 0.5, alpha = 0.7
          )
      }
    }
  }
  
  return(p)
}

all_plots_simple <- list()
counter <- 0

for (slice in slices) {
  for (feature in features) {
    counter <- counter + 1
    message("绘制: ", slice, " - ", feature, " (", counter, "/", length(slices) * length(features), ")")
    
    p <- plot_spatial_simple(plot_data, slice, feature)
    all_plots_simple[[counter]] <- p
  }
}

# 组合图形
combined_plot_simple <- ggarrange(
  plotlist = all_plots_simple,
  nrow = length(features),
  ncol = length(slices)
)

options(repr.plot.width = 16, repr.plot.height = 48)
print(combined_plot_simple)


# Figure3c and Figure 3d
combined_data <- fread('C:/Users/Ron/Desktop/Figs/Fig3/F3_cd_combined_data_celltype_boxplot.csv')

# 
library(tidyverse)
library(ggpubr)
library(ggsignif)

# 
cluster_colors <- c('#E41A1C', '#377EB8', '#4DAF4A')
options(repr.plot.width = 16, repr.plot.height = 12)

# 
core_data <- combined_data %>% filter(Region == "Core")
peri_data <- combined_data %>% filter(Region == "Peri")

# 
comparisons <- list(
  c("TLS_G1", "TLS_G2"),
  c("TLS_G1", "TLS_G3"),
  c("TLS_G2", "TLS_G3")
)

# 
create_abundance_plot <- function(data, region_name) {
  
  y_positions <- data %>%
    group_by(CellType) %>%
    summarise(
      max_val = max(Abundance, na.rm = TRUE),
      q3 = quantile(Abundance, 0.75, na.rm = TRUE) 
    ) %>%
    mutate(
      base_y = max_val * 1.05,
      step1 = base_y * 1.10, 
      step2 = base_y * 1.20, 
      step3 = base_y * 1.30 )
  
  #
  p <- ggplot(data, aes(x = TLS_cluster, y = Abundance, fill = TLS_cluster)) +
    geom_boxplot(width = 0.7, alpha = 0.8, outlier.size = 0, outlier.alpha = 0) +
    geom_jitter(width = 0.3, alpha = 0.5, size = 1.5, shape = 21, color = "black") +
    scale_fill_manual(values = cluster_colors) +
    facet_wrap(~ CellType, scales = "free", ncol = 5) +
    labs(
      title = paste(region_name, "Region: CellType Abundance"),
      subtitle = "Significance levels: *** p ≤ 0.001; ** p ≤ 0.01; * p ≤ 0.05",
      x = "TLS Groups",
      y = "Abundance"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray30"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      strip.background = element_rect(fill = "grey90", color = "grey50"),
      strip.text = element_text(face = "bold", size = 10),
      axis.title = element_text(face = "bold")
    )
  
  valid_celltypes <- data %>% group_by(CellType) %>% summarise(groups = n_distinct(TLS_cluster)) %>% filter(groups == 3) %>% pull(CellType)
  
  if(length(valid_celltypes) > 0) {
     for(ct in valid_celltypes) {
      y_pos <- y_positions %>% filter(CellType == ct)
      
      p <- p +
        # G1 vs G2
        geom_signif(
          data = filter(data, CellType == ct),
          comparisons = list(comparisons[[1]]),
          test = "wilcox.test",
          map_signif_level = TRUE, 
          textsize = 3,
          tip_length = 0.01,
          y_position = y_pos$step1,
          vjust = 0.5
        ) +
        # G1 vs G3
        geom_signif(
          data = filter(data, CellType == ct),
          comparisons = list(comparisons[[2]]),
          test = "wilcox.test",
          map_signif_level = TRUE,
          textsize = 3,
          tip_length = 0.01,
          y_position = y_pos$step2,
          vjust = 0.5
        ) +
        # G2 vs G3
        geom_signif(
          data = filter(data, CellType == ct),
          comparisons = list(comparisons[[3]]),
          test = "wilcox.test",
          map_signif_level = TRUE,
          textsize = 3,
          tip_length = 0.01,
          y_position = y_pos$step3,
          vjust = 0.5
        )
    }
  }

  return(p)
}

# 
core_plot <- create_abundance_plot(core_data, "Core")
peri_plot <- create_abundance_plot(peri_data, "Peri")

# 
print(core_plot)
print(peri_plot)

# #
# ggsave("TLS_core_celltype_abundance.pdf", core_plot,
#        width = 12, height = 8)
# ggsave("TLS_peri_celltype_abundance.pdf", peri_plot,
#        width = 12, height = 12)

# Figure 3e
library(ggplot2)
library(ggsignif)
library(rstatix)
library(dplyr)
library(tidyr)
library(data.table)
library(tibble)
library(patchwork)
library(RColorBrewer)
library(stringr)
library(pheatmap)
library(viridis)

# 1.4. 聚类
tls_zscore <- fread('F3_ef_tls_zscore.csv') %>% tibble::column_to_rownames(var='V1') %>% as.matrix()
gene_cluster <- hclust(dist(tls_zscore, method = "euclidean"), method = "ward.D2")
tls_cluster <- hclust(dist(t(tls_zscore), method = "euclidean"), method = "ward.D2")


annotation_df <- fread('F3_ef_annotation_df.csv') %>% tibble::column_to_rownames(var='V1')

# 对TLS进行ward.D2聚类并切分为3类
tls_cluster_3 <- cutree(tls_cluster, k = 3)

annotation_df$TLS_cluster <- factor(tls_cluster_3[rownames(annotation_df)],
                                   levels = 3:1, #顺序按照热图从左到右为G1-G3，所以level设定反过来，3:1,保持与ZXH的标注一致
                                   labels = paste0("TLS_G", 1:3))

# 1. 准备数据 ----------------------------------------------------------------
# 提取6种BCR指标数据和TLS聚类类别
bcr_data <- annotation_df[, c("TLS_cluster",
                             "IGL.unique_seqs_avg",
                             "IGK.unique_seqs_avg",
                             "IGH.unique_seqs_avg",

                             "IGL.avg_mu_freq_avg",
                             "IGK.avg_mu_freq_avg",
                             "IGH.avg_mu_freq_avg"),
                             ]

# 2. 转换数据为长格式（便于ggplot绘图）---------------------------------------
library(tidyr)
bcr_long <- pivot_longer(
  bcr_data,
  cols = -TLS_cluster,
  names_to = "BCR_metric",
  values_to = "Value"
)

# 3. 创建更友好的指标名称标签----------------------------------------
metric_labels <- c(
  "IGL.unique_seqs_avg" = "IGL Unique Seqs",
  "IGK.unique_seqs_avg" = "IGK Unique Seqs",
  "IGH.unique_seqs_avg" = "IGH Unique Seqs",

  "IGL.avg_mu_freq_avg" = "IGL Mutation Freq",
  "IGK.avg_mu_freq_avg" = "IGK Mutation Freq",
  "IGH.avg_mu_freq_avg" = "IGH Mutation Freq"
)

bcr_long$Metric_Label <- factor(bcr_long$BCR_metric, 
                               levels = names(metric_labels),
                               labels = metric_labels)

library(ggplot2)
library(ggsignif)
library(rstatix)

# 1. 准备数据（保持不变）
bcr_data <- annotation_df[, c("TLS_cluster", 
                             "IGL.unique_seqs_avg", 
                             "IGK.unique_seqs_avg", 
                             "IGH.unique_seqs_avg",
                             "IGL.avg_mu_freq_avg",
                             "IGK.avg_mu_freq_avg",
                             "IGH.avg_mu_freq_avg")]
bcr_long <- pivot_longer(
  bcr_data,
  cols = -TLS_cluster,
  names_to = "BCR_metric",
  values_to = "Value"
)


# 3. 创建更友好的指标名称标签----------------------------------------
metric_labels <- c(
  "IGL.unique_seqs_avg" = "IGL Unique Seqs",
  "IGK.unique_seqs_avg" = "IGK Unique Seqs",
  "IGH.unique_seqs_avg" = "IGH Unique Seqs",

  "IGL.avg_mu_freq_avg" = "IGL Mutation Freq",
  "IGK.avg_mu_freq_avg" = "IGK Mutation Freq",
  "IGH.avg_mu_freq_avg" = "IGH Mutation Freq"
)

bcr_long$Metric_Label <- factor(bcr_long$BCR_metric, 
                               levels = names(metric_labels),
                               labels = metric_labels)

# 2. 计算各组最大Y值（用于自动调整标记高度）
max_vals <- bcr_long %>% 
  group_by(BCR_metric) %>% 
  summarise(max_val = max(Value, na.rm = TRUE))

# 3. 设置比较组和标记高度
comparisons <- list(
  c("TLS_G1", "TLS_G2"),
  c("TLS_G1", "TLS_G3"),
  c("TLS_G2", "TLS_G3")
)

options(repr.plot.width = 12,repr.plot.height = 8)
# 4. 绘制箱线图并添加显著性标记
cluster_colors <- c('#E41A1C','#377EB8','#4DAF4A')
bcr_boxplot <- ggplot(bcr_long, aes(x = TLS_cluster, y = Value, fill = TLS_cluster)) +
  geom_boxplot(width = 0.7, alpha = 0.8, outlier.size = 0, outlier.alpha = 0) +
  geom_jitter(width = 0.3, alpha = 0.5) +
  scale_fill_manual(values = cluster_colors) +
  facet_wrap(~ Metric_Label, scales = "free_y", ncol = 3) +
  labs(title = "BCR Features", x = "TLS Groups", y = "Value") +
  theme_bw(base_size = 12) +
  
  # 关键修复：使用 ggsignif 添加标记 
  geom_signif(
    comparisons = comparisons,
    test = "wilcox.test",             # 非参数检验
    # map_signif_level = TRUE,          # 自动转星号
    tip_length = 0.01,                # 竖线长度
    textsize = 4,                     # 星号大小
    step_increase = 0.1,              # 避免重叠
  )

# 5. 保存结果
print(bcr_boxplot)

# Figure 3f
annotation_df <- fread('F3_ef_annotation_df.csv') %>% tibble::column_to_rownames(var='V1')

# 对TLS进行ward.D2聚类并切分为3类
tls_cluster_3 <- cutree(tls_cluster, k = 3)

annotation_df$TLS_cluster <- factor(tls_cluster_3[rownames(annotation_df)],
                                   levels = 3:1, #顺序按照热图从左到右为G1-G3，所以level设定反过来，3:1,保持与ZXH的标注一致
                                   labels = paste0("TLS_G", 1:3))

# 1. 准备数据 ----------------------------------------------------------------
# 已有包含TLS聚类结果的annotation_df
# 提取6种BCR指标数据和TLS聚类类别
bcr_data <- annotation_df[, c("TLS_cluster", 
                             "peripheral_IGL.unique_seqs_avg", 
                             "peripheral_IGK.unique_seqs_avg", 
                             "peripheral_IGH.unique_seqs_avg",

                             "peripheral_IGL.avg_mu_freq_avg",
                             "peripheral_IGK.avg_mu_freq_avg",
                             "peripheral_IGH.avg_mu_freq_avg"),
                             ]

# 2. 转换数据为长格式（便于ggplot绘图）---------------------------------------
library(tidyr)
bcr_long <- pivot_longer(
  bcr_data,
  cols = -TLS_cluster,
  names_to = "BCR_metric",
  values_to = "Value"
)

# 3. 创建更友好的指标名称标签----------------------------------------
metric_labels <- c(
  "peripheral_IGL.unique_seqs_avg" = "peripheral IGL Unique Seqs",
  "peripheral_IGK.unique_seqs_avg" = "peripheral IGK Unique Seqs",
  "peripheral_IGH.unique_seqs_avg" = "peripheral IGH Unique Seqs",

  "peripheral_IGL.avg_mu_freq_avg" = "peripheral IGL Mutation Freq",
  "peripheral_IGK.avg_mu_freq_avg" = "peripheral IGK Mutation Freq",
  "peripheral_IGH.avg_mu_freq_avg" = "peripheral IGH Mutation Freq"
)

bcr_long$Metric_Label <- factor(bcr_long$BCR_metric, 
                               levels = names(metric_labels),
                               labels = metric_labels)


# 2. 计算各组最大Y值（用于自动调整标记高度）
max_vals <- bcr_long %>% 
  group_by(BCR_metric) %>% 
  summarise(max_val = max(Value, na.rm = TRUE))

# 3. 设置比较组和标记高度
comparisons <- list(
  c("TLS_G1", "TLS_G2"),
  c("TLS_G1", "TLS_G3"),
  c("TLS_G2", "TLS_G3")
)

options(repr.plot.width = 12,repr.plot.height = 8)
# 4. 绘制箱线图并添加显著性标记
cluster_colors <- c('#E41A1C','#377EB8','#4DAF4A')
bcr_boxplot <- ggplot(bcr_long, aes(x = TLS_cluster, y = Value, fill = TLS_cluster)) +
  geom_boxplot(width = 0.7, alpha = 0.8, outlier.size = 0, outlier.alpha = 0) +
  geom_jitter(width = 0.3, alpha = 0.5) +
  scale_fill_manual(values = cluster_colors) +
  facet_wrap(~ Metric_Label, scales = "free_y", ncol = 3) +
  labs(title = "BCR Features", x = "TLS Groups", y = "Value") +
  theme_bw(base_size = 12) +
  
  # 关键修复：使用 ggsignif 添加标记
  geom_signif(
    comparisons = comparisons,
    test = "wilcox.test",             # 非参数检验
    # map_signif_level = TRUE,          # 自动转星号
    tip_length = 0.01,                # 竖线长度
    textsize = 4,                     # 星号大小
    step_increase = 0.1,              # 避免重叠
  )

# 5. 保存结果
print(bcr_boxplot)

# Figure 3h
# 
library(ggplot2)
library(ggpubr)
library(patchwork)
library(data.table)

setwd('C:/Users/Ron/Desktop/Figs/Fig3/')
# 

# 
plot_data <- fread("C:/Users/Ron/Desktop/Figs/Fig3/F3_h_GSVA_plot_data.csv")
gene_counts <- fread("C:/Users/Ron/Desktop/Figs/Fig3/F3_h_cluster_gene_counts.csv")

# 
plot_list <- list()
clusters <- unique(plot_data$Cluster)

for (cluster in clusters) {
  cluster_data <- plot_data[Cluster == cluster, ]
  gene_count <- gene_counts[Cluster == cluster, GeneCount]
  
  p <- ggplot(cluster_data, aes(x = Group, y = Score, fill = Group)) +
    geom_boxplot(width = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
    scale_fill_manual(values = c("#00BA38", "#F8766D")) +
    labs(title = paste0(cluster, "\n(", gene_count, " genes)"),
         y = "GSVA Score",
         x = NULL) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size=10),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size=10)) +
    stat_compare_means(
      method = "wilcox.test",
      label = "p.format",
      tip.length = 0.02,
      size = 5
    )
  
  plot_list[[cluster]] <- p
}

# 
combined_plot <- wrap_plots(plot_list, ncol = 4) + 
  plot_annotation(title = "GSVA score(Responder vs Non-Responder)",
                  subtitle = paste("DESeq2 filtered genes (p<0.2)"),
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
                                plot.subtitle = element_text(hjust = 0.5, size = 12)))

# 
legend_plot <- ggplot(plot_data, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#00BA38", "#F8766D"),
                    name = "Response Group",
                    labels = c("Responder (CR+PR)", "Non-Responder (SD+PD)")) +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"))

legend <- get_legend(legend_plot)

# 
final_plot <- combined_plot / legend + 
  plot_layout(heights = c(10, 1))

# 
options(repr.plot.width = 14, repr.plot.height = 6)
print(final_plot)

setwd('D:/Postdoc/Projects/Single_cell/ZWM')

# figure 3i
library(data.table)
library(ggplot2)
library(pROC)

# 
roc_curve_data <- fread("F3_i_ROC_curve_data.csv")
auc_results <- fread("F3_i_ROC_AUC_results.csv")
cluster_gene_counts <-fread('F3_i_cluster_gene_counts_ROC.csv')

# 
plot_list <- list()
unique_clusters <- unique(roc_curve_data$Cluster)

for (cluster in unique_clusters) {
  cluster_data <- roc_curve_data[Cluster == cluster, ]
  cluster_auc <- unique(cluster_data$AUC)
  gene_count <- cluster_gene_counts[Cluster == cluster, GeneCount] 
  
  p <- ggplot(cluster_data, aes(x = FPR, y = TPR)) +
    geom_line(color = "steelblue", size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    labs(title = paste0(cluster, "\n(n = ", gene_count, ")"),
         x = "1 - Specificity (False Positive Rate)",
         y = "Sensitivity (True Positive Rate)") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      panel.grid.minor = element_blank(),
      aspect.ratio = 1
    ) +
    annotate("text", 
             x = 0.75, y = 0.25, 
             label = paste0("AUC = ", cluster_auc),
             size = 5, 
             color = "red", 
             fontface = "bold")
  
  plot_list[[cluster]] <- p
}

# 6.
combined_roc_plot <- wrap_plots(plot_list, ncol = 4) + 
  plot_annotation(
    title = "ROC Curves of TLS Cluster Signatures",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

# 
options(repr.plot.width = 16, repr.plot.height = 4)
print(combined_roc_plot)

#### Figure 3j
# ---------------------------------------------------------------
library(survival)
library(survminer)
library(patchwork)
library(data.table)
library(ggplot2)

#-----------------------------------------------------
cat("读取保存的数据...\n")
survival_data <- fread("F3_j_survival_plot_data.csv")
cluster_summary <- fread("F3_j_cluster_summary.csv")
survfit_data <- fread("F3_j_survfit_data.csv")

# 
surv_type <- unique(cluster_summary$Survival_Type)
group_method <- unique(cluster_summary$Group_Method)
clusters <- unique(survival_data$Cluster)

cat(paste0("生存类型: ", surv_type, "\n"))
cat(paste0("分组方法: ", group_method, "\n"))
cat(paste0("Clusters: ", paste(clusters, collapse = ", "), "\n"))

# ------------------------------------------------
recreate_plot_list <- list()
recreate_table_list <- list()

for (cluster in clusters) {
  cat(paste0("处理cluster: ", cluster, "...\n"))
  
  # 
  cluster_data <- survival_data[Cluster == cluster, ]
  cluster_info <- cluster_summary[Cluster == cluster, ]
  
  # 
  plot_data <- data.frame(
    ID = cluster_data$ID,
    GSVA_Score = cluster_data$GSVA_Score,
    Status = cluster_data$Survival_Status,
    Time = cluster_data$Survival_Time,
    GSVA_Group = cluster_data$GSVA_Group
  )
  
  # 
  xlab_text <- ifelse(surv_type == "OS", "OS Days", "PFS Days")
  
  # 
  group_title <- paste0(cluster, "\n(", group_method, " split, n=", cluster_info$Gene_Count, ")")
  
  # 
  survfit_data_cluster <- survfit_data[Cluster == cluster, ]
  
  # 
  fit <- survfit(
    Surv(Time, Status) ~ GSVA_Group,
    data = plot_data
  )
  
  # 
  surv_plot <- ggsurvplot(
    fit,
    data = plot_data,
    pval = TRUE,
    conf.int = FALSE,
    risk.table = TRUE,
    palette = c("#E74C3C", "#3498DB"),
    legend.labs = c("High GSVA", "Low GSVA"),
    legend.title = "GSVA Group",
    xlab = xlab_text,
    ylab = paste0(surv_type, " Probability"),
    ggtheme = theme_bw(),
    title = group_title,
    font.title = c(11, "bold"),
    censor = TRUE,
    censor.size = 3,
    censor.shape = "|",
    censor.colour = "black"
  )
  
  recreate_plot_list[[cluster]] <- surv_plot$plot
  recreate_table_list[[cluster]] <- surv_plot$table
}

# -------------------------------------------------------------
cat("拼接生存曲线...\n")
surv_combined <- wrap_plots(recreate_plot_list, ncol = 4) + 
  plot_annotation(
    title = paste0("TLS GeneSets ", surv_type, " Survival Analysis (", group_method, " split)"),
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

# 
legend <- get_legend(
  recreate_plot_list[[1]] + 
    theme(legend.position = "bottom") +
    guides(color = guide_legend(nrow = 1))
)

final_surv_plot <- surv_combined / legend + 
  plot_layout(heights = c(10, 1))

# 
cat("显示生存曲线...\n")
options(repr.plot.width = 16, repr.plot.height = 4)
print(final_surv_plot)

# 
cat("拼接风险表...\n")
table_combined <- wrap_plots(recreate_table_list, ncol = 4) + 
  plot_annotation(
    title = paste0("TLS GeneSets ", surv_type, " Risk Table (", group_method, " split)"),
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

print(table_combined)

#--------------------------------
#------Figure 3g
#------------------------------
validate_markers_core <- fread(file=  "./Genelist/Genes_core_4groups_loop10_pct0.25_fc0.5_0.6perc_top500_GSVA_boxplot.csv")

clusters <- c('Non-TLS', 'TLS_G3', 'TLS_G2',  'TLS_G1')

DefaultAssay(tls_balanced) <- "spatial"

# 
for (cluster_type in clusters) {
  # 
  cluster_genelist <- validate_markers_core$gene[validate_markers_core$cluster == cluster_type]
  
  # 
  if (length(cluster_genelist) > 0) {
    # 
    temp_name <- paste0("temp_", cluster_type)
    
    # 
    tls_balanced <- AddModuleScore(
      object = tls_balanced,
      features = list(cluster_genelist),
      name = temp_name
    )
    
    # 
    new_col_name <- paste0(cluster_type, "_score")
    
    # 
    n <- ncol(tls_balanced@meta.data)
    temp_col <- tls_balanced@meta.data[, n, drop = FALSE]
    
    # 
    colnames(temp_col) <- new_col_name
    tls_balanced@meta.data <- cbind(tls_balanced@meta.data[, 1:(n-1)], temp_col)
  } else {
    warning(paste("No genes found for cluster:", cluster_type))
  }
}

# 
library(ggplot2)
library(patchwork)
library(dplyr)

# 
n <- grep('TLS_Cluster_result', colnames(tls_balanced@meta.data))
score_cols <- colnames(tls_balanced@meta.data)[(n+2):(n+5)]

# 
clean_cols <- gsub("[^[:alnum:]_]", "_", score_cols)
colnames(tls_balanced@meta.data)[(n+2):(n+5)] <- clean_cols

# 
plot_list <- list()

# 
for (col in clean_cols) {
  p <- ggplot(tls_balanced@meta.data, 
              aes(x = .data[["TLS_Cluster_result"]], 
                  y = .data[[col]], 
                  fill = .data[["TLS_Cluster_result"]])) +
    geom_violin()+
    geom_boxplot(width = 0.3, outlier.size = 0.8,fill='white') +
    labs(title = paste((col)), 
         x = "Region", 
         y = "Geneset Score") +
    theme_classic(base_size = 11) +
    theme(
      axis.text.y = element_text(size=10, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1,size=10, color = "black"),
      plot.title = element_text(size = 12, face = "bold"),
      legend.position = "none" 
    )
  
  plot_list[[col]] <- p
}
# 
combined_plot <- wrap_plots(plot_list, ncol = 4, nrow = 1)

options(repr.plot.width = 14, repr.plot.height = 4)
# 
print(combined_plot)

#### Figure 3k
library(data.table)
library(ggplot2)
library(pROC)

setwd('C:/Users/Ron/Desktop/Figs/Fig3')
# 8.1 
boxplot_data_reload <- fread( "F3_k_boxplot_data.csv")
roc_data_reload <- fread("F3_k_roc_data.csv")
auc_data_reload <- fread("F3_k_auc_data.csv")
surv_data_reload <- fread("F3_k_survival_data.csv")

# 8.2 
if (nrow(boxplot_data_reload) > 0) {
  cat("\n=== 从保存的数据重新绘制箱线图 ===\n")
  
  # 
  clusters_reload <- unique(boxplot_data_reload$Cluster)
  boxplot_list_reload <- list()
  
  for (cluster in clusters_reload) {
    # 
    cluster_data <- boxplot_data_reload[boxplot_data_reload$Cluster == cluster, ]
    gene_count <- unique(cluster_data$Genes)[1]
    
    # 
    p <- ggplot(cluster_data, aes(x = Group, y = Score, fill = Group)) +
      geom_boxplot(width = 0.6, outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
      scale_fill_manual(values = c("#00BA38", "#F8766D")) +
      labs(title = paste0(cluster, "\n(", gene_count, " genes)"),
           y = "GSVA Score", x = NULL) +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
            legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
      stat_compare_means(method = "wilcox.test", label = "p.format", 
                         tip.length = 0.02, size = 5)
    
    boxplot_list_reload[[cluster]] <- p
  }
  
  # 
  box_combined_reload <- wrap_plots(boxplot_list_reload, ncol = 4) +
    plot_annotation(
      title = "Melanoma-PRJEB23709 GSVA score(Responder(49) vs Non-Responder(42))",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
                    plot.subtitle = element_text(hjust = 0.5, size = 12))
    )
  
  print(box_combined_reload)
}

# 8.3 
if (nrow(roc_data_reload) > 0 && nrow(auc_data_reload) > 0) {
  cat("\n=== 从保存的数据重新绘制ROC曲线 ===\n")
  
  rocplot_list_reload <- list()
  clusters_roc <- unique(roc_data_reload$Cluster)
  
  for (cluster in clusters_roc) {
    # 
    cluster_roc_data <- roc_data_reload[roc_data_reload$Cluster == cluster, ]
    cluster_auc <- auc_data_reload[auc_data_reload$Cluster == cluster, ]$AUC[1]
    gene_count <- auc_data_reload[auc_data_reload$Cluster == cluster, ]$Genes[1]
    
    # 
    roc_p <- ggplot(cluster_roc_data, aes(x = FPR, y = TPR)) +
      geom_line(color = "steelblue", size = 1) +
      geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),
                   linetype = "dashed", color = "gray50") +
      labs(title = paste0(cluster, " (n = ", gene_count, ")"),
           x = "False Positive Rate (1 - Specificity)",
           y = "True Positive Rate (Sensitivity)") +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
            aspect.ratio = 1, panel.grid.minor = element_blank()) +
      annotate("text", x = 0.7, y = 0.15,
               label = paste("AUC =", round(cluster_auc, 3)),
               size = 4.5, color = "red", fontface = "bold")
    
    rocplot_list_reload[[cluster]] <- roc_p
  }
  
  # 
  roc_combined_reload <- wrap_plots(rocplot_list_reload, ncol = 4) +
    plot_annotation(
      title = "ROC Curves for Response Prediction",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
                    plot.subtitle = element_text(hjust = 0.5, size = 12))
    )
  print(roc_combined_reload)
  
  # 
  cat("重新加载的AUC数据：\n")
  print(auc_data_reload)
}

# 8.4 
if (nrow(surv_data_reload) > 0) {
  cat("\n=== 从保存的数据重新绘制生存曲线 ===\n")
  
  survplot_list_reload <- list()
  survtable_list_reload <- list()
  clusters_surv <- unique(surv_data_reload$Cluster)
  
  for (cluster in clusters_surv) {
    # 
    cluster_surv_data <- surv_data_reload[surv_data_reload$Cluster == cluster, ]
    
    # 
    median_score <- median(cluster_surv_data$GSVA_Score)
    cluster_surv_data$GSVA_Group <- ifelse(cluster_surv_data$GSVA_Score >= median_score,
                                           "High GSVA", "Low GSVA")
    group_title <- paste0(cluster, "\n(Median Split)")
    
    # 
    fit <- survfit(Surv(Time, Status) ~ GSVA_Group, data = cluster_surv_data)
    
    # 
    surv_plot <- ggsurvplot(
      fit, data = cluster_surv_data,
      pval = TRUE, conf.int = FALSE, risk.table = TRUE,
      palette = c("#E74C3C", "#3498DB"),
      legend.labs = c("High GSVA", "Low GSVA"),
      legend.title = "GSVA Group",
      xlab = "Overall Survival (Days)",
      ylab = "Survival Probability",
      ggtheme = theme_bw(base_size = 11),
      title = group_title,
      font.title = c(12, "bold"),
      censor = TRUE, censor.size = 3, censor.shape = "|",
      risk.table.height = 0.25
    )
    
    survplot_list_reload[[cluster]] <- surv_plot$plot
    survtable_list_reload[[cluster]] <- surv_plot$table
  }
  
  # 
  if (length(survplot_list_reload) > 0) {
    surv_combined_reload <- wrap_plots(survplot_list_reload, ncol = 4) + 
      plot_annotation(
        title = paste0("Melanoma TLS GeneSets Survival Analysis (", surv_type, ")"),
        subtitle = paste("Group method:", group_method),
        theme = theme(
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12)
        )
      )
    
    legend_surv_reload <- get_legend(
      survplot_list_reload[[1]] + 
        theme(legend.position = "bottom") +
        guides(color = guide_legend(nrow = 1))
    )
    
    final_surv_plot_reload <- (surv_combined_reload / legend_surv_reload) + 
      plot_layout(heights = c(10, 1))
    
    print(final_surv_plot_reload)
    
    # ggsave(paste0(outpath, "F3_k_boxplot_reloaded.pdf"), box_combined_reload, width = 14, height = 6)
    # ggsave(paste0(outpath, "F3_k_roc_reloaded.pdf"), roc_combined_reload, width = 16, height = 4)
    # ggsave(paste0(outpath, "F3_k_survival_reloaded_", surv_type, "_", group_method, ".pdf"),
    #        final_surv_plot_reload, width = 16, height = 4)
  }
}

cat("\n=== 数据重新绘制完成 ===")
