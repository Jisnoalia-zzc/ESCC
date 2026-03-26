# Figure 4b -----
# ----------------Fig 4b1-----------
######################################

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(data.table)

pair_data_long <- fread('C:/Users/Ron/Desktop/Figs/F4/Adjusted_core_paired_long.csv')

# 绘图
p <- ggplot(pair_data_long, aes(x = Group, y = Score, fill = Group)) +
  geom_violin(alpha = 0.7, width = 1.2) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.8, outlier.shape = NA) +
  geom_line(aes(group = File), color = "gray80", alpha = 0.5, size = 0.8) +
  geom_point(aes(color = Group), size = 1, alpha = 0.6) +
  scale_fill_manual(values = c("MIMER" = "#E41A1C", "Permutation\n(mean of 1000 times)" = "#377EB8")) +
  scale_color_manual(values = c("MIMER" = "#E41A1C", "Permutation\n(mean of 1000 times)" = "#377EB8")) +
  theme_bw(base_size = 11) +
  labs(
    title = "Comparison of MIMER Scores: Original vs Permutation (TF Tissue Only)",
    subtitle = paste0("Wilcoxon signed-rank test, p = ", format(p_value, digits = 3, scientific = TRUE)),
    x = "",
    y = "Adjusted Score"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold", margin = margin(b = 10)),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray30", margin = margin(b = 15)),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold", color = "black", margin = margin(t = 5)),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.y = element_text(margin = margin(r = 10)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.margin = margin(20, 20, 20, 20)
  ) +
  stat_compare_means(
    paired = TRUE,
    method = "wilcox.test",
    label = "p.format",
    label.x = 1.5,
    label.y = y_max - (y_range * 0.1),
    size = 6,
    vjust = -0.5
  )

print(p)
#####################################
# ----------------Fig 4b2-----------
######################################
# Fig 4b
# 1. 准备数据
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

pair_data <- fread('C:/Users/Ron/Desktop/Figs/F4/Moran_pair_data_TF.csv')

# 2. 转换为长格式用于绘图
pair_data_long <- pair_data %>%
  pivot_longer(
    cols = c(MIMER, Permutation),
    names_to = "Group",
    values_to = "Moran_I"
  )

# 3. 计算Wilcoxon符号秩检验
wilcox_test <- wilcox.test(pair_data$MIMER, pair_data$Permutation, paired = TRUE)
p_value <- wilcox_test$p.value

# 4. 修改x轴标签，明确显示置换次数
pair_data_long <- pair_data_long %>%
  mutate(Group = factor(Group, 
                        levels = c("MIMER", "Permutation"),
                        labels = c("MIMER", "Permutation\n(mean of 1000 times)")),
         # 确保文件顺序一致
         File = factor(File, levels = unique(pair_data$File)))

# 5. 绘制配对小提琴图
p <- ggplot(pair_data_long, aes(x = Group, y = Moran_I, fill = Group)) +
  # 小提琴图
  geom_violin(alpha = 0.7,width=1.2) +
  # 白色箱线图
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.8, outlier.shape = NA) +
  # 配对连线
  geom_line(aes(group = File), 
            color = "gray80", alpha = 0.5, size = 0.8) +
  # 点
  geom_point(aes(color = Group), size = 1, alpha = 0.6) +
  # 颜色
  scale_fill_manual(values = c("MIMER" = "#E41A1C", "Permutation\n(mean of 1000 times)" = "#377EB8")) +
  scale_color_manual(values = c("MIMER" = "#E41A1C", "Permutation\n(mean of 1000 times)" = "#377EB8")) +
  # 主题
  theme_bw(base_size = 11) +  # 增加基础字体大小
  labs(
    title = "Comparison of Moran's I: MIMER vs Permutation",
    subtitle = paste0("Wilcoxon signed-rank test, p = ", format(p_value, digits = 3, scientific = TRUE)),
    x = "",
    y = "Moran's I"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold", margin = margin(b = 10)),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray30", margin = margin(b = 15)),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold", color = "black", margin = margin(t = 5)),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.y = element_text(margin = margin(r = 10)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.margin = margin(20, 20, 20, 20)
  )

# 6. 添加统计标注
# 计算y轴最大值用于放置p值标签
y_max <- max(pair_data_long$Moran_I, na.rm = TRUE)
y_min <- min(pair_data_long$Moran_I, na.rm = TRUE)
y_range <- y_max - y_min

p <- p + stat_compare_means(
  paired = TRUE,
  method = "wilcox.test",
  label = "p.format",
  label.x = 1.5,
  label.y = y_max - (y_range * 0.1),  # 在最高点上方10%处
  size = 6,
  vjust = -0.5
)

# 7. 显示图形
print(p)


# Figure 4e -----
load("cpdb_result.Rdata")
ggplot(cpdb,aes(x = LR,y = gene_pair, fill = mean,size = logPval)) +
  geom_point(shape=21,color='black') +
  scale_fill_gradientn(colors = color) +
  scale_y_discrete(limits=rev(levels(plot_df$subC))) +
  scale_size_continuous(range = c(1,8),name = 'Proportion') +
  facet_grid(~type,space = 'free_x',scale = 'free_x',labeller = label_wrap_gen(width=5)) +
  theme_bw() + xlab('') + ylab('') +
  theme(axis.text.x = element_text(angle = 60,hjust = 1),axis.text = element_text(size = 14,color = 'black'),strip.text.x = element_text(size = 10))
ligand_receptor_df = data.frame(ligand  = c("WNT5A","VEGFA","VEGFA","TNFSF4","TIMP1","TGFB1","SPP1","SPP1","SPP1","SPP1","SPP1","PVR","PDCD1","NECTIN2","MDK","MDK","LGALS9","IGF2","FN1","FN1","CXCL8","COL5A1",'COL4A2','COL4A1',"COL3A1","COL1A2",'COL1A1','COL1A1','COL1A1',"ANXA1"),
                                receptor = c("ANTXR1","KDR",'FLT1','TNFRSF4','FGFR2','TGFBR1',"PTGER4","CD44",'CCR8','ITGA9','ITGA4',"TIGIT","FAM3C",'TIGIT','SORL1','LRP1','HAVCR2','IGF2R','ITGA5','ITGA4','ACKR1','ITGA2','ITGA2','ITGA2','ITGA2','ITGA2','ITGA2','ITGA1','ITGA10',"FPR1"),
                                receptor2 = c("","","","","","TGFBR2","","","","ITGB1","ITGB1",'','','','',"","","","ITGB1","ITGB1",'',"ITGB1","ITGB1","ITGB1","ITGB1","ITGB1","ITGB1","ITGB1","ITGB1",""))
s1 = subset(sp,development =='unknown',invert=T)
  # slide = df[df$file == p,"image"]
  # # sample = df[df$file == p,"image"]
  # s1@images = s1@images[slide]
  
  ST.exp = s1@assays$spatial@data %>% as.data.frame()
  spot_abundance_category = s1@meta.data %>% as.matrix()
  
  ligand_receptor_df1 =  ligand_receptor_df %>% filter(receptor2 == "")
  ligand_receptor_df1$receptor2 = NULL
  
  LigandR_mean.ls = lapply(1:nrow(ligand_receptor_df1), function(i){
    ct = ligand_receptor_df1[i,] %>% as.character()
    tmp = apply(ST.exp[ct, ], 2, mean)
    return(tmp)
  })
  ligand_receptor_df1$lrpair <- paste0(ligand_receptor_df1$ligand, "_", ligand_receptor_df1$receptor)
  LigandR_mean.m <- do.call(rbind, LigandR_mean.ls)
  ligandReceptor.ls <- split(ligand_receptor_df1, seq(nrow(ligand_receptor_df1)))
  rownames(LigandR_mean.m) <- unlist(lapply(ligandReceptor.ls, function(x) paste0(x[1], "_", x[2]))) # step 2
  rownames(LigandR_mean.m) -> turn1
  avgExp <- as.data.frame(t(LigandR_mean.m)) %>% rownames_to_column(var = "cellID") # step 3
  data.plot <- left_join(avgExp, as.data.frame(spot_abundance_category[, c("mimer_score_spot")]) %>% rownames_to_column(var = "cellID"), by = "cellID")
  
  colnames(data.plot)[ncol(data.plot)] <- "mimer_score_spot"
  data_long1 <- reshape2::melt(data.plot, id.vars = c("cellID", "mimer_score_spot"), measure.vars = ligand_receptor_df1$lrpair, variable.name = "LigandReceptor", value.name = "LRmean")
  data_long1[, "mimer_score_spot"] <- factor(data_long1[, "mimer_score_spot"])
  
  
  
  ligand_receptor_df2 =  ligand_receptor_df %>% filter(receptor2 != "")
  # ligand_receptor_df2$receptor2 = NULL
  
  LigandR_mean.ls = lapply(1:nrow(ligand_receptor_df2), function(i){
    ct = ligand_receptor_df2[i,] %>% as.character()
    tmp = apply(ST.exp[ct, ], 2, mean)
    
    return(tmp)
  })
  ligand_receptor_df2$lrpair <- paste0(ligand_receptor_df2$ligand, "_", ligand_receptor_df2$receptor,"_",ligand_receptor_df2$receptor2)
  LigandR_mean.m <- do.call(rbind, LigandR_mean.ls)
  ligandReceptor.ls <- split(ligand_receptor_df2, seq(nrow(ligand_receptor_df2)))
  rownames(LigandR_mean.m) <- unlist(lapply(ligandReceptor.ls, function(x) paste0(x[1], "_", x[2],"_",x[3]))) # step 2
  rownames(LigandR_mean.m) -> turn2
  avgExp <- as.data.frame(t(LigandR_mean.m)) %>% rownames_to_column(var = "cellID") # step 3
  data.plot <- left_join(avgExp, as.data.frame(spot_abundance_category[, c("mimer_score_spot")]) %>% rownames_to_column(var = "cellID"), by = "cellID")
  
  colnames(data.plot)[ncol(data.plot)] <- "mimer_score_spot"
  data_long2 <- reshape2::melt(data.plot, id.vars = c("cellID", "mimer_score_spot"), measure.vars = ligand_receptor_df2$lrpair, variable.name = "LigandReceptor", value.name = "LRmean")
  data_long2[, "mimer_score_spot"] <- factor(data_long2[, "mimer_score_spot"])
  
  data_long = rbind(data_long1,data_long2)
  turn = gsub("\\_$","",
       paste0(ligand_receptor_df$ligand,"_",ligand_receptor_df$receptor,"_",ligand_receptor_df$receptor2))
  data_long$LigandReceptor = factor(data_long$LigandReceptor,
                                    levels = turn)
  g <- ggplot(data_long,aes(x=LigandReceptor,y= LRmean,fill=mimer_score_spot))+ 
    geom_boxplot(width = 0.4, outlier.shape = NA, lwd = 0.3, position = position_dodge(width = 0.7)) +
    # facet_grid(~LigandReceptor) + 
    scale_fill_manual(values = c("mimer"= 'brown','other' = 'grey'))+
    theme_bw() +
    xlab("") +
    ylab("Normalized ligand-receptor\naverage co-expression") +
    ggtitle("") +
    stat_compare_means(aes(group=mimer_score_spot),label.y = 3,show.legend = F,label = "p.signif",size=2)+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5, size = 8),
          panel.grid = element_blank(),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 8))+ylim(0,3.5)
  g
s1 = subset(sp,ME =='unknown',invert=T)
  slide = df[df$file == p,"image"]
  # # sample = df[df$file == p,"image"]
  # s1@images = s1@images[slide]
  
  ST.exp = s1@assays$spatial@data %>% as.data.frame()
  spot_abundance_category = s1@meta.data %>% as.matrix()
  
  ligand_receptor_df1 =  ligand_receptor_df %>% filter(receptor2 == "")
  ligand_receptor_df1$receptor2 = NULL
  
  LigandR_mean.ls = lapply(1:nrow(ligand_receptor_df1), function(i){
    ct = ligand_receptor_df1[i,] %>% as.character()
    tmp = apply(ST.exp[ct, ], 2, mean)
    return(tmp)
  })
  ligand_receptor_df1$lrpair <- paste0(ligand_receptor_df1$ligand, "_", ligand_receptor_df1$receptor)
  LigandR_mean.m <- do.call(rbind, LigandR_mean.ls)
  ligandReceptor.ls <- split(ligand_receptor_df1, seq(nrow(ligand_receptor_df1)))
  rownames(LigandR_mean.m) <- unlist(lapply(ligandReceptor.ls, function(x) paste0(x[1], "_", x[2]))) # step 2
  avgExp <- as.data.frame(t(LigandR_mean.m)) %>% rownames_to_column(var = "cellID") # step 3
  data.plot <- left_join(avgExp, as.data.frame(spot_abundance_category[, c("mimer_score_spot")]) %>% rownames_to_column(var = "cellID"), by = "cellID")
  
  colnames(data.plot)[ncol(data.plot)] <- "mimer_score_spot"
  data_long1 <- reshape2::melt(data.plot, id.vars = c("cellID", "mimer_score_spot"), measure.vars = ligand_receptor_df1$lrpair, variable.name = "LigandReceptor", value.name = "LRmean")
  data_long1[, "mimer_score_spot"] <- factor(data_long1[, "mimer_score_spot"])
  
  
  
  ligand_receptor_df2 =  ligand_receptor_df %>% filter(receptor2 != "")
  # ligand_receptor_df2$receptor2 = NULL
  
  LigandR_mean.ls = lapply(1:nrow(ligand_receptor_df2), function(i){
    ct = ligand_receptor_df2[i,] %>% as.character()
    tmp = apply(ST.exp[ct, ], 2, mean)
    return(tmp)
  })
  ligand_receptor_df2$lrpair <- paste0(ligand_receptor_df2$ligand, "_", ligand_receptor_df2$receptor,"_",ligand_receptor_df2$receptor2)
  LigandR_mean.m <- do.call(rbind, LigandR_mean.ls)
  ligandReceptor.ls <- split(ligand_receptor_df2, seq(nrow(ligand_receptor_df2)))
  rownames(LigandR_mean.m) <- unlist(lapply(ligandReceptor.ls, function(x) paste0(x[1], "_", x[2],"_",x[3]))) # step 2
  avgExp <- as.data.frame(t(LigandR_mean.m)) %>% rownames_to_column(var = "cellID") # step 3
  data.plot <- left_join(avgExp, as.data.frame(spot_abundance_category[, c("mimer_score_spot")]) %>% rownames_to_column(var = "cellID"), by = "cellID")
  
  colnames(data.plot)[ncol(data.plot)] <- "mimer_score_spot"
  data_long2 <- reshape2::melt(data.plot, id.vars = c("cellID", "mimer_score_spot"), measure.vars = ligand_receptor_df2$lrpair, variable.name = "LigandReceptor", value.name = "LRmean")
  data_long2[, "mimer_score_spot"] <- factor(data_long2[, "mimer_score_spot"])
  
  data_long = rbind(data_long1,data_long2)
  
  turn = gsub("\\_$","",
              paste0(ligand_receptor_df$ligand,"_",ligand_receptor_df$receptor,"_",ligand_receptor_df$receptor2))
  data_long$LigandReceptor = factor(data_long$LigandReceptor,
                                    levels = turn)
  
  g <- ggplot(data_long,aes(x=LigandReceptor,y= LRmean,fill=mimer_score_spot))+ 
    geom_boxplot(width = 0.4, outlier.shape = NA, lwd = 0.3, position = position_dodge(width = 0.7)) +
    # facet_grid(~LigandReceptor) + 
    scale_fill_manual(values = c("mimer"= 'brown','other' = 'grey'))+
    theme_bw() +
    xlab("") +
    ylab("Normalized ligand-receptor\naverage co-expression") +
    ggtitle("") +
    stat_compare_means(aes(group=mimer_score_spot),label.y = 3,show.legend = F,label = "p.signif",size=2)+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5, size = 8),
          panel.grid = element_blank(),
          axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 8))+ylim(0,3.5)
  g

#------------ Fig 4d---------
##############################
library(dplyr)
library(ggplot2)
library(dittoSeq)  # 添加dittoSeq包用于配色
library(tidyr)     # 用于数据整理

mimer_cc <- fread('C:/Users/Ron/Desktop/Figs/F4/mimer_cc.csv')
# 1. 数据预处理
cell_table <- mimer_cc %>%
  select(Slide_ID = file, CN_Label = cc_15, Is_MIMER = mimer) %>%
  mutate(Is_MIMER = ifelse(Is_MIMER == "MIMER", 1, 0))

# 2. 新增第一步：样本级过滤 —— 排除整个样本中MIMER+总数少于10的样本
cat("原始数据行数:", nrow(cell_table), "\n")
cat("原始样本数:", length(unique(cell_table$Slide_ID)), "\n")

# 计算每个样本的MIMER+细胞总数
sample_mimer_counts <- cell_table %>%
  group_by(Slide_ID) %>%
  summarise(total_mimer_in_sample = sum(Is_MIMER), .groups = "drop")

# 识别MIMER+总数 >= 10 的有效样本
valid_samples <- sample_mimer_counts %>%
  filter(total_mimer_in_sample >= 10) %>%
  select(Slide_ID)

# 首次过滤：只保留有效样本中的数据
cell_table_filtered_by_sample <- cell_table %>%
  inner_join(valid_samples, by = "Slide_ID")

# 3. 第二步：细胞类型-样本组合过滤 —— 只保留细胞类型计数≥10的组合
# 计算每个样本中每种细胞类型的数量
cell_counts <- cell_table_filtered_by_sample %>%
  group_by(Slide_ID, CN_Label) %>%
  summarise(cn_count = n(), .groups = "drop")

# 识别有效的组合
valid_combinations <- cell_counts %>%
  filter(cn_count >= 10) %>%
  select(Slide_ID, CN_Label)

# 第二次过滤：只保留有效组合的数据
cell_table_double_filtered <- cell_table_filtered_by_sample %>%
  inner_join(valid_combinations, by = c("Slide_ID", "CN_Label"))

# 重新获取过滤后的细胞类型列表
all_ccs <- unique(cell_table_double_filtered$CN_Label) %>% sort()

# 使用dittoColors创建颜色映射
color_mapping <- setNames(dittoColors()[1:length(all_ccs)], all_ccs)

# 4. 核心分析函数（精简版）- 使用双重过滤后的数据
run_mimer_analysis <- function(cell_table, n_perm = 50) {
  # 计算观测值
  obs_stats <- cell_table %>%
    group_by(Slide_ID, CN_Label) %>%
    summarise(
      cn_size = n(),
      obs_mimer_count = sum(Is_MIMER == 1),
      obs_prop = obs_mimer_count / cn_size,
      .groups = "drop"
    )
  
  # 置换检验
  set.seed(123)
  perm_results <- lapply(1:n_perm, function(i) {
    cell_table %>%
      group_by(Slide_ID) %>%
      mutate(Is_MIMER = sample(Is_MIMER)) %>%
      group_by(Slide_ID, CN_Label) %>%
      summarise(perm_prop = sum(Is_MIMER == 1) / n(), .groups = "drop")
  }) %>% bind_rows()
  
  # 统计显著性
  final_stats <- perm_results %>%
    group_by(Slide_ID, CN_Label) %>%
    summarise(
      null_mean = mean(perm_prop),
      null_sd = sd(perm_prop),
      .groups = "drop"
    ) %>%
    left_join(obs_stats, by = c("Slide_ID", "CN_Label")) %>%
    rowwise() %>%
    mutate(
      p_perm = (sum(perm_results$perm_prop[perm_results$CN_Label == CN_Label & 
                                             perm_results$Slide_ID == Slide_ID] >= obs_prop) + 1) / (n_perm + 1),
      z_score = (obs_prop - null_mean) / (null_sd + 1e-6),
      log2FC = log2((obs_prop + 0.001) / (null_mean + 0.001))
    ) %>%
    ungroup() %>%
    mutate(p_adj = p.adjust(p_perm, method = "BH"))
  
  return(final_stats)
}

# 5. 修改后的可视化函数 - 按中位数从大到小排序
plot_results <- function(results, color_mapping) {
  # 计算每个细胞类型的z_score中位数，用于排序
  zscore_medians <- results %>%
    group_by(CN_Label) %>%
    summarise(median_z = median(z_score, na.rm = TRUE)) %>%
    arrange(median_z)  # 按中位数从大到小排序
  
  # 按中位数排序细胞类型
  results$CN_Label <- factor(results$CN_Label, levels = zscore_medians$CN_Label)
  
  # 重新排序颜色映射以匹配因子水平
  ordered_colors <- color_mapping[levels(results$CN_Label)]
  
  # 总体汇总图 - 按z_score中位数排序
  p1 <- ggplot(results, aes(x = CN_Label, y = z_score)) +
    geom_boxplot(aes(fill = CN_Label), outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
    geom_hline(yintercept = 0, linetype = "solid", color = "#000000", linewidth = 0.5) +
    geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed",
               color = c("blue", "red"), alpha = 0.7) +
    scale_fill_manual(values = ordered_colors) +  # 使用重新排序的颜色
    coord_flip() +
    theme_classic() +
    labs(x = "CC_15", y = "Enrichment Z-score", 
         subtitle = "Double Filtered: Sample MIMER+ >= 10 & CC per sample >= 10\nSorted by median z-score") +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      plot.subtitle = element_text(size = 9, color = "gray50", hjust = 0.5)
    )
  
  # 显著性气泡图 - 也按相同的顺序排序
  p2 <- ggplot(results, aes(x = CN_Label, y = Slide_ID)) +
    geom_point(aes(size = -log10(p_adj), color = z_score)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    scale_x_discrete(limits = levels(results$CN_Label)) +  # 使用相同的排序
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, 
                                 color = ordered_colors[levels(results$CN_Label)]),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    ) +
    labs(x = "CN Label", y = "Sample", color = "Z-score", size = "-log10(FDR)",
         subtitle = "Double Filtered: Sample MIMER+ >= 10 & CC per sample >=10")
  
  return(list(summary_plot = p1, bubble_plot = p2, 
              zscore_medians = zscore_medians, ordered_colors = ordered_colors))
}

# 6. MIMER点占CC点比例的横向箱线图 - 按中位数从大到小排序
plot_prop_distribution <- function(results, color_mapping) {
  # 创建绘图数据
  plot_data <- results %>%
    select(Slide_ID, CN_Label, obs_prop, cn_size, obs_mimer_count) %>%
    distinct()  # 确保每个样本-细胞类型组合只出现一次
  
  # 计算每个细胞类型的观测比例中位数，用于排序
  median_props <- plot_data %>%
    group_by(CN_Label) %>%
    summarise(median_prop = median(obs_prop, na.rm = TRUE)) %>%
    arrange(median_prop)  # 按中位数从大到小排序
  
  # 按照中位数排序细胞类型
  plot_data$CN_Label <- factor(plot_data$CN_Label, 
                               levels = median_props$CN_Label)
  
  # 创建颜色向量，确保与排序后的细胞类型对应
  ordered_colors <- color_mapping[levels(plot_data$CN_Label)]
  
  # 创建横向箱线图
  p3 <- ggplot(plot_data, aes(x = CN_Label, y = obs_prop)) +
    geom_boxplot(aes(fill = CN_Label), outlier.shape = NA, width = 0.7) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1, color = "gray30") +
    scale_fill_manual(values = ordered_colors) +
    scale_y_continuous(
      limits = c(-0.0001, 1.0001),
      expand = expansion(mult = c(0, 0.05))
    ) +
    coord_flip() +
    theme_classic() +
    labs(
      x = "Cell Type (CC_15)",
      y = "Proportion of MIMER+ cells",
      title = "MIMER+ Cell Proportion by Cell Type",
      subtitle = paste0("Double Filtered: Sample MIMER+ ≥10 & CC per sample ≥10\n")
    ) +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray50"),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor.x = element_line(color = "grey95", linewidth = 0.2)
    )
  
  return(list(plot = p3, median_props = median_props, ordered_colors = ordered_colors))
}

# 7. 运行分析 - 使用双重过滤后的数据
results <- run_mimer_analysis(cell_table_double_filtered, n_perm = 1000)

# 8. 查看显著结果
significant_results <- results %>% 
  filter(z_score > 1.96, p_adj < 0.05) %>%
  arrange(desc(z_score))

if (nrow(significant_results) > 0) {
  print(significant_results)
} else {
  cat("没有发现显著富集的结果。\n")
}

# 设置图形参数
options(repr.plot.width = 8, repr.plot.height = 6)

# 10. 绘图
plots <- plot_results(results, color_mapping)
print(plots$summary_plot)  # 显示汇总图


# Figure 4d -----
# ==========================================
# 最终确认版：按分类打分并配对绘图脚本
# ==========================================
load('~/workspace/proj_ESCC_STW_ZWM_2022_01/liuliqiu/final_analysis/codex_datasets.Rdata')
load('~/workspace/proj_ESCC_STW_ZWM_2022_01/liuliqiu/final_analysis/classification_results.Rdata')
codex
head(classification_results)
mimer_mapping <- data.frame(
  MIMER = c(
    "CD4_C7_0X40", "CD4_C7_0X40",
    "CD8_C6_CD39", "CD8_C6_CD39", "CD8_C6_CD39",
    "Endo_C3_RGCC", "Endo_C3_RGCC",
    "FB_C3_COL1A1", "FB_C3_COL1A1",
    "Mac_C2_SPP1", "Mac_C2_SPP1"
  ),
  CODEX = c(
    "CD4T-aTreg", "Tcells.cc",
    "PDCD1+Tex", "LAG3+Tex", "PDCD1+Tex2",
    "PLVAP+Endo", "CollagenIV+Endo",
    "COL1A1+CAF", "POSTN+CAF",
    "ISG15+Mac", "SPP1+Mac"
  )
)
colnames(codex@meta.data)
sub_codex = subset(codex,subCelltype.x %in% mimer_mapping$CODEX )
mimer_id = classification_results %>% filter(final_classification == 'MIMER_group')
sub_codex$Mimer_group = sub_codex$subCelltype.x
# 更简洁的单行写法
sub_codex$Mimer_group <- ifelse(
  rownames(sub_codex@meta.data) %in% mimer_id$cell_id,
  paste0(sub_codex$Mimer_group, '_MIMER'),
  sub_codex$subCelltype.x
)
table(sub_codex$Mimer_group)
save(sub_codex,file = 'data_for_codex_score.Rdata')
# --- 1. 定义映射关系 (建立细胞大类与通路的专属纽带) ---
# 确保 Tex 相关的得分只在 Tex 类细胞中计算
pathway_defs <- list(
  # CD8_C6_CD39 相关 (Tex)
  Tex_ExhaustionScore = c("PD-1", "TOX", "LAG3"),
  Tex_TCF1 = c("TCF-1"),
  Tex_EffectorScore = c("GZMB", "IFNG", "CD57"),
  Tex_ProliferationScore = c("Ki67", "PCNA", "HistoneH3-pSer28"),
  
  # CD4_C7_0X40 相关 (Treg)
  Treg_CD39 = c("CD39"),
  Treg_ActivationScore = c("ICOS", "OX40"),
  Treg_ProliferationScore = c("Ki67", "PCNA", "HistoneH3-pSer28"),
  
  # Mac_C2_SPP1 相关 (TAM)
  TAM_ImmunosuppressionScore = c("PD-L1", "IDO1", "VISTA", "CD39"),
  TAM_M2LikeScore = c("CD163", "CD206"),
  TAM_AntigenPresentationScore = c("HLA-DR"),
  
  # Endo_C3_RGCC 相关 (Endo)
  Endo_PermeabilityScore = c("PLVAP-PV-1", "Caveolin"),
  Endo_ImmuneRegulationScore = c("CD39", "PD-L1", "IDO1", "HLA-DR"),
  Endo_ProliferationScore = c("Ki67", "PCNA", "HistoneH3-pSer28"),
  
  # FB_C3_COL1A1 相关 (CAF)
  CAF_SMA = c("SMA"),
  CAF_ECMRemodelingScore = c("Periostin", "COL-1", "CollagenIV"),
  CAF_ImmuneRegulationScore = c("PD-L1", "IDO1"),
  Stromal_ISG15_ResponseScore = c("ISG15")
)


mimer_to_axes_list <- list(
  "CD8_C6_CD39"  = c("Tex_ExhaustionScore", "Tex_TCF1", "Tex_EffectorScore", "Tex_ProliferationScore"),
  "CD4_C7_0X40"  = c("Treg_CD39", "Treg_ActivationScore", "Treg_ProliferationScore"),
  "Mac_C2_SPP1"  = c("TAM_ImmunosuppressionScore", "TAM_M2LikeScore", "TAM_AntigenPresentationScore"),
  "Endo_C3_RGCC" = c("Endo_PermeabilityScore", "Endo_ImmuneRegulationScore", "Endo_ProliferationScore"),
  "FB_C3_COL1A1" = c("CAF_SMA", "CAF_ECMRemodelingScore", "CAF_ImmuneRegulationScore", "Stromal_ISG15_ResponseScore")
)

# --- 2. 核心打分函数 (只在对应分类下打分) ---
run_exclusive_scoring <- function(obj, mapping_df, axes_map, pathway_defs, cofactor = 5) {
  cat("Step 1: Preparing data and mapping...\n")
  
  # 提取元数据并建立基础亚型与 Condition
  meta <- obj@meta.data[, c("Patient", "Mimer_group")]
  meta$base_mimer <- gsub("_MIMER$", "", meta$Mimer_group)
  meta$condition <- ifelse(grepl("_MIMER$", meta$Mimer_group), "MIMER", "Original")
  
  # 关联 MIMER 大类标签 (如 CD8_C6_CD39)
  # 假设 mapping_df 包含 mimer (亚型) 和 MIMER (大类) 两列
  meta <- merge(meta, mapping_df, by.x = "base_mimer", by.y = "mimer", all.x = TRUE)
  
  # 提取表达矩阵 (背景扣除与 Arcsinh)
  all_markers <- unique(unlist(pathway_defs))
  expr_mat <- GetAssayData(obj, assay = "Akoya", layer = "data")
  available_markers <- intersect(all_markers, rownames(expr_mat))
  
  data_mat <- t(as.matrix(expr_mat[available_markers, ]))
  df_full <- cbind(as.data.frame(data_mat), meta)
  
  # 按 Patient 做背景校正
  df_list <- split(df_full, df_full$Patient)
  df_processed <- do.call(rbind, lapply(df_list, function(sub_df) {
    for (m in available_markers) {
      bg <- quantile(sub_df[[m]], 0.01, na.rm = TRUE)
      sub_df[[m]] <- asinh(pmax(0, sub_df[[m]] - bg) / cofactor)
    }
    return(sub_df)
  }))
  
  # --- 关键：只在同一个细胞特征下进行打分 ---
  cat("Step 2: Calculating category-exclusive scores...\n")
  final_results <- list()
  
  for (m_category in names(axes_map)) {
    # 只筛选属于该大类的细胞 (例如所有 Tex 亚型)
    cat_cells <- df_processed[df_processed$MIMER == m_category & !is.na(df_processed$MIMER), ]
    if (nrow(cat_cells) == 0) next
    
    target_axes <- axes_map[[m_category]]
    for (axis in target_axes) {
      markers <- intersect(pathway_defs[[axis]], colnames(cat_cells))
      if (length(markers) > 0) {
        # 计算该大类下每个细胞的得分
        cat_cells[[axis]] <- rowMeans(cat_cells[, markers, drop = FALSE], na.rm = TRUE)
        
        # 准备聚合数据：保留原始亚型名称 base_mimer
        # 使用 aggregate 避免 tidyverse 报错
        agg_data <- aggregate(formula(paste(axis, "~ Patient + condition + base_mimer + MIMER")), 
                              data = cat_cells, FUN = mean, na.rm = TRUE)
        colnames(agg_data)[ncol(agg_data)] <- "mean_score"
        agg_data$Pathway <- axis
        final_results[[paste0(m_category, "_", axis)]] <- agg_data
      }
    }
  }
  return(do.call(rbind, final_results))
}

# --- 4. 运行实例 ---
scored_plot_df <- run_exclusive_scoring(sub_codex, mimer_mapping2, mimer_to_axes_list, pathway_defs)

# ==========================================
# 针对已聚合的 sample-level 数据进行配对绘图
# ==========================================

plot_exclusive_results_direct <- function(plot_df, output_dir = "plots_pformat") {
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  # 强制统一列名，确保引用正确
  # 假设列顺序为: Patient, condition, base_mimer, MIMER, mean_score, Pathway
  # 如果您的列名已经是这些，可以跳过重命名
  
  library(ggplot2)
  library(ggpubr)
  
  # 1. 遍历每个大类 (MIMER)
  unique_mimers <- unique(plot_df$MIMER)
  
  for (m_cat in unique_mimers) {
    # 提取当前大类的数据
    sub_data <- plot_df[plot_df$MIMER == m_cat & !is.na(plot_df$MIMER), ]
    
    # 2. 关键：筛选配对病人 (必须在同一个 Pathway + base_mimer 下同时有 Original 和 MIMER)
    # 即使是 sample-level，也需要剔除那些没有对照的“孤儿样本”
    # 使用 base R 的 split-apply 逻辑检查配对性
    check_pair <- aggregate(condition ~ Patient + base_mimer + Pathway, 
                            data = sub_data, FUN = length)
    paired_keys <- check_pair[check_pair$condition == 2, c("Patient", "base_mimer", "Pathway")]
    
    # 执行过滤
    sub_data_paired <- merge(sub_data, paired_keys, by = c("Patient", "base_mimer", "Pathway"))
    
    if (nrow(sub_data_paired) == 0) {
      cat("Skipping Category:", m_cat, "- No paired samples found.\n")
      next
    }
    
    # 3. 绘图
    p <- ggplot(sub_data_paired, aes(x = condition, y = mean_score, color = condition)) +
      # 连线：连接同一个 Patient 的两组数据
      geom_line(aes(group = Patient), color = "gray80", size = 0.5) +
      geom_boxplot(aes(fill = condition), 
                   width = 0.3,            # 减小宽度（例如从 0.5 调到 0.3），柱子间距会显著增大
                   alpha = 0.2, 
                   outlier.shape = NA,
                   position = position_dodge(0.8)) +
      geom_point(size = 2, alpha = 0.8) +
      # 分面：横向为原始亚型 (base_mimer)，纵向为该大类的专属通路 (Pathway)
      facet_grid(Pathway ~ base_mimer, scales = "free_y") + 
      # 配对 Wilcoxon 检验
      stat_compare_means(
        paired = TRUE,
        method = "wilcox.test", label = "p.format") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      labs(title = paste("Category Exclusive Analysis:", m_cat),
           subtitle = "Sample-level paired comparison",
           y = "Mean Functional Score", x = "") +
      scale_color_manual(values = c("other" = "#66c2a5", "MIMER" = "#fc8d62")) +
      scale_fill_manual(values = c("other" = "#66c2a5", "MIMER" = "#fc8d62"))
    
    # 保存
    ggsave(file.path(output_dir, paste0(m_cat, "_exclusive_paired_plot.pdf")), 
           p, width = 10, height = 10)
    
    cat("Successfully generated plot for category:", m_cat, "\n")
  }
}

# --- 调用 ---
head(scored_plot_df)
table(scored_plot_df$con)
scored_plot_df$condition = gsub('Original','other',
                                scored_plot_df$condition)
plot_exclusive_results_direct(scored_plot_df)
scored_plot_df %>% write.csv(.,file = 'codex_function_data.csv',quote = F)
head(scored_plot_df)

scored_plot_df$Mimer_group = ifelse(scored_plot_df$condition == 'MIMER',
                                    paste0(scored_plot_df$base_mime,"_MIMER"),
                                    scored_plot_df$base_mime)
head(scored_plot_df)
sub_codex@meta.data %>%
  dplyr::group_by(Patient, Mimer_group) %>%
  dplyr::summarise(number = n()) %>%
  dplyr::left_join(scored_plot_df,by = c('Patient','Mimer_group')) %>%
  write.csv(.,file = 'codex_function_data_with_count.csv',quote = F)



load('~/workspace/proj_ESCC_STW_ZWM_2022_01/codex/codex_final.Rdata')
# 载入必要的库
library(tidyverse)
library(ggpubr)
library(reshape2)

# --- 1. 数据模拟准备 (请将其替换为你自己的 Cell Table) ---
# 假设数据框名为 cell_table
# 必须包含列: Slide_ID, CN_Label, Is_MIMER (1为是, 0为否)
# 示例：cell_table <- read.csv("your_codex_data.csv")
cell_table = codex2@meta.data %>%
  dplyr::select(-cell_id) %>%
  tibble::rownames_to_column('cell_id') %>%
  dplyr::select(cell_id,neighborhood10_new) %>%
  dplyr::left_join(classification_results,by = 'cell_id')
# 确保已经安装了 dplyr, tibble, ggplot2 (绘图用)
# install.packages(c("dplyr", "tibble", "ggplot2"))
cell_table$Is_MIMER = ifelse(cell_table$final_classification =='MIMER_group',1,0)
cell_table$CN_Label = cell_table$neighborhood10_new
cell_table$Slide_ID = cell_table$Sample
# 核心分析流程函数
# cell_table: 包含 Slide_ID, CN_Label, Is_MIMER (0/1) 的数据框
# 核心分析流程函数
# cell_table: 包含 Slide_ID, CN_Label, Is_MIMER (0/1) 的数据框
run_mimer_spatial_analysis <- function(cell_table, n_perm = 1000) {
  
  # --- Step A: 计算观测到的真实分布 ---
  message("正在计算原始观测数据...")
  obs_stats <- cell_table %>%
    dplyr::group_by(Slide_ID, CN_Label) %>%
    dplyr::summarise(
      cn_size = dplyr::n(),
      obs_mimer_count = sum(Is_MIMER == 1, na.rm = TRUE),
      obs_prop = obs_mimer_count / cn_size,
      .groups = "drop"
    )
  
  # --- Step B: 执行置换检验 (Permutation Test) ---
  message(paste0("正在执行 ", n_perm, " 次置换检验..."))
  
  # 预分配存储列表
  perm_results_list <- vector("list", n_perm)
  
  set.seed(42) # 保证结果可重复
  for(i in 1:n_perm) {
    # 核心逻辑：保持空间结构(CN)不变，在每个Slide内部随机打乱MIMER标签
    perm_data <- cell_table %>%
      dplyr::group_by(Slide_ID) %>%
      dplyr::mutate(Is_MIMER = sample(Is_MIMER)) %>%
      dplyr::group_by(Slide_ID, CN_Label) %>%
      dplyr::summarise(perm_prop = sum(Is_MIMER == 1) / dplyr::n(), .groups = "drop")
    
    perm_results_list[[i]] <- perm_data
    if(i %% 100 == 0) message(paste("已完成:", i, "/", n_perm))
  }
  
  # 合并所有随机结果
  all_perms <- dplyr::bind_rows(perm_results_list)
  
  # --- Step C: 统计显著性计算 ---
  message("正在计算统计显著性指标...")
  final_stats <- all_perms %>%
    dplyr::group_by(Slide_ID, CN_Label) %>%
    dplyr::summarise(
      null_mean = mean(perm_prop),
      null_sd = sd(perm_prop),
      # 经验P值：随机分布中比例大于等于观测比例的频率
      .groups = "drop"
    ) %>%
    # 合并观测值
    dplyr::left_join(obs_stats, by = c("Slide_ID", "CN_Label")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      # 计算 P 值 (基于置换分布)
      p_perm = (sum(all_perms$perm_prop[all_perms$CN_Label == CN_Label & 
                                          all_perms$Slide_ID == Slide_ID] >= obs_prop) + 1) / (n_perm + 1)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      # 多重检验校正
      p_adj = stats::p.adjust(p_perm, method = "BH"),
      # Z-score: 衡量偏离随机分布的程度
      z_score = (obs_prop - null_mean) / (null_sd + 1e-6),
      # 富集倍数 (Log2 Fold Change)
      log2FC = log2((obs_prop + 0.001) / (null_mean + 0.001))
    )
  
  return(final_stats)
}

library(ggplot2)

# --- A. 单张片子详细富集图 (气泡图) ---
plot_mimer_single_slide <- function(results, slide_id) {
  df_sub <- results %>% dplyr::filter(Slide_ID == slide_id)
  
  ggplot(df_sub, aes(x = stats::reorder(CN_Label, z_score), y = z_score)) +
    geom_point(aes(size = -log10(p_adj), color = z_score)) +
    geom_hline(yintercept = 1.96, linetype = "dashed", color = "#E64B35") + # Z=1.96 对应 p=0.05
    geom_hline(yintercept = 0, color = "black") +
    scale_color_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0) +
    coord_flip() +
    theme_bw() +
    labs(title = paste("MIMER Enrichment in Sample:", slide_id),
         subtitle = "Dashed line: Significant enrichment (Z > 1.96)",
         x = "Cellular Neighborhood (CN)", y = "Enrichment Z-score",
         size = "-log10(FDR)", color = "Z-score")
}

# --- B. 所有片子的一致性总结图 (箱线图) ---
plot_mimer_summary_box <- function(results) {
  ggplot(results, aes(x = stats::reorder(CN_Label, z_score, FUN = median), y = z_score)) +
    geom_boxplot(aes(fill = CN_Label), outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.2, size = 1) +
    geom_hline(yintercept = 0, linetype = "solid") +
    geom_hline(yintercept = 1.96, linetype = "dotted", color = "red",linewidth = 1) +
    coord_flip() +
    theme_classic() +
    theme(legend.position = "none") +
    labs(title = "MIMER Consistency across All Samples",
         subtitle = "Each point is one Slide; Dashed line: Z = 1.96",
         x = "Cellular Neighborhood", y = "Enrichment Z-score")
}

# 假设你的数据叫 cell_metadata
# 1. 运行统计 (可能需要几分钟，取决于细胞量)
final_results <- run_mimer_spatial_analysis(cell_table, n_perm = 1000)

# 2. 保存结果
utils::write.csv(final_results, "MIMER_CN_Enrichment_Results.csv", row.names = FALSE)

# 3. 画一张特定片子的图
p1 <- plot_mimer_single_slide(final_results, "B01")
print(p1)

# 4. 画全样本总结图
cols_20 <- c("#4E79A7", "#A0CBE8", "#90ee90", "#FFBE7D", 
             "#f28e2b", "#8CD17D", "#B6992D", "#F1CE63", "#ff0000", 
             "#86BCB6", "#E15759", "#FF9D9A", "#79706E", "#BAB0AC", 
             "#D37295", "#FABFD2", "#B07AA1", "#D4A6C8", "#9D7660", 
             "#D7B5A6")
niche_cols = cols_20[1:16]
names(niche_cols) = paste0('CN',1:16)
p2 <- plot_mimer_summary_box(final_results)+
  scale_fill_manual(values = c(niche_cols,"other"='#eae8dc'))
print(p2)
ggplot2::ggsave("MIMER_Enrichment_total.pdf", p2, width = 6, height = 4)


plot_mimer_combined <- function(results) {
  # 确保 p_adj 不为 0，否则取 -log10 会报错，设置一个极小值
  results <- results %>%
    dplyr::mutate(log_p = -log10(p_adj + 1e-10))
  
  ggplot2::ggplot(results, 
                  ggplot2::aes(x = stats::reorder(CN_Label, z_score), y = z_score)) +
    # 核心：点的大小由显著性决定，颜色由 Z-score 决定
    ggplot2::geom_point(ggplot2::aes(size = log_p, color = z_score)) +
    
    # 显著性参考线 (Z = 1.96)
    ggplot2::geom_hline(yintercept = 1.96, linetype = "dashed", color = "#E64B35", alpha = 0.7) +
    ggplot2::geom_hline(yintercept = 0, color = "black", size = 0.3) +
    
    # 颜色梯度：蓝色(排斥) -> 白色 -> 红色(富集)
    ggplot2::scale_color_gradient2(low = "#3C5488FF", mid = "white", high = "#E64B35FF", midpoint = 0) +
    
    # 分面：按照 Slide_ID 自动拼图，ncol=3 表示每行画3个样
    ggplot2::facet_wrap(~Slide_ID, ncol = 10, scales = "free_y") +
    
    # 坐标轴翻转，让 CN 名称横向显示，易于阅读
    ggplot2::coord_flip() +
    
    # 样式美化
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "gray90"), # 标题栏颜色
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 8)
    ) +
    
    ggplot2::labs(
      title = "MIMER Enrichment Analysis Across Multiple Samples",
      subtitle = "Point size: -log10(Adjusted P-value); Red line: Z-score = 1.96 (Significance threshold)",
      x = "Cellular Neighborhood (CN)",
      y = "Enrichment Z-score",
      size = "-log10(FDR)",
      color = "Z-score"
    )
}

# 使用方法：
p_final <- plot_mimer_combined(final_results)
ggplot2::ggsave("MIMER_Enrichment_Combined.pdf", p_final, width = 40, height = 40)

# Fig 4F
# 加载必要的包
library(ggplot2)
library(scales)
library(data.table)
library(dplyr)
library(tidyr)  # 用于数据重塑

# 设置分析类型：'development' 或 'ME'
analysis_type <- "development"  # 可更改为 "development" 或 "ME"

# 根据分析类型处理数据
if (analysis_type == "ME") {
  # 使用ME列
  type_col <- "ME"
  all_types <- c("Nor-ME", "Hyp-ME", "MiD-ME", "MoD-ME", "SD&CA-ME", "ICA-ME", "MCA-ME")
  # 定义Nor组
  nor_group <- "Nor-ME"
} else {
  # 使用development列
  type_col <- "development"
  all_types <- c("Nor", "Hyp", "MiD", "MoD", "SD&CA", "ICA", "MCA")
  # 定义Nor组
  nor_group <- "Nor"
}

# 过滤掉unknown和MCA类型
filtered_types <- all_types[!all_types %in% c("unknown", 
                                              ifelse(analysis_type == "ME", "MCA-ME", "MCA"))]
cat(sprintf("\n分析类型: %s\n", analysis_type))
cat(sprintf("使用的列: %s\n", type_col))
cat("过滤后的类型顺序:\n")
print(filtered_types)

# 定义要使用的阈值类型
thresholds <- c("top5_perc_mimer", "top10_perc_mimer", "top20_perc_mimer")
threshold_names <- c("Top 5%", "Top 10%", "Top 20%")
threshold_colors <- c(
  "Top 5%" = "#FF7F0E",    # 橙色
  "Top 10%" = "#3fa02c",   # 绿色
  "Top 20%" = "#1F77B4"    # 蓝色
)

# 初始化存储结果
all_threshold_stats <- list()
comparison_results <- list()  # 存储比较结果

# 过滤掉unknown和MCA类型的数据
base_data <- fread('C:/Users/Ron/Desktop/Figs/F4/Development_ME_dashdot_base_data.csv')

# 计算每个样本每个类型的总细胞数
type_total_counts <- base_data[, 
                               .(total_cells = .N), 
                               by = .(file, get(type_col))
]
setnames(type_total_counts, "get", "type")

# 设置类型为因子
type_total_counts$type <- factor(type_total_counts$type, levels = filtered_types)

# 计算每个类型实际出现在哪些文件中
type_files <- base_data[, 
                        .(files = unique(file)), 
                        by = get(type_col)
]
setnames(type_files, "get", "type")

# 计算每个类型的文件数（样本数）
type_n_files <- type_files[, .(n_files = .N), by = type]
type_n_files$type <- factor(type_n_files$type, levels = filtered_types)
type_n_files <- type_n_files[order(type)]

cat("\n每个类型出现在多少个文件中:\n")
print(type_n_files)

# 创建带样本数的x轴标签
x_labels <- sapply(filtered_types, function(t) {
  n_files <- type_n_files[type == t, n_files]
  if (length(n_files) > 0 && n_files > 0) {
    return(paste0(t, "\n(n=", n_files, ")"))
  } else {
    return(t)
  }
})

# 循环处理每个阈值
for (i in seq_along(thresholds)) {
  threshold <- thresholds[i]
  threshold_name <- threshold_names[i]
  
  cat(sprintf("\n处理阈值: %s\n", threshold_name))
  
  # 过滤掉unknown和MCA类型的数据
  threshold_data <- base_data[get(type_col) %in% filtered_types]
  threshold_data[, mimer_status := get(threshold)]
  # 添加类型列
  threshold_data[, type := get(type_col)]
  threshold_data$type <- factor(threshold_data$type, levels = filtered_types)
  
  # 对于每个类型，只使用包含该类型的文件
  sample_stats_list <- list()
  
  for (current_type in filtered_types) {
    # 获取包含当前类型的所有文件
    files_with_type <- unique(threshold_data[type == current_type, file])
    
    if (length(files_with_type) > 0) {
      # 对于每个包含当前类型的文件，计算MIMER比例
      for (current_file in files_with_type) {
        # 获取该文件当前类型的所有细胞
        type_cells <- threshold_data[file == current_file & type == current_type]
        
        # 计算MIMER细胞数
        mimer_cells <- sum(type_cells$mimer_status == "MIMER", na.rm = TRUE)
        
        # 计算总细胞数
        total_cells <- nrow(type_cells)
        
        # 计算比例
        prop <- ifelse(total_cells > 0, mimer_cells / total_cells, 0)
        
        # 保存结果
        sample_stats_list[[paste(current_type, current_file, sep="_")]] <- data.frame(
          file = current_file,
          type = current_type,
          mimer_cells = mimer_cells,
          total_cells = total_cells,
          prop = prop,
          threshold = threshold_name
        )
      }
    } else {
      cat(sprintf("  警告: 类型 %s 在数据中未找到\n", current_type))
    }
  }
  
  # 合并所有结果
  if (length(sample_stats_list) > 0) {
    sample_stats <- rbindlist(sample_stats_list)
  } else {
    sample_stats <- data.table(
      file = character(),
      type = factor(levels = filtered_types),
      mimer_cells = integer(),
      total_cells = integer(),
      prop = numeric(),
      threshold = character()
    )
  }
  
  # 设置类型为因子
  sample_stats$type <- factor(sample_stats$type, levels = filtered_types)
  
  # 计算每个类型的统计量
  mean_stats <- sample_stats[, 
                             .(mean_prop = mean(prop, na.rm = TRUE),
                               se_prop = sd(prop, na.rm = TRUE) / sqrt(.N),
                               median_prop = median(prop, na.rm = TRUE),
                               sd_prop = sd(prop, na.rm = TRUE),
                               n_samples = .N,
                               n_nonzero = sum(prop > 0)), 
                             by = type
  ]
  
  # 对于没有数据的类型，添加0行
  missing_types <- setdiff(filtered_types, mean_stats$type)
  if (length(missing_types) > 0) {
    for (mt in missing_types) {
      mean_stats <- rbindlist(list(
        mean_stats,
        data.table(
          type = mt,
          mean_prop = 0,
          se_prop = 0,
          median_prop = 0,
          sd_prop = 0,
          n_samples = 0,
          n_nonzero = 0
        )
      ), fill = TRUE)
    }
  }
  
  # 添加零值百分比
  mean_stats[, zero_percent := ifelse(n_samples > 0, 
                                      round((n_samples - n_nonzero) / n_samples * 100, 1), 
                                      0)]
  mean_stats[, threshold := threshold_name]
  
  # 确保类型顺序
  mean_stats$type <- factor(mean_stats$type, levels = filtered_types)
  mean_stats <- mean_stats[order(type)]
  
  # 存储结果
  all_threshold_stats[[threshold_name]] <- list(
    prop_data = sample_stats,
    mean_stats = mean_stats
  )
  
  # 进行Wilcoxon检验 - 所有组与Nor组比较
  # 提取Nor组数据
  nor_data <- sample_stats[type == nor_group, prop]
  n_nor <- length(nor_data)
  cat(sprintf("  Nor组样本数: %d\n", n_nor))
  
  # 为每个非Nor组进行检验
  comparison_results[[threshold_name]] <- list()
  
  for (current_type in filtered_types[filtered_types != nor_group]) {
    # 提取当前组数据
    current_data <- sample_stats[type == current_type, prop]
    n_current <- length(current_data)
    
    cat(sprintf("  %s样本数: %d\n", current_type, n_current))
    
    # 与Nor组比较
    if (n_nor >= 2 && n_current >= 2) {
      # 进行Wilcoxon检验
      test_result <- wilcox.test(current_data, nor_data, exact = FALSE)
      p_value <- test_result$p.value
    } else {
      p_value <- NA
      cat(sprintf("  警告: %s vs Nor比较样本不足，至少需要每组2个样本\n", current_type))
    }
    
    # 存储比较结果
    comparison_results[[threshold_name]][[current_type]] <- list(
      p.value = p_value,
      n_current = n_current,
      n_nor = n_nor
    )
    
    cat(sprintf("  %s vs %s: p = %.4e\n", current_type, nor_group, 
                ifelse(is.na(p_value), NA, p_value)))
  }
  
  # 输出每个类型的样本数
  cat("  每个类型的样本数（文件数）:\n")
  print(mean_stats[, .(type, n_samples)][order(factor(type, levels = filtered_types))])
}

# 合并所有阈值的数据
all_prop_data <- rbindlist(lapply(all_threshold_stats, function(x) x$prop_data))
all_mean_stats <- rbindlist(lapply(all_threshold_stats, function(x) x$mean_stats))

# 将threshold转换为因子，确保正确的顺序
all_mean_stats$threshold <- factor(all_mean_stats$threshold, levels = threshold_names)
if (nrow(all_prop_data) > 0) {
  all_prop_data$threshold <- factor(all_prop_data$threshold, levels = threshold_names)
}

# 确保类型顺序
all_mean_stats$type <- factor(all_mean_stats$type, levels = filtered_types)
if (nrow(all_prop_data) > 0) {
  all_prop_data$type <- factor(all_prop_data$type, levels = filtered_types)
}

# 准备比较结果的标注数据
annotation_data_list <- list()

for (thresh in threshold_names) {
  p_values <- comparison_results[[thresh]]
  
  # 为每个非Nor组创建标注数据
  for (current_type in filtered_types[filtered_types != nor_group]) {
    if (!is.null(p_values[[current_type]]) && !is.na(p_values[[current_type]]$p.value)) {
      # 获取当前阈值和类型下的统计信息
      type_stats <- all_mean_stats[threshold == thresh & type == current_type]
      
      if (nrow(type_stats) > 0) {
        # 获取均值
        mean_prop <- type_stats$mean_prop
        se_prop <- type_stats$se_prop
        
        # 获取p值
        p_value <- p_values[[current_type]]$p.value
        
        # 计算x位置（类型的索引）
        x_pos <- which(filtered_types == current_type)
        
        # 创建标注数据
        annotation_data_list[[paste(thresh, current_type, sep="_")]] <- data.frame(
          threshold = thresh,
          type = current_type,
          x_pos = x_pos,
          mean_prop = mean_prop,
          se_prop = se_prop,
          p_value = p_value
        )
      }
    }
  }
}

# 合并所有标注数据
if (length(annotation_data_list) > 0) {
  annotation_data <- rbindlist(annotation_data_list)
} else {
  annotation_data <- data.frame(
    threshold = character(),
    type = character(),
    x_pos = numeric(),
    mean_prop = numeric(),
    se_prop = numeric(),
    p_value = numeric()
  )
}

# 格式化为保留2位有效数字的p值
format_p_value <- function(p) {
  if (p < 0.0001) {
    # 对于非常小的p值，使用科学计数法
    formatC(p, format = "e", digits = 2)
  } else if (p < 0.001) {
    # 对于小于0.001的p值，显示4位小数
    format(round(p, 4), nsmall = 4, scientific = FALSE)
  } else if (p < 0.01) {
    # 对于小于0.01的p值，显示3位小数
    format(round(p, 3), nsmall = 3, scientific = FALSE)
  } else if (p < 0.1) {
    # 对于小于0.1的p值，显示2位小数
    format(round(p, 2), nsmall = 2, scientific = FALSE)
  } else {
    # 对于较大的p值，显示2位小数
    format(round(p, 2), nsmall = 2, scientific = FALSE)
  }
}

# 添加显著性标注和格式化后的p值标签
if (nrow(annotation_data) > 0) {
  annotation_data[, significance := ifelse(p_value < 0.001, "***",
                                           ifelse(p_value < 0.01, "**",
                                                  ifelse(p_value < 0.05, "*", "ns")))]
  # 创建标签，包含格式化后的p值和显著性标注
  annotation_data[, formatted_p := sapply(p_value, format_p_value)]
  annotation_data[, label := paste0("p=", formatted_p, " ", significance)]
  
  # 将threshold转换为因子
  annotation_data$threshold <- factor(annotation_data$threshold, levels = threshold_names)
  
  # 将type转换为因子
  annotation_data$type <- factor(annotation_data$type, levels = filtered_types)
}

# 设置图表标题
if (analysis_type == "ME") {
  chart_title <- "MIMER+ Spot Proportion Across ME Types"
  x_label <- "ME Type"
} else {
  chart_title <- "MIMER+ Spot Proportion Across Epi/Cancer Types"
  x_label <- "Epi/Cancer Type"
}

# 计算y轴上限
y_max <- max(all_mean_stats$mean_prop + all_mean_stats$se_prop, na.rm = TRUE) * 1.3
# 确保y轴有足够的空间显示标注
if (nrow(annotation_data) > 0) {
  y_max <- max(y_max, max(annotation_data$mean_prop + annotation_data$se_prop) * 1.3)
}

# 创建主图
p_main <- ggplot(all_mean_stats, aes(x = type, y = mean_prop, group = threshold, color = threshold)) +
  # 添加均值点
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  # 添加连接线
  geom_line(linewidth = 1, position = position_dodge(width = 0.3)) +
  # 添加误差棒
  geom_errorbar(
    aes(ymin = mean_prop - se_prop, ymax = mean_prop + se_prop),
    width = 0.2,
    linewidth = 0.8,
    position = position_dodge(width = 0.3)
  )

# 只在有显著结果时添加显著性标记
if (nrow(annotation_data) > 0) {
  # 计算dodge位置
  dodge_width <- 0.3
  n_thresholds <- length(threshold_names)
  
  # 创建偏移位置
  offsets <- seq(from = -dodge_width*(n_thresholds-1)/2, 
                 to = dodge_width*(n_thresholds-1)/2, 
                 length.out = n_thresholds)
  
  # 合并dodge偏移
  annotation_data[, dodge_offset := offsets[as.numeric(factor(threshold, levels = threshold_names))]]
  annotation_data[, x_dodged := as.numeric(type) + dodge_offset]
  
  # 计算y位置（在均值+标准误的上方再加一些空间）
  annotation_data[, y_pos := mean_prop + se_prop + 0.01 * y_max]
  
  p_main <- p_main + 
    geom_text(
      data = annotation_data,
      aes(x = x_dodged, y = y_pos, label = label, color = threshold),
      size = 3.5,
      fontface = "bold",
      show.legend = FALSE
    )
}

# 完成图表设置
p_main <- p_main +
  # 应用颜色方案
  scale_color_manual(values = threshold_colors) +
  # 设置x轴标签，包括样本数
  scale_x_discrete(limits = filtered_types, labels = x_labels) +
  # 设置y轴
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, y_max),
    expand = expansion(mult = c(0.1, 0.1))
  ) +
  # 标签和标题
  labs(
    title = chart_title,
    subtitle = "Mean ± SEM for different thresholds\nWilcoxon test vs Nor group: *p<0.05, **p<0.01, ***p<0.001",
    x = x_label, 
    y = "Proportion of MIMER+ Spots",
    color = "Threshold"
  ) +
  # 主题设置
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    plot.title = element_text(
      hjust = 0.5, 
      face = "bold", 
      size = 16,
      margin = margin(b = 5)
    ),
    plot.subtitle = element_text(
      hjust = 0.5, 
      size = 10,
      color = "gray40",
      margin = margin(b = 10)
    ),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_line(color = "gray95", linewidth = 0.3),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  )

# 显示主图
print(p_main)

# Figure 4g -----
library(dplyr)
library(ggplot2)
library(stringr)
library(tibble)

dir <- dirname(rstudioapi::getSourceEditorContext()$path)
patient_scores <- readRDS(paste0(dir, "/paired_patient_level_scores.rds"))
paired_stats <- readRDS(paste0(dir, "/paired_stats_20.rds"))

sig_stats <- paired_stats[paired_stats$sig_fdr=="TRUE"|paired_stats$pathway=="Tex_CD57"|paired_stats$pathway=="Treg_CD39"|paired_stats$pathway=="Treg_OX40",]
  
plot_df <- patient_scores %>%
  dplyr::inner_join(
    sig_stats %>%
      dplyr::transmute(
        mimer_supertype,
        pathway,
        pathway_label = paste0(pathway, "\nFDR=", scales::pvalue(p_adj, accuracy = 0.001))
      ),
  by = c("mimer_supertype", "pathway")
)

for (this_super in unique(plot_df$mimer_supertype)) {
  df_sub <- plot_df %>%
    dplyr::filter(mimer_supertype == this_super)

  if (length(unique(df_sub$pathway_label)) < 4) {
    df_sub$pathway_label <- factor(df_sub$pathway_label, levels = c(unique(df_sub$pathway_label), "", " ", "  "))
  }
  p <- ggplot(df_sub, aes(x = condition, y = mean_score, color = condition)) +
    geom_line(aes(group = patient_id), color = "grey80", linewidth = 0.5) +
    geom_boxplot(aes(fill = condition), width = 0.35, alpha = 0.18, outlier.shape = NA) +
    geom_point(size = 2.2, alpha = 0.85) +
    facet_wrap(~ pathway_label, nrow = 2, ncol = 3,scales = "free_y",drop = FALSE) +
    scale_color_manual(values = c("other" = "#4E79A7", "MIMER" = "#E15759")) +
    scale_fill_manual(values = c("other" = "#4E79A7", "MIMER" = "#E15759")) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "#F2F2F2", color = "#D9D9D9"),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold")
    ) +
    labs(
      title = paste0("CODEX functional comparison - ", this_super),
      x = NULL,
      y = "Patient-level mean pathway score (within-patient z-scored markers)"
    )
  
  ggplot2::ggsave(
    filename = file.path(dir,paste0("Boxplot_", this_super, "_FDR_sig.pdf")),
    plot = p,
    width = 8,
    height = 6.5
  )
}


# Figure 4m
library(survival)
library(survminer)
pdf(file="hazard_ratio.pdf", width=8, height=6)
data <- read.table("data1.txt",header=T)

model <- coxph(Surv(Survival_time, Survival_status) ~ Age_2 + Grade_2 + AJCC_2 + MIMER, data = data)
ggforest(model, data=data, main='Hazard ratio', refLabel = '1', noDigits =3)

fit <- survfit(Surv(Survival_time, Survival_status) ~MIMER, data=data)
ggsurvplot(fit, data = data, conf.int = TRUE, pval = TRUE, fun = "pct", xlim = c(0, 80), risk.table = TRUE, size = 1, linetype = "strata", palette = c("#E64B35CC", "#4DBBD5CC"), legend = "top",legend.labs = c("MIMER high", "MIMER low"), legend.title = "")


# Figure 4l -----
codex2 = codex
names(codex2@images) = unique(codex2$Sample)
unique_types_cols = rep('grey',length(unique_types))
names(unique_types_cols) = unique_types
Mimer = c("SPP1+Mac",'ISG15+Mac','Ki-67+Mac',
          'CD4T-aTreg','Tcells.cc',
          'PDCD1+Tex','PDCD1+Tex2','LAG3+Tex',
          'COL1A1+CAF','POSTN+CAF',
          'PLVAP+Endo','Vimentin+Endo')
unique_types_cols[names(unique_types_cols) %in% Mimer] = 'red'


unique_types <- unique(Mimer) 

CNS = unique(codex2$neighborhood10_new)[unique(codex2$neighborhood10_new)!='other']

type_dist_results_list <- list()
type_dist_df_CN_list <- list()
type_dist_df_CN_sample_list <- list()
Samples = names(codex2@images) 

for (sa in Samples) {
  location = GetTissueCoordinates(codex2@images[[sa]])
  location %>%
    column_to_rownames("cell") -> pos
  dist_mat <- as.matrix(dist(pos))
  dist_mat <- (dist_mat + t(dist_mat))/2
  diag(dist_mat) <- 0
  colnames(dist_mat) <- rownames(pos)
  rownames(dist_mat) <- rownames(pos)
  distance_matrix_df = dist_mat
  dist_mat[1:4,1:4]
  
  type_combinations <- expand.grid(type1 = unique_types, type2 = unique_types) %>%
    filter(as.numeric(type1) <= as.numeric(type2)) # 确保只处理一次 (A-B 和 B-A 算一次)
  cell_ids = codex2@meta.data %>%
    mutate(cell_id = rownames(.),
           cell_type = subCelltype) %>% 
    filter(Sample==sa) %>%
    dplyr::select(cell_id,cell_type,neighborhood10)
  # 用于存储结果的列表
    cell_id_to_type = cell_ids

    # 遍历所有细胞类型组合
    for (i in 1:nrow(type_combinations)) {
      type1 <- type_combinations$type1[i] %>% as.character()
      type2 <- type_combinations$type2[i] %>% as.character()
      
      
      # 获取属于当前类型的细胞ID
      cells_type1_ids <- cell_id_to_type %>%
        filter(cell_type == type1) %>%
        pull(cell_id)
      
      cells_type2_ids <- cell_id_to_type %>%
        filter(cell_type == type2) %>%
        pull(cell_id)
      
      # 确保有细胞属于这些类型
      if (length(cells_type1_ids) == 0 || length(cells_type2_ids) == 0) {
        next # 跳过没有细胞的类型
      }
      
      # 从距离矩阵中提取所有 Type1 和 Type2 之间两两距离的子集
      # 使用行名和列名进行子集选择
      sub_matrix <- distance_matrix_df[cells_type1_ids, cells_type2_ids, drop = FALSE]
      
      # 将子矩阵展平为向量
      flat_distances <- as.vector(sub_matrix)
      
      # 如果是同一种类型 (Type1 == Type2)，我们需要处理对角线 (自距离为0)
      if (type1 == type2) {
        # 移除值为0的自连接距离
        flat_distances <- flat_distances[flat_distances != 0]
      }
      
      # 计算各种距离统计量
      if (length(flat_distances) > 0) {
        mean_dist <- mean(flat_distances)
        min_dist <- min(flat_distances)
        q25 = quantile(flat_distances,0.25)
        q20 = quantile(flat_distances,0.20)
        max_dist <- max(flat_distances)
      } else {
        mean_dist <- NA
        min_dist <- NA
        max_dist <- NA
        q25 <- NA
        q20 <- NA
      }
      
      # 存储结果
      type_dist_results_list[[i]] <- data.frame(
        Type1 = type1,
        Type2 = type2,
        Mean_Distance = mean_dist,
        Min_Distance = min_dist,
        Max_Distance = max_dist
      )
    }
    
    # 将结果列表合并为数据框
    type_dist_df <- do.call(rbind, type_dist_results_list)
    
    type_dist_df$Sample = sa
  type_dist_df_CN_sample_list[[sa]] = type_dist_df
  print(paste0(i,CN,sa,"done!!!"))
}

type_dist_df_CN_sample_list_df_all  = do.call(rbind,type_dist_df_CN_sample_list)
save(type_dist_df_CN_sample_list_df_all,file = "all_region_distinct.Rdata")
head(type_dist_df_CN_sample_list_df_all)
codex@meta.data %>%
  distinct(Sample,.keep_all = T) %>%
  dplyr::select(Sample,Tissue,Patient_ID) -> meta_tissue
theme_niwot <- function(){
  theme(
    legend.key=element_blank(),   # 图例键为空
    legend.text = element_text(color="black",size=10), # 定义图例文本
    legend.spacing.x=unit(0.1,'cm'), # 定义文本书平距离
    legend.key.width=unit(0.5,'cm'), # 定义图例水平大小
    legend.key.height=unit(0.5,'cm'), # 定义图例垂直大小
    legend.background=element_blank()) 
}
tissue_cols <- c("Adj" = "#377EB8", "Tumor" = "#E41A1C") 
type_dist_df_CN_sample_list_df_all %>%
  group_by(Sample,) %>%
  summarise(mean_distance = mean(Min_Distance)) %>%
  left_join(meta_tissue,by = "Sample") %>%
  # filter(mean_distance <25) %>%
  ggplot(.,aes(x=Tissue,y=mean_distance)) +
  geom_boxplot(aes(fill = Tissue),lwd = 0.3,width=0.5,outliers = F, position = position_dodge(width = 0.2))+
  geom_jitter(aes(color=Tissue),
             width = 0.2)+
  stat_compare_means(comparisons = list(c("Adj","Tumor")),method = 'wilcox.test')+
  scale_fill_manual(values = tissue_cols)+
  scale_color_manual(values = alpha(tissue_cols,0.5))+labs(y='The average minimum distance\n between Mimer members')+
  theme_bw()
ggsave(filename = 'Average_minimum_distance.pdf',width = 3,height = 4)

type_dist_df_CN_sample_list_df_all %>%
  group_by(Sample,) %>%
  summarise(mean_distance = mean(Mean_Distance)) %>%
  left_join(meta_tissue,by = "Sample") %>%
  ggplot(.,aes(x=Tissue,y=mean_distance)) +
  geom_boxplot(aes(fill = Tissue),lwd = 0.3,width=0.5,outliers = F, position = position_dodge(width = 0.2))+
  geom_jitter(aes(color=Tissue),
              width = 0.2)+
  stat_compare_means(comparisons = list(c("Adj","Tumor")),method = 'wilcox.test')+
  scale_fill_manual(values = tissue_cols)+
  scale_color_manual(values = alpha(tissue_cols,0.5))+labs(y='The average mean distance\n between Mimer members')+
  theme_bw()


type_dist_df_CN_sample_list_df_all %>%
  filter(!is.na(Mean_Distance)) %>%
  group_by(Sample) %>%
  summarise(mean_distance = mean(Mean_Distance)) -> mean_min_distance 

mean_min_distance = mean_min_distance %>%
  left_join(meta_tissue,by='Sample') %>%
  filter(Tissue =='Tumor')
survival_info <- readxl::read_xlsx("~/ESCC_codex/CN/data.xlsx",sheet = 2)
survival_info <- survival_info %>%
  filter(Patient_ID %in% mean_min_distance$Patient_ID)%>%
  left_join(mean_min_distance,by='Patient_ID')

cn_var= "mean_distance"
group_var = "mean_distance_gorup"
threshold <- mean(survival_info[[cn_var]], na.rm = TRUE)
survival_info$time = as.numeric(survival_info$`Survival_time (months)`)
survival_info$status = survival_info$`Survival_status (1, dead; 0, alive)`
value <- surv_cutpoint(survival_info, time = "time", event = "status", variables = "mean_distance",minprop = 0.1) 
cut_off <- as.numeric(value[["cutpoint"]][1, 1])

survival_info[[group_var]] <- ifelse(survival_info[[cn_var]] > cut_off, "Far", "Near")
survival_info[[group_var]] <- factor(survival_info[[group_var]], levels = c("Near", "Far"))

library(survival)
# 生存分析
head(survival_info)

fit <- survfit(as.formula(paste("Surv(time, status) ~", group_var)), data = survival_info)  
df <- data.frame(
  time = fit$time,
  surv = fit$surv,
  lower = fit$lower,
  upper = fit$upper,
  strata = rep(names(fit$strata), fit$strata)
)  

cox_formula <- as.formula(paste("Surv(time, status) ~", group_var))
cox_fit <- coxph(cox_formula, data = survival_info)
cox_summary <- summary(cox_fit)

hr <- round(cox_summary$coefficients[1, "exp(coef)"], 3)
p_hr <- signif(cox_summary$coefficients[1, "Pr(>|z|)"], 3)
ci_lower <- round(cox_summary$conf.int[1, "lower .95"], 3)
ci_upper <- round(cox_summary$conf.int[1, "upper .95"], 3)

surv_diff <- survdiff(as.formula(paste("Surv(time, status) ~", group_var)), data = survival_info)
p_logrank <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
p_logrank_label <- signif(p_logrank, 3)

# 合并标签
hr_label <- paste0(
  "n(Far) = 50","\n",
  "n(Near) = 20","\n",
  "Logrank p = ", p_logrank_label, "\n",
  "p(HR) = ", p_hr , "\n",
  "HR (high/low) = ", hr, " (", ci_lower, "-", ci_upper, ")\n"
)
library(survival)
library(survminer)
p <- ggsurvplot(
  fit, 
  data = survival_info, 
  conf.int.style = 'step',
  pval = F, 
  conf.int = TRUE,
  risk.table = TRUE,
  palette = c("red", "#5887ec"),
  legend.title = '',
  risk.table.y.text = TRUE ,
  risk.table.height = 0.25,
  risk.table.title = "", 
  title = paste0("KM Curve for ", cn_var))
p
max_time <- max(survival_info$time, na.rm = TRUE)
p$plot <- p$plot + annotate("text", x = max_time, y = 0.8, label = hr_label, size = 4, hjust=1) + 
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
p

# -------------Fig 4k-----------------
#######################################
library(data.table)
library(tidyverse)
library(ggplot2)
library(patchwork)

tls_analysis <- fread('C:/Users/Ron/Desktop/Figs/F4/mimer_pos_10_TLS_mark.csv')
annotation_df <- fread('C:/Users/Ron/Desktop/Figs/F4/TLS_G1G2G3_heatmap_annotation_df_add_geo_mean_mimer.csv')
annotation_df_mark_mimer_pos_G3 <- left_join(annotation_df,tls_analysis[,c('TLS_ID','is_mimer_pos_tls')],by='TLS_ID')

# 获取需要绘图的变量
all_cols <- names(annotation_df_mark_mimer_pos_G3)
start_idx <- which(all_cols == "cluster_size")
end_idx <- which(all_cols == "peripheral_Tumorsuppression_avg")
plot_vars <- all_cols[start_idx:end_idx]

# 筛选需要比较的两组数据
plot_data <- annotation_df_mark_mimer_pos_G3 %>%
  select(TLS_ID, is_mimer_pos_tls, all_of(plot_vars)) %>%
  filter(is_mimer_pos_tls %in% c("other_G3", "mimer_G3"))

# 创建函数来生成单个变量的图形
create_single_plot <- function(var_name, data) {
  # 提取当前变量的数据
  var_data <- data %>%
    select(TLS_ID, is_mimer_pos_tls, value = all_of(var_name)) %>%
    filter(!is.na(value))
  
  # 计算统计检验
  p_value <- NA
  sig_label <- "NA"
  
  if (length(unique(var_data$is_mimer_pos_tls)) == 2) {
    test_result <- tryCatch(
      wilcox.test(value ~ is_mimer_pos_tls, data = var_data),
      error = function(e) NULL
    )
    
    if (!is.null(test_result)) {
      p_value <- test_result$p.value
      sig_label <- case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    }
  }
  
  # 计算y轴范围
  y_range <- range(var_data$value, na.rm = TRUE)
  y_max <- y_range[2]
  y_min <- y_range[1]
  y_span <- y_max - y_min
  
  # 调整y轴上限，为标注留出空间
  y_limit_upper <- y_max + 0.2 * y_span
  
  # 创建当前变量的箱线图
  p <- ggplot(var_data, aes(x = is_mimer_pos_tls, y = value, fill = is_mimer_pos_tls)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.2, show.legend = FALSE) +
    geom_text(
      x = 1.5,
      y = y_max + 0.1 * y_span,
      label = sig_label,
      size = 5,
      # fontface = "bold",
      inherit.aes = FALSE
    ) +
    labs(
      title = var_name,
      y = "Value",
      x = ""
    ) +
    scale_fill_brewer(palette = "Set1") +
    coord_cartesian(ylim = c(y_min, y_limit_upper)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
    )
  
  return(p)
}

# 为每个变量创建图形并存储在列表中
plot_list <- lapply(plot_vars, function(var) {
  create_single_plot(var, plot_data)
})

# 使用patchwork拼图
n_cols <- 6
combined_plot <- wrap_plots(plot_list, ncol = n_cols) +
  plot_annotation(
    title = "Comparison of TLS Features: other_G3 vs mimer_G3",
    caption = "*** p<0.001, ** p<0.01, * p<0.05, ns not significant"
  )

# 设置图形尺寸并显示
options(repr.plot.width = 18, repr.plot.height = 3 * ceiling(length(plot_vars)/n_cols))
print(combined_plot)
