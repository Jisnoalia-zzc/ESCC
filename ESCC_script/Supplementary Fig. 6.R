#### S6A
library(ggplot2)
library(ggpubr)
library(patchwork)
library(data.table)

setwd('C:/Users/Ron/Desktop/Figs/S6/')
plot_data <- fread("C:/Users/Ron/Desktop/Figs/S6/S6_a_GSVA_plot_data.csv")
gene_counts <- fread("C:/Users/Ron/Desktop/Figs/S6/S6_a_cluster_gene_counts.csv")

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

combined_plot <- wrap_plots(plot_list, ncol = 4) + 
  plot_annotation(title = "GSVA score(Responder vs Non-Responder)",
                  subtitle = paste("DESeq2 filtered genes (p<0.2)"),
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
                                plot.subtitle = element_text(hjust = 0.5, size = 12)))

legend_plot <- ggplot(plot_data, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#00BA38", "#F8766D"),
                    name = "Response Group",
                    labels = c("Responder (CR+PR)", "Non-Responder (SD+PD)")) +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"))

legend <- get_legend(legend_plot)

final_plot <- combined_plot / legend + 
  plot_layout(heights = c(10, 1))

options(repr.plot.width = 14, repr.plot.height = 6)
print(final_plot)


# fig S6B
library(data.table)
library(ggplot2)
library(pROC)

setwd('C:/Users/Ron/Desktop/Figs/S6/')
roc_curve_data <- fread("S6_b_ROC_curve_data.csv")
auc_results <- fread("S6_b_ROC_AUC_results.csv")
cluster_gene_counts <-fread('S6_b_cluster_gene_counts_ROC.csv')

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

combined_roc_plot <- wrap_plots(plot_list, ncol = 4) + 
  plot_annotation(
    title = "ROC Curves of TLS Cluster Signatures",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

options(repr.plot.width = 16, repr.plot.height = 4)
print(combined_roc_plot)


#### S6C-D
library(survival)
library(survminer)
library(patchwork)
library(data.table)
library(ggplot2)

cat("读取保存的数据...\n")
survival_data <- fread("S6_cd_survival_plot_data.csv")
cluster_summary <- fread("S6_cd_cluster_summary.csv")
survfit_data <- fread("S6_cd_survfit_data.csv")

surv_type <- unique(cluster_summary$Survival_Type)
group_method <- unique(cluster_summary$Group_Method)
clusters <- unique(survival_data$Cluster)

cat(paste0("生存类型: ", surv_type, "\n"))
cat(paste0("分组方法: ", group_method, "\n"))
cat(paste0("Clusters: ", paste(clusters, collapse = ", "), "\n"))

recreate_plot_list <- list()
recreate_table_list <- list()

for (cluster in clusters) {
  cat(paste0("处理cluster: ", cluster, "...\n"))
  
  cluster_data <- survival_data[Cluster == cluster, ]
  cluster_info <- cluster_summary[Cluster == cluster, ]
  
  plot_data <- data.frame(
    ID = cluster_data$ID,
    GSVA_Score = cluster_data$GSVA_Score,
    Status = cluster_data$Survival_Status,
    Time = cluster_data$Survival_Time,
    GSVA_Group = cluster_data$GSVA_Group
  )
  
  xlab_text <- ifelse(surv_type == "OS", "OS Days", "PFS Days")
  
  group_title <- paste0(cluster, "\n(", group_method, " split, n=", cluster_info$Gene_Count, ")")
  
  survfit_data_cluster <- survfit_data[Cluster == cluster, ]
  
  fit <- survfit(
    Surv(Time, Status) ~ GSVA_Group,
    data = plot_data
  )
  
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

cat("拼接生存曲线...\n")
surv_combined <- wrap_plots(recreate_plot_list, ncol = 4) + 
  plot_annotation(
    title = paste0("TLS GeneSets ", surv_type, " Survival Analysis (", group_method, " split)"),
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

legend <- get_legend(
  recreate_plot_list[[1]] + 
    theme(legend.position = "bottom") +
    guides(color = guide_legend(nrow = 1))
)

final_surv_plot <- surv_combined / legend + 
  plot_layout(heights = c(10, 1))

cat("显示生存曲线...\n")
options(repr.plot.width = 16, repr.plot.height = 4)
print(final_surv_plot)

cat("拼接风险表...\n")
table_combined <- wrap_plots(recreate_table_list, ncol = 4) + 
  plot_annotation(
    title = paste0("TLS GeneSets ", surv_type, " Risk Table (", group_method, " split)"),
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )
print(table_combined)


#--------------------------------
#------Sup-Fig.6e
#------------------------------
#----------------------------------------------------------------
library(GSVA)
library(limma)
library(BiocParallel)
library(ggpubr)
library(patchwork)
library(data.table)
library(pROC)
library(DESeq2)

# 1.---------------------------------------------------
combined_data <- fread('./ImmTherapy_cohort/Processed_immcohort_therapy_response_survival_gene_expr.csv')

# 
combined_data$ResponseGroup <- ifelse(
  combined_data$Best_overall_response %in% c("CR", "PR"), "Responder",
  "Non-Responder"
)
combined_data <- combined_data[!is.na(ResponseGroup), ]
combined_data$ResponseGroup <- factor(
  combined_data$ResponseGroup,
  levels = c("Responder", "Non-Responder")
)

gene_start_col <- which(colnames(combined_data) == "TSPAN6")
gene_end_col <- which(colnames(combined_data) == "ResponseGroup") - 1
expr_matrix <- as.matrix(combined_data[, gene_start_col:gene_end_col, with = FALSE])
mode(expr_matrix) <- "integer"
rownames(expr_matrix) <- combined_data$ID

group <- combined_data$ResponseGroup

stopifnot(
  is.numeric(expr_matrix[1, 1]),
  nrow(expr_matrix) == length(group),
  all(!is.na(group))
)
message(paste("有效样本数:", nrow(expr_matrix), "| 基因数:", ncol(expr_matrix)))

# 32. GSVA ---------------------------------------------------------------
maker_list <- list(
  TLS_9 = c("PTGDS", "RBP5", "EIF1AY", "CETP", "SKAP1", "LAT", "CCR6", "CD1D", "CD79B"),
  
  TLS_12 = c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11", "CXCL13"),
  
  TLS_50 = c("FDCSP","CR2","CXCL13","LTF","CD52","MS4A1","CCL19","LINC00926","LTB","CORO1A",
             "CD79B","TXNIP","CD19","LIMD2","CD37","ARHGAP45","BLK","TMC8","CCL21","PTPN6","ATP2A3",
             "IGHM","SPIB","TMSB4X","CXCR4","NCF1","CD79A","ARHGAP9","DEF6","EVL","TBC1D10C","RASAL3",
             "INPP5D","RNASET2","RASGRP2","TNFRSF13C","RAC2","CD22","ARHGEF1","AC103591.3","TRAF3IP3",
             "HLA-DQB1","CD53","ARHGAP4","TRBC2","POU2AF1","TRAF5","OGA","FCRL3","HLA-DQA1")
)

# clusters <- c('Non-TLS', 'peri_TLS_G3', 'TLS_G3', 'peri_TLS_G2', 'TLS_G2', 'peri_TLS_G1', 'TLS_G1')
clusters <- c('TLS_9', 'TLS_12',  'TLS_50')

# Deseq2 diff gene
dds <- readRDS(file = './Genelist/Immcohort_deseq2_dds.RDS')
res <- results(dds, alpha = 0.1) 
# diff_genes <- rownames(subset(res, padj < 0.1 & abs(log2FoldChange) > 0.5))
diff_genes <- rownames(subset(res, pvalue < 0.2))

# 
plot_list <- list()
auc_results <- data.frame()

for (cluster in clusters) {
  valid_genes <- unique(maker_list[[cluster]])
  
  if (length(valid_genes) < 3) {
    message(paste0("跳过 ", cluster, "：有效基因数不足 (", length(valid_genes), ")"))
    next
  }
  
  # 
  gsva_param <- gsvaParam(
    exprData = t(expr_matrix),
    geneSets = list(cluster = valid_genes),
    kcdf = "Poisson" )
  gsva_scores <- gsva(gsva_param, verbose = FALSE, BPPARAM = SerialParam())
  scores <- gsva_scores[1, ]
  
  plot_data <- data.frame(
    Sample = names(scores),
    Score = scores,
    Group = group,
    Cluster = cluster
  )
  
  p <- ggplot(plot_data, aes(x = Group, y = Score, fill = Group)) +
    geom_boxplot(width = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
    scale_fill_manual(values = c("#00BA38", "#F8766D")) +
    labs(title = paste0(cluster, "\n(", length(valid_genes), " genes)"),
         y = "GSVA Score",
         x = NULL) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold",size=10),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1,size=10)) +
    stat_compare_means(
      method = "wilcox.test", 
      label = "p.format",
      tip.length = 0.02,
      size = 5
    )
  
  plot_list[[cluster]] <- p
  
  binary_scores <- plot_data$Score
  binary_labels <- ifelse(plot_data$Group == "Responder", 1, 0)
  
  roc_obj <- roc(response = binary_labels, predictor = binary_scores, quiet = TRUE)
  auc_value <- auc(roc_obj)
  
  auc_results <- rbind(auc_results, data.frame(
    Cluster = cluster,
    AUC = round(auc_value, 3),
    Genes = length(valid_genes)
  ))
}

# 4. -------------------------------------------------------
if (length(plot_list) > 0) {
  combined_plot <- wrap_plots(plot_list, ncol = 4) + 
    plot_annotation(title = "GSVA score(Responder vs Non-Responder)",
                    subtitle = paste(length(valid_genes), "genes"),
                    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
                                  plot.subtitle = element_text(hjust = 0.5, size = 12)))
  
 legend_plot <- ggplot(plot_data, aes(x = Group, y = Score, fill = Group)) +
    geom_boxplot() +
    scale_fill_manual(values = c("#00BA38", "#F8766D"),
                      name = "Response Group",
                      labels = c("Responder (CR+PR)", "Non-Responder (SD+PD)")) +
    theme(legend.position = "bottom",
          legend.title = element_text(face = "bold"))
  
  legend <- get_legend(legend_plot)
  
  final_plot <- combined_plot / legend + 
    plot_layout(heights = c(10, 1))
  
  options(repr.plot.width = 14, repr.plot.height = 6)
  print(final_plot)
  
  #   ggsave("GSVA_Analysis_Plot_TwoGroups.pdf", final_plot, width = 14, height = 16)
  #   fwrite(auc_results, "GSVA_AUC_Results_TwoGroups.csv")
} else {
  warning("没有足够的有效基因集生成图表，请检查差异基因筛选和基因集匹配")
}

options(repr.plot.width = 5, repr.plot.height = 5)

ggplot(auc_results, aes(x='AUC',y=Cluster, fill=AUC)) +
  geom_tile() + 
  geom_text(
    aes(label = sprintf("%.2f", AUC)),
    color = "black", 
    size = 8
  ) +
  scale_fill_gradient2(low="white", high="red", midpoint=0.5)+
  theme(axis.text = element_text(size=12,),
        axis.text.x = element_text(angle = 45,hjust = 1)
  )

ggsave(plot = final_plot,filename = './Genelist/TLS_9_12_50_immcohort_responder_boxplot.pdf',width = 14,height = 6)


#--------------------------------
#------Sup-Fig6.f 
#------------------------------
#----------------------------------------------------------------
library(GSVA)
library(limma)
library(BiocParallel)
library(ggpubr)
library(patchwork)
library(data.table)
library(pROC)
library(DESeq2)
library(ggplot2)

# 1.---------------------------------------------------
combined_data <- fread('./ImmTherapy_cohort/Processed_immcohort_therapy_response_survival_gene_expr.csv')

# 
combined_data$ResponseGroup <- ifelse(
  combined_data$Best_overall_response %in% c("CR", "PR"), "Responder",
  "Non-Responder"
)
combined_data <- combined_data[!is.na(ResponseGroup), ]
combined_data$ResponseGroup <- factor(
  combined_data$ResponseGroup,
  levels = c("Responder", "Non-Responder")
)

gene_start_col <- which(colnames(combined_data) == "TSPAN6")
gene_end_col <- which(colnames(combined_data) == "ResponseGroup") - 1
expr_matrix <- as.matrix(combined_data[, gene_start_col:gene_end_col, with = FALSE])
mode(expr_matrix) <- "integer"
rownames(expr_matrix) <- combined_data$ID

group <- combined_data$ResponseGroup

stopifnot(
  is.numeric(expr_matrix[1, 1]),
  nrow(expr_matrix) == length(group),
  all(!is.na(group))
)
message(paste("有效样本数:", nrow(expr_matrix), "| 基因数:", ncol(expr_matrix)))

# 2. DESeq2---------------------------------------
# colData <- data.frame(
#   row.names = rownames(expr_matrix),
#   condition = group
# )

# dds <- DESeqDataSetFromMatrix(
#   countData = t(expr_matrix),
#   colData = colData,
#   design = ~ condition
# )

# dds <- DESeq(dds)
dds <- readRDS(file = './Genelist/Immcohort_deseq2_dds.RDS')
res <- results(dds, alpha = 0.1) 
diff_genes <- rownames(subset(res, pvalue < 0.2))

# 3. --------------------------------------------------

maker_list <- list(
  TLS_9 = c("PTGDS", "RBP5", "EIF1AY", "CETP", "SKAP1", "LAT", "CCR6", "CD1D", "CD79B"),
  
  TLS_12 = c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11", "CXCL13"),
  
  TLS_50 = c("FDCSP","CR2","CXCL13","LTF","CD52","MS4A1","CCL19","LINC00926","LTB","CORO1A",
             "CD79B","TXNIP","CD19","LIMD2","CD37","ARHGAP45","BLK","TMC8","CCL21","PTPN6","ATP2A3",
             "IGHM","SPIB","TMSB4X","CXCR4","NCF1","CD79A","ARHGAP9","DEF6","EVL","TBC1D10C","RASAL3",
             "INPP5D","RNASET2","RASGRP2","TNFRSF13C","RAC2","CD22","ARHGEF1","AC103591.3","TRAF3IP3",
             "HLA-DQB1","CD53","ARHGAP4","TRBC2","POU2AF1","TRAF5","OGA","FCRL3","HLA-DQA1")
)

# clusters <- c('Non-TLS', 'peri_TLS_G3', 'TLS_G3', 'peri_TLS_G2', 'TLS_G2', 'peri_TLS_G1', 'TLS_G1')
clusters <- c('TLS_9',  'TLS_12', 'TLS_50')

plot_list <- list()
auc_results <- data.frame()

for (cluster in clusters) {
  valid_genes <- unique(maker_list[[cluster]])
  
  if (length(valid_genes) < 3) {
    message(paste0("跳过 ", cluster, "：有效基因数不足 (", length(valid_genes), ")"))
    next
  }
  
  gsva_param <- gsvaParam(
    exprData = t(expr_matrix),  # 仅差异基因
    geneSets = list(cluster = valid_genes),
    kcdf = "Poisson"  # Count数据用Poisson
  )
  gsva_scores <- gsva(gsva_param, verbose = FALSE, BPPARAM = SerialParam())
  scores <- gsva_scores[1, ]
  
  binary_labels <- ifelse(group == "Responder", 1, 0)
  roc_obj <- roc(response = binary_labels, predictor = scores, quiet = TRUE)
  auc_value <- auc(roc_obj)
  
  roc_data <- data.frame(
    Specificity = roc_obj$specificities,
    Sensitivity = roc_obj$sensitivities
  )
  
  p <- ggplot(roc_data, aes(x = 1 - Specificity, y = Sensitivity)) +
    geom_line(color = "steelblue", size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    labs(title = paste0(cluster, "\n(n = ", length(valid_genes), ")"),
         x = "1 - Specificity (False Positive Rate)",
         y = "Sensitivity (True Positive Rate)") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      panel.grid.minor = element_blank(),
      aspect.ratio = 1
    ) +  annotate("text", 
             x = 0.75, y = 0.25, 
             label = paste0("AUC = ", round(auc_value, 3)),
             size = 5, 
             color = "red", 
             fontface = "bold")
  
  plot_list[[cluster]] <- p
  
  auc_results <- rbind(auc_results, data.frame(
    Cluster = cluster,
    AUC = round(auc_value, 3),
    Genes = length(valid_genes)
  ))
}

# 4. 
combined_roc_plot <- wrap_plots(plot_list, ncol = 4) + 
  plot_annotation(
    title = "ROC Curves of TLS Cluster Signatures",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

options(repr.plot.width = 16, repr.plot.height = 4)
pdf(file = "./Genelist/TLS_9_12_50_ROC_Curves.pdf", width = 16, height = 4)
print(combined_roc_plot)
dev.off()

# 5. 
# fwrite(auc_results, "./Genelist/core_loop10_AUC_Results.csv")


#### S6g 6h 6i
library(data.table)
library(ggplot2)
library(pROC)

setwd('C:/Users/Ron/Desktop/Figs/S6')
# 
boxplot_data_reload <- fread( "S6_ghi_boxplot_data.csv")
roc_data_reload <- fread("S6_ghi_roc_data.csv")
auc_data_reload <- fread("S6_ghi_auc_data.csv")
surv_data_reload <- fread("S6_ghi_survival_data.csv")

# 
if (nrow(boxplot_data_reload) > 0) {
  cat("\n=== 从保存的数据重新绘制箱线图 ===\n")
  
  # 
  clusters_reload <- unique(boxplot_data_reload$Cluster)
  boxplot_list_reload <- list()
  
  for (cluster in clusters_reload) {
    cluster_data <- boxplot_data_reload[boxplot_data_reload$Cluster == cluster, ]
    gene_count <- unique(cluster_data$Genes)[1]
    
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
    cluster_roc_data <- roc_data_reload[roc_data_reload$Cluster == cluster, ]
    cluster_auc <- auc_data_reload[auc_data_reload$Cluster == cluster, ]$AUC[1]
    gene_count <- auc_data_reload[auc_data_reload$Cluster == cluster, ]$Genes[1]
    
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
  
  roc_combined_reload <- wrap_plots(rocplot_list_reload, ncol = 4) +
    plot_annotation(
      title = "ROC Curves for Response Prediction",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
                    plot.subtitle = element_text(hjust = 0.5, size = 12))
    )
  print(roc_combined_reload)
  
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
    cluster_surv_data <- surv_data_reload[surv_data_reload$Cluster == cluster, ]
    
    median_score <- median(cluster_surv_data$GSVA_Score)
    cluster_surv_data$GSVA_Group <- ifelse(cluster_surv_data$GSVA_Score >= median_score,
                                           "High GSVA", "Low GSVA")
    group_title <- paste0(cluster, "\n(Median Split)")
    
    fit <- survfit(Surv(Time, Status) ~ GSVA_Group, data = cluster_surv_data)
    
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
    
    # 保存重新绘制的图形（可选）
    # ggsave(paste0(outpath, "S6_ghi_boxplot_reloaded.pdf"), box_combined_reload, width = 14, height = 6)
    # ggsave(paste0(outpath, "S6_ghi_roc_reloaded.pdf"), roc_combined_reload, width = 16, height = 4)
    # ggsave(paste0(outpath, "S6_ghi_survival_reloaded_", surv_type, "_", group_method, ".pdf"),
    #        final_surv_plot_reload, width = 16, height = 4)
  }
}

cat("\n=== 数据重新绘制完成 ===")
