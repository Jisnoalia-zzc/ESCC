#Extended Data Fig4 -----
#Extended Data Fig4a -----
fc <- read_csv("~/ESCC_codex/CN/codex_CN_subCelltype_k_16_windows_20_without_epi.csv")
neighborhood_mat <- as.matrix(fc[2:ncol(fc)])
rownames(neighborhood_mat) <- paste0("CN",rownames(fc))
palette_length <- 100
my_color <- colorRampPalette(c("#2696f2","#58a3e8", "white","#ff8f6b","#d9100b"))(palette_length)
my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05,3, length.out=floor(palette_length/2)))

pdf('~/ESCC_codex/3-codex_without_epi_CN_k16.pdf',width = 12,height = 6)
pheatmap::pheatmap(neighborhood_mat,scale = 'none',clustering_method='average',
                   color = my_color,
                   # display_numbers = T,
                   breaks = my_breaks,
                   cellwidth = 10,
                   cellheight = 10,
                   fontsize_row = 9,
                   fontsize_col = 9,
                   clustering_distance_rows = "euclidean",
                   clustering_distance_cols = "euclidean")
dev.off()

#Extended Data Fig4b -----
load( "/BGFS1/home/zhangzc/ESCC_codex/cor/codex_avg.Rdata")

load("/BGFS1/home/zhangzc/ESCC_codex/cor/sp_avg.Rdata")
p2r = openxlsx::read.xlsx('/BGFS1/home/zhangzc/ESCC_codex/cor/rna_protein_df_edit.xlsx')
codex <- NormalizeData(object = codex, normalization.method = "CLR", margin = 2)
codex  = ScaleData(codex)
codex3 = subset(codex,neighborhood10_new=='other',invert=T)
Idents(codex3) = 'neighborhood10_new'
avg = AverageExpression(codex3,layer='data')
avg = avg$Akoya
codex_avg = avg 
codex_avg = as.data.frame(codex_avg)
head(p2r)
table(p2r$Protein %in% rownames(codex_avg))
p2r$Protein[p2r$Protein %in% rownames(codex_avg)]
setdiff(p2r$Protein,rownames(codex_avg))
setdiff(rownames(codex_avg),p2r$Protein)
p2r$Protein = gsub("Pan.Cytokeratin","Pan-Cytokeratin",
                   gsub("HistoneH3.pSer28","HistoneH3-pSer28",
                        gsub("Keratin8.18","Keratin8-18",
                             gsub("b-Catenin","b-Catenin1",
                                  gsub('COL.1',"COL-1",
                                       gsub("PLVAP.PV.1","PLVAP-PV-1",
                                            gsub("Keratin8.18","Keratin8-18",
                                                 gsub("E.cadherin","E-cadherin",
                                                      gsub("Bcl.2","Bcl-2",p2r$Protein)))))))))
sp_avg
sp_avg %>%
  mutate(RNA = rownames(.)) %>%
  left_join(p2r,by = "RNA") %>%
  filter(RNA !='KRT18') %>%
  column_to_rownames("Protein") %>%
  dplyr::select(-RNA) -> sp_avg2


all(rownames(sp_avg2)==rownames(codex_avg))
setdiff(rownames(sp_avg2),rownames(codex_avg))
setdiff(rownames(codex_avg),rownames(sp_avg2))

sp_avg2 = sp_avg2[match(rownames(codex_avg),rownames(sp_avg2)),]
sp_avg3 = t(scale(t(sp_avg2)))
codex_avg2 = t(scale(t(codex_avg)))
# codex_avg2 = codex_avg
num_cols_sp <- ncol(sp_avg3)
num_cols_codex <- ncol(codex_avg2)
cross_correlation_matrix <- matrix(NA, nrow = num_cols_sp, ncol = num_cols_codex)
rownames(cross_correlation_matrix) <- colnames(sp_avg3)
colnames(cross_correlation_matrix) <- colnames(codex_avg2)

for (i in 1:num_cols_sp) {
  for (j in 1:num_cols_codex) {
    # 计算 matrix_A 的第 i 列和 matrix_B 的第 j 列之间的相关性
    cross_correlation_matrix[i, j] <- cor(sp_avg3[, i], codex_avg2[, j],method = 'pearson')
  }
}
cross_correlation_matrix2 = as.data.frame(cross_correlation_matrix)
cross_correlation_matrix2 = cross_correlation_matrix2[!rownames(cross_correlation_matrix2) %in% c("CC3","CC5",'CC9',"CC10","CC11","CC12"),]

palette_length <- 100
my_color <- colorRampPalette(c("#2696f2","#58a3e8", "white","#ff8f6b","#d9100b"))(palette_length)
my_breaks <- c(seq(-0.5, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05,0.8, length.out=floor(palette_length/2)))

pdf('sp_codex_cor_pheatmap.pdf',width = 8,height = 10 )
pheatmap::pheatmap(as.matrix(cross_correlation_matrix2),scale = 'none',
                   clustering_distance_rows = "euclidean",
                   clustering_distance_cols = "euclidean",
                   color = my_color,
                   # annotation_col = anno_cols_codex,
                   # annotation_row = anno_cols_scRNA,
                   # annotation_colors = anno_cor_colors,
                   breaks = my_breaks,
                   fontsize_row = 6,fontsize_col = 6,
                   cellwidth  = 10,cellheight  = 10)
dev.off()

library(reshape2)
cross_correlation_matrix %>% melt() %>% group_by(Var2) %>%
  top_n(1,value) ->tmp



num_cols_scRNA <- ncol(sp_avg2)
num_cols_codex <- ncol(codex_avg)
cross_correlation_matrix <- matrix(NA, nrow = num_cols_scRNA, ncol = num_cols_codex)
rownames(cross_correlation_matrix) <- colnames(sp_avg2)
colnames(cross_correlation_matrix) <- colnames(codex_avg)


for (i in 1:num_cols_scRNA) {
  for (j in 1:num_cols_codex) {
    # 计算 matrix_A 的第 i 列和 matrix_B 的第 j 列之间的相关性
    cross_correlation_matrix[i, j] <- cor(sp_avg2[, i], codex_avg[, j],method = 'kendall')
  }
}



pdf('scRNA_codex_cor_pheatmap_data.pdf',width = 8,height = 10 )
pheatmap::pheatmap(cross_correlation_matrix,scale = 'none',cluster_cols = T,width = 4,height = 4)
dev.off()
#Extended Data Fig4c -----
library(ggplot2)
library(dplyr)

theme_niwot <- function(){
  theme(
    legend.key=element_blank(),   # 图例键为空
    legend.text = element_text(color="black",size=8), # 定义图例文本
    legend.spacing.x=unit(0.1,'cm'), # 定义文本书平距离
    legend.key.width=unit(0.5,'cm'), # 定义图例水平大小
    legend.key.height=unit(0.5,'cm'), # 定义图例垂直大小
    legend.background=element_blank()) 
}
neighborhoods$neighborhood10 <- gsub("^CM", "CN", neighborhoods$neighborhood10)

data <- read.csv("/home/users/zhangzhichao/workspace/esca/codex_CN/cn_meta_info_subCelltype_k_16_without_epi.csv", row.names = 1)
cell_color <- setNames(cell_color$color, cell_color$ct)

ratio_subtype_sample <- neighborhoods %>%
  group_by(Sample, subCelltype) %>%
  summarise(count = n()) %>%
  mutate(Ratio = count / sum(count)) %>%
  left_join(neighborhoods %>% select(Sample, Tissue) %>% distinct(), by = "Sample")


ratio_subtype_tissue <- neighborhoods %>%
  group_by(Tissue, subCelltype) %>%
  summarise(count = n()) %>%
  mutate(Ratio = count / sum(count)) %>%
  left_join(neighborhoods %>% select(Sample, Tissue) %>% distinct(), by = "Sample")

ratio_neighbor_sample <- neighborhoods %>%
  group_by(Sample, neighborhood10) %>%
  summarise(count = n()) %>%
  mutate(Ratio = count / sum(count)) %>%
  left_join(neighborhoods %>% select(Sample, Tissue) %>% distinct(), by = "Sample")

ratio_neighbor_tissue <- neighborhoods %>%
  group_by(Tissue, neighborhood10) %>%
  summarise(count = n()) %>%
  mutate(Ratio = count / sum(count)) 


p <- ggplot(ratio_neighbor_tissue ,aes(x=Tissue,y = Ratio,fill = neighborhood10))+
  geom_bar(stat = "identity", position = "fill",  width = 0.9
  ) +
  scale_fill_manual(values = niche_cols)+
  #scale_y_continuous(labels = scales::percent_format()) +
  #scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0))+
  #coord_polar()+
  #theme_minimal()+
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),
    
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size=10),
    #plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.direction = 'vertical',
    #panel.grid.major=element_blank(), 
    #panel.grid.minor=element_blank(), 
    panel.background = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  guides(fill = guide_legend(ncol = 2))+
  theme_niwot()
library(cowplot)

legend <- cowplot::get_legend(p + theme(legend.position = "right"))
plot_no_legend <- p + theme(legend.position = "none")

p <- plot_grid(plot_no_legend, legend, rel_widths = c(1, 1))  # 3:1比例

ggsave(
  filename = "Tissue_neighbor.pdf",  # 支持 PDF/PNG/TIFF 等格式
  plot = p,
  device = "pdf",            # 保存为 PDF
  width = 5,               # A4 宽度 (mm)
  height = 5,              # A4 高度 (mm)
  #units = "mm",              # 单位设为毫米
  dpi = 300                  # 分辨率（默认 300 DPI）
)

ggplot(ratio_subtype_sample, aes(x = Tissue, y = Ratio, fill = Tissue)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(0.8)) +
  geom_jitter(aes(color = "black"), width = 0.2, size = 0.1, alpha = 0.8) +
  facet_wrap(~ subCelltype, scales = "free_y") +   
  scale_fill_manual(values = tissue_cols) +
  theme_void() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) 



mat <- ratio_neighbor_sample %>%  
  select(-count)%>%
  select(-Tissue)%>%
  pivot_wider(
    names_from = neighborhood10,
    values_from = Ratio,
    values_fill = 0
  ) %>% column_to_rownames("Sample")
hc <- hclust(vegan::vegdist(mat, method = 'bray'), method = 'average')
sample_order <- hc$labels[hc$order]
names(niche_cols) <- gsub("^CM", "CN", names(niche_cols))


ratio_neighbor_sample$Sample <- factor(ratio_neighbor_sample$Sample, levels = sample_order)
ratio_neighbor_sample$neighborhood10 <- factor(ratio_neighbor_sample$neighborhood10, levels = paste0("CN", 1:16))
p <- ggplot(ratio_neighbor_sample,aes(x=Sample, y = Ratio, fill = neighborhood10))+
  geom_bar(stat = "identity", position = "fill", width = 0.9
  ) +
  scale_fill_manual(values = niche_cols)+
  
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0))+
  theme_minimal()+
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=8),
    legend.position = "right",
    legend.direction = 'vertical',
    panel.grid.major=element_blank(), 
    panel.grid.minor=element_blank(), 
    panel.background = element_blank(),
    plot.margin = margin(20, 10, 10, 10)
  ) +  theme_niwot()

tissue_bar <- ratio_neighbor_sample %>%
  distinct(Sample, Tissue) %>%
  mutate(group = "Tissue")
tissue_bar$Sample <- factor(tissue_bar$Sample, levels = sample_order)

tissue_cols <- c("Adj" = "#377EB8", "Tumor" = "#E41A1C") 
tissue_plot <- ggplot(tissue_bar, aes(x = Sample, y = group, fill = Tissue)) +
  geom_tile(width = 1, height = 0.1, position = position_dodge(width = 1)) +
  scale_fill_manual(values = tissue_cols) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme(axis.title = element_blank(),
        # axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background = element_blank()) +
  theme_niwot()

final_plot <- p %>% insert_bottom(tissue_plot, height = 0.03)


ggsave(
  filename = "Sample_neighborhood10_cluster_label.pdf",  # 支持 PDF/PNG/TIFF 等格式
  plot = final_plot,
  device = "pdf",            # 保存为 PDF
  width = 20,               # A4 宽度 (mm)
  height = 5,              # A4 高度 (mm)
  #units = "mm",              # 单位设为毫米
  dpi = 300                  # 分辨率（默认 300 DPI）
)
write.csv(neighborhoods, "neighborhoods.csv", row.names = FALSE)

#Extended Data Fig4d -----
ratio_neighbor_sample$Tissue <- factor(ratio_neighbor_sample$Tissue, levels = c("Tumor", "Adj"))
my_comparisons <- list(c("Tumor", "Adj"))
p <- ggplot(ratio_neighbor_sample, aes(x = neighborhood10, y = Ratio, fill = Tissue)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(0.8)) +
  geom_jitter(color = 'black',  size = 0.5, alpha = 0.8,position = position_dodge(0.8)) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1))+
  scale_fill_manual(values = tissue_cols) +
  stat_compare_means(method = "wilcox.test", label = "p.format", hide.ns = FALSE, size=3, label.y.npc = 0.95) +
  theme_bw(base_size =10) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid = element_blank()
  )

ggsave(
  filename = "Sample_neighborhood10_boxplot_format.pdf",  # 支持 PDF/PNG/TIFF 等格式
  plot = p,
  device = "pdf",            # 保存为 PDF
  width = 15,               # A4 宽度 (mm)
  height = 5,              # A4 高度 (mm)
  #units = "mm",              # 单位设为毫米
  dpi = 300                  # 分辨率（默认 300 DPI）
)

### 4e-g
# =========================================================
# Full KMeans pipeline for Tumor CN analysis
# Steps:
# 1) Read input workbook
# 2) CLR transform for compositional CN data
# 3) KMeans k-selection (k = 2:6)
# 4) Final KMeans clustering using selected k
# 5) Clinical association analysis
# 6) Survival analysis
# 7) Plotting
# 8) Export tables and figures
# =========================================================

required_pkgs <- c(
  "readxl", "dplyr", "tidyr", "ggplot2", "cluster", "mclust",
  "survival", "survminer", "openxlsx", "pheatmap"
)

missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    paste0(
      "Missing required packages: ",
      paste(missing_pkgs, collapse = ", "),
      ". Please install them first."
    )
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cluster)
  library(mclust)
  library(survival)
  library(survminer)
  library(openxlsx)
  library(pheatmap)
  library(broom)
})

input_file <- "~/Downloads/cn_kmeans/tumor_cn_kmeans3_input.xlsx"
out_prefix <- "~/Downloads/cn_kmeans/outputs5/kmeans_full_pipeline"

# ---------------------------------------------------------
# Read input
# ---------------------------------------------------------
clinical_data <- read_excel(input_file, sheet = "clinical_data")
clinical_data <- clinical_data[clinical_data$Sample %in% need_sample,]
tumor_cn_ratio <- read_excel(input_file, sheet = "tumor_cn_ratio")
tumor_cn_ratio <- tumor_cn_ratio[tumor_cn_ratio$Sample %in% need_sample,]
merged_input <- read_excel(input_file, sheet = "merged_input")
merged_input <- merged_input[merged_input$Sample %in% need_sample,]

cn_cols <- grep("^CN[0-9]+$", colnames(tumor_cn_ratio), value = TRUE)
cn_cols <- cn_cols[order(as.numeric(sub("CN", "", cn_cols)))]

# ---------------------------------------------------------
# CLR transform
# ---------------------------------------------------------
X <- as.matrix(tumor_cn_ratio[, cn_cols])
positive_min <- min(X[X > 0], na.rm = TRUE)
X_pseudo <- X + positive_min / 2
X_pseudo <- X_pseudo / rowSums(X_pseudo)

clr_transform <- function(x) {
  lx <- log(x)
  lx - mean(lx)
}
X_clr <- t(apply(X_pseudo, 1, clr_transform))
rownames(X_clr) <- tumor_cn_ratio$Sample

# ---------------------------------------------------------
# Helper: multi-group log-rank p-value
# ---------------------------------------------------------
calc_logrank_p <- function(time, status, group) {
  dat <- data.frame(time = time, status = status, group = group)
  dat <- dat[complete.cases(dat), , drop = FALSE]
  if (length(unique(dat$group)) < 2) return(NA_real_)
  fit <- survdiff(Surv(time, status) ~ group, data = dat)
  1 - pchisq(fit$chisq, df = length(fit$n) - 1)
}

# ---------------------------------------------------------
# K-selection for KMeans only
# Criteria:
# 1) silhouette
# 2) subsample stability by ARI
# 3) survival separation by log-rank p
# 4) minimum cluster size penalty
# Lower composite score = better
# ---------------------------------------------------------
set.seed(0)
k_grid <- 2:6
#k_grid <- 3
benchmark_list <- list()

for (k in k_grid) {
  km0 <- kmeans(X_clr, centers = k, nstart = 50)
  labels0 <- km0$cluster

  sil <- cluster::silhouette(labels0, dist(X_clr))
  sil_mean <- mean(sil[, "sil_width"], na.rm = TRUE)

  surv_p <- calc_logrank_p(
    time = merged_input$time,
    status = merged_input$status,
    group = labels0
  )

  cluster_sizes <- as.integer(table(labels0))
  min_cluster_size <- min(cluster_sizes)

  set.seed(100 + k)
  ari_vec <- c()
  n <- nrow(X_clr)
  for (b in 1:20) {
    idx <- sort(sample(seq_len(n), size = floor(0.8 * n), replace = FALSE))
    km_sub <- kmeans(X_clr[idx, , drop = FALSE], centers = k, nstart = 20)
    ari_vec <- c(ari_vec, mclust::adjustedRandIndex(labels0[idx], km_sub$cluster))
  }
  stability_ari <- mean(ari_vec, na.rm = TRUE)
  wcss <- km0$tot.withinss

  benchmark_list[[as.character(k)]] <- data.frame(
    Method = "KMeans_CLR",
    k = k,
    Silhouette = sil_mean,
    Logrank_p = surv_p,
    Stability_ARI = stability_ari,
    Min_cluster_size = min_cluster_size,
    Cluster_sizes = paste(cluster_sizes, collapse = ", "),
    Wcss = wcss,
    stringsAsFactors = FALSE
  )
}

benchmark_df <- bind_rows(benchmark_list) %>%
  mutate(
    Rank_silhouette = rank(Silhouette, ties.method = "min"),
    Rank_survival   = rank(Logrank_p, ties.method = "min"),
    Rank_stability  = rank(-Stability_ARI, ties.method = "min"),
    Size_penalty    = case_when(
      Min_cluster_size >= 10 ~ 0,
      Min_cluster_size >= 7  ~ 1,
      TRUE ~ 2
    ),
    Composite_score = Rank_silhouette + Rank_survival + Rank_stability + Size_penalty
  ) %>%
  arrange(Composite_score, Logrank_p, desc(Stability_ARI))

selected_k <- benchmark_df$k[1]

# ---------------------------------------------------------
# Final KMeans clustering using selected k
# ---------------------------------------------------------
set.seed(0)
km_final <- kmeans(X_clr, centers = selected_k, nstart = 50)

cluster_df <- tumor_cn_ratio[, c("Sample")]
cluster_df$Cluster_raw <- km_final$cluster

cluster_means <- merged_input %>%
  mutate(Cluster_raw = km_final$cluster) %>%
  group_by(Cluster_raw) %>%
  summarise(across(all_of(cn_cols), mean, na.rm = TRUE), .groups = "drop")

make_subtype_name <- function(v, cluster_id) {
  top_cn <- names(sort(v, decreasing = TRUE))[1]
  paste0("S", cluster_id, "_", top_cn, "-dominant")
}

cluster_names <- sapply(seq_len(nrow(cluster_means)), function(i) {
  make_subtype_name(as.numeric(cluster_means[i, cn_cols]), cluster_means$Cluster_raw[i])
})

name_map <- setNames(cluster_names, cluster_means$Cluster_raw)

cluster_df$Final_Tumor_subtype <- name_map[as.character(cluster_df$Cluster_raw)]

result_df <- merged_input %>%
  left_join(cluster_df, by = "Sample")

result_df$Subtype <- paste0("Ecotype_", result_df$Cluster_raw)
result_df$Subtype <- factor(result_df$Subtype, levels=c("Ecotype_1", "Ecotype_2", "Ecotype_3"))
# ---------------------------------------------------------
# CN feature summary
# ---------------------------------------------------------
feature_summary <- result_df %>%
  group_by(Final_Tumor_subtype) %>%
  summarise(across(all_of(cn_cols), mean, na.rm = TRUE), .groups = "drop")

cn_diff_list <- lapply(cn_cols, function(v) {
  dat <- result_df[, c("Final_Tumor_subtype", v)]
  dat <- dat[complete.cases(dat), , drop = FALSE]
  kw <- kruskal.test(dat[[v]] ~ dat$Final_Tumor_subtype)
  out <- data.frame(CN = v, P_value = kw$p.value, stringsAsFactors = FALSE)
  for (stype in unique(dat$Final_Tumor_subtype)) {
    vv <- dat[dat$Final_Tumor_subtype == stype, v, drop = TRUE]
    out[[paste0("Mean_", stype)]] <- mean(vv, na.rm = TRUE)
    out[[paste0("Median_", stype)]] <- median(vv, na.rm = TRUE)
  }
  out
}) %>% bind_rows() %>%
  mutate(FDR = p.adjust(P_value, method = "BH")) %>%
  arrange(P_value)

# ---------------------------------------------------------
# Clinical association analysis
# ---------------------------------------------------------
clinical_vars <- c(
  "Gender", "Age", "Grade", "Grade_2", "Tumor_size (cm)",
  "T", "T_stage", "T_stage_2", "N", "N_stage",
  "LNM (0,no;1,yes)", "AJCC_Stage", "AJCC_Stage_2", "pStage",
  "mean_distance", "mean_distance_gorup"
)

clinical_assoc <- lapply(clinical_vars, function(v) {
  if (!v %in% colnames(result_df)) return(NULL)
  dat <- result_df[, c("Final_Tumor_subtype", v)]
  dat <- dat[complete.cases(dat), , drop = FALSE]
  if (nrow(dat) == 0 || dplyr::n_distinct(dat[[v]]) < 2) return(NULL)

  if (is.numeric(dat[[v]]) && dplyr::n_distinct(dat[[v]]) > 5) {
    kw <- kruskal.test(dat[[v]] ~ dat$Final_Tumor_subtype)
    data.frame(
      Clinical_variable = v,
      Test = "Kruskal-Wallis",
      N = nrow(dat),
      P_value = kw$p.value,
      stringsAsFactors = FALSE
    )
  } else {
    tab <- table(dat$Final_Tumor_subtype, dat[[v]])
    chi <- suppressWarnings(chisq.test(tab))
    data.frame(
      Clinical_variable = v,
      Test = "Chi-square",
      N = nrow(dat),
      P_value = chi$p.value,
      stringsAsFactors = FALSE
    )
  }
}) %>% bind_rows()

clinical_assoc <- clinical_assoc %>%
  mutate(FDR = p.adjust(P_value, method = "BH")) %>%
  arrange(P_value)

# ---------------------------------------------------------
# Survival analysis
# ---------------------------------------------------------
surv_df <- result_df %>%
  filter(!is.na(time), !is.na(status))
surv_df$Age_2 <- as.numeric(surv_df$Age > 64)
surv_df$Tumor_subtype <- as.numeric(surv_df$Cluster_raw)
surv_df$AJCC_Stage_2 <- as.numeric(surv_df$AJCC_Stage_2)
km_fit <- survfit(Surv(time, status) ~ Final_Tumor_subtype, data = surv_df)
surv_diff <- survdiff(Surv(time, status) ~ Final_Tumor_subtype, data = surv_df)
global_logrank_p <- 1 - pchisq(surv_diff$chisq, df = length(surv_diff$n) - 1)
model <- coxph(Surv(time, status) ~ Age_2 + Grade_2 + AJCC_Stage_2 + Tumor_subtype, data = surv_df) 

cox_fit <- coxph(Surv(time, status) ~ Final_Tumor_subtype, data = surv_df)
cox_summary <- summary(cox_fit)

cox_table <- data.frame(
  Variable = rownames(cox_summary$coefficients),
  HR = cox_summary$coefficients[, "exp(coef)"],
  CI95_low = cox_summary$conf.int[, "lower .95"],
  CI95_high = cox_summary$conf.int[, "upper .95"],
  P_value = cox_summary$coefficients[, "Pr(>|z|)"],
  stringsAsFactors = FALSE
)

surv_counts <- surv_df %>%
  group_by(Final_Tumor_subtype) %>%
  summarise(
    N_survival = n(),
    Events = sum(status, na.rm = TRUE),
    .groups = "drop"
  )

cn_cols <- paste0("CN", 1:16)
stopifnot(all(cn_cols %in% colnames(surv_df)))
cnsurv_list <- list()
for (col in cn_cols) {
  temp_df <- surv_df %>%
    select(time, status, !!sym(col)) %>%
    na.omit()
  if (nrow(temp_df) == 0 || sd(temp_df[[col]], na.rm = TRUE) == 0) {
    warning(paste("Skipping", col, ": no valid data"))
    next
  }
  med_val <- median(temp_df[[col]], na.rm = TRUE)
  temp_df$group <- as.numeric(temp_df[[col]] > med_val)
  fit <- coxph(Surv(time, status) ~ group, data = temp_df)
  tidy_res <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == "group")  %>%
    select(
      HR = estimate,
      CI95_low = conf.low,
      CI95_high = conf.high,
      P.value = p.value
    )
  tidy_res$CN <- col
  cnsurv_list[[col]] <- tidy_res
}

cnsurv_table <- bind_rows(cnsurv_list) %>%
  mutate(
    P.value = round(P.value, 4),
    HR = round(HR, 3),
    CI95_low = round(CI95_low, 3),
    CI95_high = round(CI95_high, 3)
  ) %>%
  arrange(P.value)

cnsurv_table <- cnsurv_table[,c("CN", "HR", "CI95_low", "CI95_high", "P.value")]

# ---------------------------------------------------------
# Plots
# ---------------------------------------------------------
p_sil <- ggplot(benchmark_df, aes(x = k, y = Silhouette)) +
  geom_line() +
  geom_point(size = 3) +
  geom_vline(xintercept = selected_k, linetype = 2) +
  scale_x_continuous(breaks = k_grid) +
  labs(
    title = "KMeans k selection by silhouette",
    x = "k",
    y = "Silhouette score"
  ) +
  theme_bw(base_size = 12)
ggsave(filename=paste0(out_prefix, "_silhouette.pdf"), plot=p_sil, width = 7, height = 5)

p_elbow <- ggplot(benchmark_df, aes(x = k, y = Wcss)) +
  geom_line() + 
  geom_point(size = 3) +
  geom_vline(xintercept = selected_k, linetype = 2) +
  scale_x_continuous(breaks = k_grid) +
  labs(title = "Elbow Method for Optimal k",
       x = "k",
       y = "Within-Cluster Sum of Squares (WCSS)") +
  theme_bw(base_size = 12)
ggsave(filename=paste0(out_prefix, "_elbow.pdf"), plot=p_elbow, width = 7, height = 5)

p_survsel <- ggplot(benchmark_df, aes(x = k, y = -log10(Logrank_p))) +
  geom_line() +
  geom_point(size = 3) +
  geom_vline(xintercept = selected_k, linetype = 2) +
  scale_x_continuous(breaks = k_grid) +
  labs(
    title = "KMeans k selection by survival separation",
    x = "k",
    y = "-log10(log-rank p)"
  ) +
  theme_bw(base_size = 12)
ggsave(filename=paste0(out_prefix, "_survival_selection.pdf"), plot=p_survsel, width = 7, height = 5)

pca <- prcomp(X_clr, center = TRUE, scale. = FALSE)
pca_df <- data.frame(
  Sample = tumor_cn_ratio$Sample,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
) %>% left_join(cluster_df, by = "Sample")

pca_df$Subtype <- paste0("Ecotype_", pca_df$Cluster_raw)
pca_df$Subtype <- factor(pca_df$Subtype, levels=c("Ecotype_1", "Ecotype_2", "Ecotype_3"))

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = Subtype)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(
    title = paste0("Final KMeans CLR subtype PCA (k = ", selected_k, ")"),
    x = paste0("PC1 (", round(summary(pca)$importance[2,1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca)$importance[2,2] * 100, 1), "%)")
  ) +
  theme_bw(base_size = 12)
ggsave(filename=paste0(out_prefix, "_pca.pdf"), plot=p_pca, width = 7, height = 6)

km_plot <- ggsurvplot(
  km_fit,
  data = surv_df,
  pval = TRUE,
  conf.int = TRUE, 
  fun = "pct",
  risk.table = TRUE,
  risk.table.height = 0.25,
  linetype = "strata",
  legend = "top",
  legend.labs = c("Ecotype_1", "Ecotype_2", "Ecotype_3"),
  legend.title = "Subtype",
  title = paste0("Final subtype survival (selected k = ", selected_k,
                 ", log-rank p = ", signif(global_logrank_p, 3), ")")
)
pdf(paste0(out_prefix, "_km.pdf"), width = 8, height = 6)
#ggsave(filename=paste0(out_prefix, "_km.pdf"), plot=km_plot, width = 8, height = 8)
print(km_plot)
dev.off()

p_for <- ggforest(model,refLabel = '1', noDigits =3)
ggsave(filename=paste0(out_prefix, "_forest.pdf"), plot=p_for, width = 7, height = 8)

heat_df <- result_df %>%
  select(Sample, Subtype, all_of(cn_cols)) %>%
  arrange(Subtype)

heat_mat <- as.matrix(heat_df[, cn_cols])
rownames(heat_mat) <- heat_df$Sample
ann <- data.frame(Subtype = heat_df$Subtype)
rownames(ann) <- heat_df$Sample
ann_colors = list(Subtype = c(Ecotype_1 = "#EB9187",Ecotype_2 = "#6DC266", Ecotype_3 ="#84A9F9"))
pdf(paste0(out_prefix, "_heatmap.pdf"), width = 16, height = 8)
pheatmap(
  t(heat_mat),
  annotation_col = ann,
  annotation_colors = ann_colors,
  scale = "none",
  cluster_cols = FALSE,
  main = paste0("Tumor CN composition heatmap (KMeans, k = ", selected_k, ")")
)
dev.off()

comp_df <- result_df %>%
  group_by(Subtype) %>%
  summarise(across(all_of(cn_cols), mean, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(cols = all_of(cn_cols), names_to = "CN", values_to = "Mean_ratio")

mycol <- c("#4E79A7","#A0CBE8","#90EE90","#FFBE7D","#F28E2B","#8CD17D","#B6992D","#F1CE63","#FF0000","#86BCB6","#FF5759","#FF9D9A","#79706E","#BAB0AC", "#D37295", "#FABFD2")

comp_df$CN <- factor(comp_df$CN, levels=unique(comp_df$CN))
p_comp <- ggplot(comp_df, aes(Subtype, Mean_ratio * 100, fill = CN)) +
  geom_col() +
  labs(
    title = paste0("Average CN composition of final KMeans subtypes (k = ", selected_k, ")"),
    x = "Subtype",
    y = "Mean CN proportion (%)"
  ) +
  theme_bw(base_size = 12) + scale_fill_manual(values = mycol)
ggsave(filename=paste0(out_prefix, "_composition.pdf"), plot=p_comp, width = 9, height = 5)

# ---------------------------------------------------------
# Export results
# ---------------------------------------------------------
wb <- createWorkbook()

addWorksheet(wb, "README")
writeData(wb, "README", data.frame(
  Item = c(
    "Input file",
    "Method",
    "Selected k",
    "Selection rule",
    "Global log-rank p"
  ),
  Value = c(
    input_file,
    "KMeans on CLR-transformed CN composition",
    selected_k,
    "Composite ranking using silhouette + stability + survival + cluster-size penalty",
    signif(global_logrank_p, 4)
  )
))

addWorksheet(wb, "Benchmark_k_selection")
writeData(wb, "Benchmark_k_selection", benchmark_df)

addWorksheet(wb, "Subtype_assignments")
writeData(wb, "Subtype_assignments", result_df)

addWorksheet(wb, "Subtype_feature_summary")
writeData(wb, "Subtype_feature_summary", feature_summary)

addWorksheet(wb, "Subtype_CN_differences")
writeData(wb, "Subtype_CN_differences", cn_diff_list)

addWorksheet(wb, "Clinical_association")
writeData(wb, "Clinical_association", clinical_assoc)

addWorksheet(wb, "Survival_summary")
writeData(wb, "Survival_summary", surv_counts)

addWorksheet(wb, "Cox_results")
writeData(wb, "Cox_results", cox_table)

addWorksheet(wb, "Cox_perCN")
writeData(wb, "Cox_perCN", cnsurv_table)

saveWorkbook(wb, paste0(out_prefix, "_results.xlsx"), overwrite = TRUE)

# ---------------------------------------------------------
# Console summary
# ---------------------------------------------------------
cat("Full KMeans pipeline finished.\n")
cat("Selected k =", selected_k, "\n")
cat("Global log-rank p =", signif(global_logrank_p, 4), "\n\n")
cat("Main outputs:\n")
cat(" - ", paste0(out_prefix, "_results.xlsx"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_silhouette.png"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_survival_selection.png"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_pca.png"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_km.png"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_heatmap.png"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_composition.png"), "\n", sep = "")


