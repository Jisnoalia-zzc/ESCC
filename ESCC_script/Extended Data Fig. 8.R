#Extended  Data Fig8a -----
library(scRNAtoolVis)
load("CD8_C6_CD39.Rdata");deg$cluster<-'CD8_C6_PDCD1';deg$gene<-rownames(deg);com<-deg
load("CD4_C7_OX40.Rdata");deg$cluster<-'CD4_C7_OX40';deg$gene<-rownames(deg);com<-rbind(com,deg)
load("FB_C3_COL1A1.Rdata");deg$cluster<-'FB_C3_COL1A1';deg$gene<-rownames(deg);com<-rbind(com,deg)
load("Mac_C2_SPP1.Rdata");deg$cluster<-'Mac_C2_SPP1';deg$gene<-rownames(deg);com<-rbind(com,deg)
load("Endo_C3_RGCC.Rdata");deg$cluster<-'Endo_C3_RGCC';deg$gene<-rownames(deg);com<-rbind(com,deg)


mygene <- c('TOX','PDCD1','LAYN','GEM','ENTPD1','DUSP4','CXCL13',
            'CXCL8','TGFB1','TGFB3','WNT2','WNT5A','FAP','MMP14','COL12A1','COL5A1','COL11A1',
            'COL3A1','COL1A1','MMP1','MMP11','POSTN','ANTXR1','CXCL12','TGFBR3','IL6',
            'SPP1','IL1RN','OLR1','MARCO','EREG','CCL20','MMP7','VEGFA','MMP12','MMP9','CXCL8','CCL2','AREG',
            'C1QB','C1QA','C1QC','CCL18','IL18',
            'FOXP3','TNFRSF18','TNFRSF4','TNFRSF9','TNFRSF1B','IL32','TIGIT','BATF','DUSP4','LAIR2','CTLA4','IL2RA',
            "PLVAP", "COL4A1", "COL4A2", "HSPG2", "VWF", "IGFBP7", "PECAM1", "SERPINE1", "SPARC", "INSR"
            )
mygene<-mygene[!duplicated(mygene)]
com = com[com$avg_log2FC >0,]
com$cluster <- factor(com$cluster,
                      levels = c("CD8_C6_PDCD1","CD4_C7_OX40","FB_C3_COL1A1","Endo_C3_RGCC","Mac_C2_SPP1"))
plot<-jjVolcano(diffData = com, myMarkers = mygene,size  = 4,
                fontface = 'italic',log2FC.cutoff = 0.25,
                
                tile.col = c("#5aa1a3","#20b2aa","#ba6c35","#7d5599","#c9c193"),
                aesCol = c('blue','red'),pSize = 1)
plot


#Extended  Data Fig8b -----
rm(list=ls());gc()
library(GseaVis)
library(clusterProfiler)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))
library(fgsea)
library(msigdbr)
msigdbr_collections()
human_H<-msigdbr(species = "human", category = "H")
human_C2<-msigdbr(species = "human", category = "C2")
human_C5<-msigdbr(species = "human", category = "C5",subcategory = "BP")

human_C2<-human_C2[human_C2$gs_subcat%in%c('CP:BIOCARTA','CP:KEGG','CP:REACTOME'),]
human<-rbind(human_C2,human_H,human_C5)
gmt<-human %>% dplyr::select(gs_name,gene_symbol)
deg<- deg_fb
rank<-deg[,c('Gene','avg_log2FC')]
geneList<-rank$avg_log2FC
names(geneList)=rank$Gene 
geneList=sort(geneList,decreasing = TRUE) 
gseaRes <- GSEA(geneList = geneList,TERM2GENE = gmt, minGSSize = 0, maxGSSize = 1000, pvalueCutoff = 1, pAdjustMethod = "BH",verbose = FALSE)
geneSetID = c('KEGG_ECM_RECEPTOR_INTERACTION','KEGG_FOCAL_ADHESION','REACTOME_SIGNALING_BY_MET','REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX')
gs<-gseaRes[gseaRes@result$Description%in%geneSetID,]
gs->gs_fb
rank<-deg[,c('Gene','avg_log2FC')]
geneList<-rank$avg_log2FC
names(geneList)=rank$Gene 
geneList=sort(geneList,decreasing = TRUE) 
gseaRes <- GSEA(geneList = geneList,TERM2GENE = gmt, minGSSize = 0, maxGSSize = 1000, pvalueCutoff = 1, pAdjustMethod = "BH",verbose = FALSE)
geneSetID = c('GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY','GOBP_LEUKOCYTE_DIFFERENTIATION','GOBP_LEUKOCYTE_MEDIATED_IMMUNITY')
gs<-gseaRes[gseaRes@result$Description%in%geneSetID,]
gs->gs_cd8
rank<-deg[,c('gene','avg_log2FC')]
geneList<-rank$avg_log2FC
names(geneList)=rank$gene 
geneList=sort(geneList,decreasing = TRUE) 
gseaRes <- GSEA(geneList = geneList,TERM2GENE = gmt, minGSSize = 0, maxGSSize = 1000, pvalueCutoff = 1, pAdjustMethod = "BH",verbose = FALSE)
geneSetID = c('GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY','GOBP_REGULATION_OF_T_CELL_DIFFERENTIATION')
gs<-gseaRes[gseaRes@result$Description%in%geneSetID,]
gs->gs_cd4
deg<-deg_mac
rank<-deg[,c('gene','avg_log2FC')]
geneList<-rank$avg_log2FC
names(geneList)=rank$gene 
geneList=sort(geneList,decreasing = TRUE) 
gseaRes <- GSEA(geneList = geneList,TERM2GENE = gmt, minGSSize = 0, maxGSSize = 1000, pvalueCutoff = 1, pAdjustMethod = "BH",verbose = FALSE)
geneSetID = c('HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION','HALLMARK_ANGIOGENESIS')
gs<-gseaRes[gseaRes@result$Description%in%geneSetID,]
gs->gs_mac
deg<-deg_endo
rank<-deg[,c('gene','avg_log2FC')]
geneList<-rank$avg_log2FC
names(geneList)=rank$gene 
geneList=sort(geneList,decreasing = TRUE) 
gseaRes <- GSEA(geneList = geneList,TERM2GENE = gmt, minGSSize = 0, maxGSSize = 1000, pvalueCutoff = 1, pAdjustMethod = "BH",verbose = FALSE)
geneSetID = c("REACTOME_VEGF_LIGAND_RECEPTOR_INTERACTIONS","REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS","REACTOME_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES")
gs<-gseaRes[gseaRes@result$Description%in%geneSetID,]
gs->gs_endo
#####################
gs_endo$class<-"Endo_C3_RGCC"
gs_cd4$class<-"CD4_C7_OX40"
gs_cd8$class<-"CD8_C6_CDPD1"
gs_mac$class<-"Mac_C2_SPP1"
gs_fb$class<-"FB_C3_COL1A1"

dat<-rbind(gs_cd8,gs_cd4,gs_endo,gs_fb,gs_mac)
dat<-dat[,c('ID','NES','pvalue','class')]

library(forcats)
library(ggplot2)
dat$ID <- as.factor(dat$ID)
dat$ID <- fct_inorder(dat$ID)
dat$class <- as.factor(dat$class)
dat$class <- fct_inorder(dat$class)
dat$log10p <- -log(dat$pvalue+0.001,10)

ggplot(dat, aes(class, ID)) +
  geom_point(aes(fill=NES, size=-log10(pvalue)),shape=21,color='black')+theme_bw()+
  scale_fill_gradientn(colours = rev(getPalette(10)),limits=c(1,2.7),)+

  theme(axis.text.x = element_text(size=14,angle = 90,hjust=1,colour = "black"),panel.background = element_blank(),
        axis.text.y = element_text(size=14,colour = "black"),panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.7),legend.text = element_text(size = 14),legend.title = element_text(size=10),
        legend.background = element_blank(),legend.key = element_blank(),title = element_text(size = 14))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=1))
ggsave(file="/realspace/project/proj_ESCC_STW_ZWM_2022_01/liuliqiu/figs1/function.pdf",height = 6,width = 9)

#### E8C1
library(ggplot2)
library(dplyr)
library(data.table)

setwd('C:/Users/Ron/Desktop/Figs/E8/')
file_stats <- fread('C:/Users/Ron/Desktop/Figs/E8/E8_c1_file_stats.csv')
tissue_colors <- c(
  "TF" = "#8B0000",
  "LpF" = "#CD5C5C", 
  "AF" = "#034272",
  "NF" = "#0d6cc0",
  "LnF" = "#87CEEB"
)

present_tissues <- unique(file_stats$tissue)
available_colors <- tissue_colors[names(tissue_colors) %in% present_tissues]

file_stats <- file_stats %>%
  mutate(
    p_adj = ifelse(p_value == 0, .Machine$double.xmin, p_value),
    neg_log_p = -log10(p_adj)
  )

finite_p <- file_stats$neg_log_p[is.finite(file_stats$neg_log_p)]
y_max <- if(length(finite_p) > 0) max(finite_p, na.rm = TRUE) else 10
y_max <- if(is.infinite(y_max) | is.na(y_max)) 10 else y_max
y_limit <- y_max * 1.2

file_stats <- file_stats %>%
  mutate(neg_log_p_plot = ifelse(is.infinite(neg_log_p), y_limit * 0.9, neg_log_p))

p <- ggplot(file_stats, aes(x = diff, y = neg_log_p_plot, color = tissue)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
             color = "red", alpha = 0.5, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", 
             color = "gray50", alpha = 0.5, linewidth = 0.8) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = available_colors, name = "Tissue Type") +
  scale_y_continuous(
    limits = c(0, y_limit),
    expand = expansion(mult = c(0, 0.1)),
    name = "-log10(p-value)"
  ) +
  labs(
    title = "MIMER vs Random Combinations",
    subtitle = paste0(sum(file_stats$Significance == "Yes"), "/", nrow(file_stats), 
                     " files (", round(mean(file_stats$Significance == "Yes") * 100, 1), 
                     "%) significant at p < 0.05"),
    x = "Difference (MIMER mean - Random mean)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  )

significant_points <- file_stats %>%
  filter(p_value < 0.05, abs(diff) > quantile(abs(diff), 0.75, na.rm = TRUE))

if(nrow(significant_points) > 0) {
  p <- p + 
    geom_text(data = significant_points,
              aes(label = file), vjust = -0.5, hjust = 0.5, 
              size = 3, color = "black", check_overlap = TRUE)
}

print(p)

####
# E8C2
library(data.table)
library(tibble)
library(dplyr)
library(ggplot2)

setwd('C:/Users/Ron/Desktop/Figs/E8/')
plot_data <- fread('C:/Users/Ron/Desktop/Figs/E8/E8_c2_Moran_plotdata.csv')
meta.data.tissue <- fread('C:/Users/Ron/Desktop/Figs/E8/E8_c2_meta.data.tissue.csv')

file_tissue_mark <- meta.data.tissue %>% 
  unique() %>%
  mutate(tissue = ifelse(tissue %in% c('TF', 'NF', 'LnF', 'LpF', 'AF'), 
                         tissue, "mixed"))

p_threshold = 0.05
point_size = 3
alpha = 0.7

tissue_colors <- c(
"TF" = "#8B0000","LpF" = "#CD5C5C", "AF" = "#034272","NF" = "#0d6cc0","LnF" = "#87CEEB","mixed" = "#808080")

present_tissues <- unique(plot_data$tissue)
available_colors <- tissue_colors[names(tissue_colors) %in% present_tissues]

p <- ggplot(plot_data, aes(x = diff_obs_perm, y = neg_log10_p)) +
geom_hline(yintercept = -log10(p_threshold), 
            linetype = "dashed", color = "red", alpha = 0.5) +
geom_vline(xintercept = 0, 
            linetype = "dashed", color = "gray", alpha = 0.5) +
geom_point(aes(color = tissue),  shape = 19,  size = point_size, alpha = alpha) + scale_color_manual(values = available_colors) +
labs(
    title = "Moran's I: Observed vs Permutation",
    x = "Difference (Observed - Permutation Mean)",
    y = "-log10(p-value)",
    color = "Tissue Type"
) + theme_minimal() +
theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
)

stats_df <- plot_data %>% group_by(tissue) %>% summarise(
    total = n(),
    significant = sum(Significant == "Yes"),
    nonsignificant = sum(Significant == "No"),
    .groups = "drop"
) %>% mutate( label = sprintf("%s: Yes=%d, No=%d (Total=%d)", tissue, significant, nonsignificant, total))
p



#### E8de
library(data.table)
library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)
library(tibble)

options(repr.plot.width = 12, repr.plot.height = 9)

setwd('C:/Users/Ron/Desktop/Figs/E8/')
mimer_median_enrich_development_ME <- fread('C:/Users/Ron/Desktop/Figs/E8/E8_def_mimer_median_enrich_development_ME.csv')

cell_type_mapping <- c(
  "CD8_C6_CD39_above_median" = "CD8_C6_CD39",
  "CD4_C7_OX40_above_median" = "CD4_C7_OX40", 
  "Mac_C2_SPP1_above_median" = "Mac_C2_SPP1",
  "FB_C3_COL1A1_above_median" = "FB_C3_COL1A1",
  "Endo_C3_RGCC_above_median" = "Endo_C3_RGCC"
)

tissue_map <- c(A = "Adj", N = "N", T = "T", Ln = "LnN", Lp = "LnP")
tissue_levels <- c("N", "Adj", "T", "LnP", "LnN")
tissue_colors <- c(
  "N" = "#6888B9",
  "Adj" = "#9AB292", 
  "T" = "#C09DB0",
  "LnP" = "#CC834F",
  "LnN" = "#DCC13A"
)

cell_type_levels <- c("MIMER", "Mixed", "CD8_C6_CD39" ,"CD4_C7_OX40", "Mac_C2_SPP1",
                      "FB_C3_COL1A1", "Endo_C3_RGCC","Other")

cell_type_colors <- c(
  "MIMER" = "#A52A2A",
  "Mixed" = "#FF6B6B",
  "CD8_C6_CD39" = "#569c9d",
  "CD4_C7_OX40" = "#1fb1aa",
  "Mac_C2_SPP1" = "#c2bc8d",
  "FB_C3_COL1A1" = "#b96c35",
  "Endo_C3_RGCC" = "#775095",
  "Other" = "#bebebe"
)

stage_colors <- c("Nor" = "#1B9E77", "Hyp" = "#D95F02",  "MiD" = "#666666",  "MoD" = "#E7298A",  "SD&CA" = "#66A61E",
  "ICA" = "#E6AB02",  "MCA" = "#A6761D",  "Nor-ME" = "#1B9E77",  "Hyp-ME" = "#D95F02",  "MiD-ME" = "#666666",  "MoD-ME" = "#E7298A",
  "SD&CA-ME" = "#66A61E",  "ICA-ME" = "#E6AB02",  "MCA-ME" = "#A6761D")

define_cell_type <- function(dt) {
  dt[, {
    if (mimer == "MIMER") {
      cell_type <- "MIMER"
    } else {
      markers <- c(CD8_C6_CD39_above_median, CD4_C7_OX40_above_median,
                  Mac_C2_SPP1_above_median, FB_C3_COL1A1_above_median,
                  Endo_C3_RGCC_above_median)
      true_count <- sum(markers)
      
      if (true_count >= 2) {
        cell_type <- "Mixed"
      } else if (true_count == 1) {
        true_idx <- which(markers)
        marker_names <- c("CD8_C6_CD39_above_median", "CD4_C7_OX40_above_median",
                         "Mac_C2_SPP1_above_median", "FB_C3_COL1A1_above_median",
                         "Endo_C3_RGCC_above_median")
        cell_type <- cell_type_mapping[marker_names[true_idx]]
      } else {
        cell_type <- "Other"
      }
    }
    .(cell_type = cell_type)
  }, by = .(cellID, file, patient, tissue_SPT)]
}

analyze_and_plot <- function(dt_subset, title, filename_prefix, group_var) {
  dt_subset[, tissue := tissue_map[tissue_SPT]]
  dt_subset <- dt_subset[tissue %in% tissue_levels]
  
  cell_types <- define_cell_type(dt_subset)
  dt_subset <- merge(dt_subset, cell_types, by = c("cellID", "file", "patient", "tissue_SPT"))
  dt_subset[, tissue := factor(tissue, levels = tissue_levels)]
  dt_subset[, cell_type := factor(cell_type, levels = cell_type_levels)]
  
  sample_comp <- dt_subset[, .(count = .N), by = .(file, tissue, patient, cell_type)]
  sample_totals <- sample_comp[, .(total = sum(count)), by = .(file, tissue, patient)]
  sample_comp <- merge(sample_comp, sample_totals, by = c("file", "tissue", "patient"))
  sample_comp[, percentage := count / total * 100]
  
  sample_comp[, sample_id := paste(patient, file, tissue, sep = "_")]
  mimer_percent <- sample_comp[cell_type == "MIMER", .(sample_id, mimer_percent = percentage)]
  sample_comp <- merge(sample_comp, mimer_percent, by = "sample_id", all.x = TRUE)
  sample_comp[is.na(mimer_percent), mimer_percent := 0]
  
  stage_comp <- dt_subset[, .(count = .N), by = .(file, tissue, patient, get(group_var))]
  setnames(stage_comp, "get", "stage")
  stage_totals <- stage_comp[, .(total = sum(count)), by = .(file, tissue, patient)]
  stage_comp <- merge(stage_comp, stage_totals, by = c("file", "tissue", "patient"))
  stage_comp[, percentage := count / total * 100]
  stage_comp[, sample_id := paste(patient, file, tissue, sep = "_")]
  
  stage_comp <- merge(stage_comp, unique(sample_comp[, .(sample_id, mimer_percent)]), by = "sample_id")
  
  stage_levels <- names(stage_colors)[names(stage_colors) %in% unique(stage_comp$stage)]
  if(length(stage_levels) > 0) {
    stage_comp[, stage := factor(stage, levels = stage_levels)]
  }
  
  sample_order <- unique(sample_comp[order(-mimer_percent)]$sample_id)
  sample_comp[, sample_id := factor(sample_id, levels = sample_order)]
  stage_comp[, sample_id := factor(sample_id, levels = sample_order)]
  
  sample_comp <- sample_comp[!is.na(sample_id)]
  stage_comp <- stage_comp[!is.na(sample_id)]
  
  sample_info <- unique(sample_comp[, .(sample_id, tissue)])
  sample_info[, tissue_color := tissue_colors[as.character(tissue)]]
  
  p_cell <- ggplot(sample_comp, aes(x = sample_id, y = percentage, fill = cell_type)) +
    geom_bar(stat = "identity", position = "stack", width = 0.8) +
    scale_fill_manual(values = cell_type_colors, name = "Cell Type") +
    labs(title = paste(title, "- Cell Type Composition"),
         x = "", y = "Percentage (%)") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100.1))
  
  p_stage <- ggplot(stage_comp, aes(x = sample_id, y = percentage, fill = stage)) +
    geom_bar(stat = "identity", position = "stack", width = 0.8) +
    scale_fill_manual(values = stage_colors, name = group_var) +
    labs(title = paste(title, "-", group_var, "Type Composition"),
         x = "", y = "Percentage (%)") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100.1))
  
  p_tissue <- ggplot(sample_info, aes(x = sample_id, y = 1, fill = tissue)) +
    geom_tile(color = "black", size = 0.2) +
    scale_fill_manual(values = tissue_colors) +
    labs(x = "Sample (Patient_File_Tissue)", y = "Tissue") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
      axis.text.y = element_text(size = 10, face = "bold"),
      axis.title.y = element_text(size = 12, angle = 0, vjust = 0.5, face = "bold"),
      legend.position = "none",
      panel.grid = element_blank()
    )
  
  combined_plot <- p_cell / p_stage / p_tissue + plot_layout(heights = c(8, 0.53, 0.7))
  
  print(combined_plot)
    ggsave(paste0(filename_prefix, "_", group_var, ".pdf"),
         plot = combined_plot, 
         width = 16, height = 6)
  
    wide_data_cell <- dcast(sample_comp, patient + file + tissue ~ cell_type, 
                          value.var = "percentage", fill = 0)
  wide_data_stage <- dcast(stage_comp, patient + file + tissue ~ stage, 
                           value.var = "percentage", fill = 0)
  
#   fwrite(wide_data_cell, paste0(filename_prefix, "_cell_type_summary.csv"))
#   fwrite(wide_data_stage, paste0(filename_prefix, "_", group_var, "_summary.csv"))
  
  return(list(
    combined_plot = combined_plot,
    dt_subset = dt_subset,
    group_var = group_var
  ))
}

plot_polar_coord <- function(dt_subset, group_var, title_suffix) {
  stage_counts <- dt_subset[, .(count = .N), by = .(get(group_var), cell_type)]
  setnames(stage_counts, "get", "stage")
  
  stage_totals <- stage_counts[, .(total = sum(count)), by = stage]
  stage_counts <- merge(stage_counts, stage_totals, by = "stage")
  stage_counts[, ratio := count / total * 100]
  
  stage_levels <- c('Nor', 'Hyp', 'MiD', 'MoD', 'SD&CA', 'ICA', 'MCA')
  if (group_var == "ME") {
    stage_levels <- paste0(stage_levels, "-ME")
  }
  
  stage_levels <- stage_levels[stage_levels %in% unique(stage_counts$stage)]
  stage_counts[, stage := factor(stage, levels = stage_levels)]
  
  stage_counts[, cell_type := factor(cell_type, levels = cell_type_levels)]
  
  polar_cols <- cell_type_colors
  
  p_polar <- ggplot(stage_counts, aes(x = stage, y = ratio)) +
    geom_hline(
      aes(yintercept = y), 
      data.frame(y = c(0:4) * 25),
      color = "lightgrey"
    ) +
    geom_col(aes(fill = cell_type), size = 2, width = 0.7) +
    geom_segment(
      aes(
        x = stage,
        y = 0,
        xend = stage,
        yend = 100
      ),
      linetype = "dashed",
      color = "gray12"
    ) +
    scale_fill_manual(values = polar_cols, name = "Cell Type") +
    coord_polar(start = 0) +
    annotate(
      "text",
      x = rep(0.5, 5),
      y = c(0, 25, 50, 75, 100),
      label = c("0%", "25%", "50%", "75%", "100%"),
      size = 3
    ) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(color = "gray12", size = 12),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "bottom",
      panel.grid = element_blank()
    ) +
    labs(title = paste("Cell Type Composition in", title_suffix, "\n(Polar Coordinate)"))
  
  stage_counts_no_other <- stage_counts[cell_type != "Other"]
  
  if (nrow(stage_counts_no_other) > 0) {
    y_breaks <- seq(0, 100, by = 10)#seq(0, ceiling(max_ratio/10)*10, by = 10)
    
    p_polar_no_other <- ggplot(stage_counts_no_other, aes(x = stage, y = ratio)) +
      geom_hline(
        aes(yintercept = y), 
        data.frame(y = y_breaks),
        color = "lightgrey"
      ) +
      geom_col(aes(fill = cell_type), size = 2, width = 0.7) +
      geom_segment(
        aes(
          x = stage,
          y = 0,
          xend = stage,
          yend = max(y_breaks)
        ),
        linetype = "dashed",
        color = "gray12"
      ) +
      scale_fill_manual(values = polar_cols[names(polar_cols) != "Other"], name = "Cell Type") +
      coord_polar(start = 0) +
      annotate(
        "text",
        x = rep(0.5, length(y_breaks)),
        y = y_breaks,
        label = paste0(y_breaks, "%"),
        size = 3
      ) +
      theme_minimal() +
      theme(
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "gray12", size = 12),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "bottom",
        panel.grid = element_blank()
      ) +
      labs(title = paste("Cell Type Composition in", title_suffix, "\n(Without 'Other', Polar Coordinate)"))
  } else {
    p_polar_no_other <- NULL
  }
    print(p_polar)
    print(p_polar_no_other)
  return(list(
    polar_plot = p_polar,
    polar_plot_no_other = p_polar_no_other
  ))
}


dt <- mimer_median_enrich_development_ME

dt_me <- dt[ME != "unknown"]
plot_me <- analyze_and_plot(
  dt_me, 
  "Cell Type Composition in Microenvironment",
  "ME_cell_type_composition",
  "ME")

if (!is.null(plot_me$dt_subset)) {
  polar_plots_me <- plot_polar_coord(plot_me$dt_subset, "ME", "ME Regions")
  
  ggsave("ME_cell_type_polar.pdf", 
         plot = polar_plots_me$polar_plot, 
         width = 8, height = 8)
  
  if (!is.null(polar_plots_me$polar_plot_no_other)) {
    ggsave("ME_cell_type_polar_no_other.pdf", 
           plot = polar_plots_me$polar_plot_no_other, 
           width = 8, height = 8)
  }
}


# 
dt_dev <- dt[development != "unknown"]
plot_dev <- analyze_and_plot(
  dt_dev,
  "Cell Type Composition in Epi/Tumor",
  "development_cell_type_composition",
  "development")

if (!is.null(plot_dev$dt_subset)) {
  polar_plots_dev <- plot_polar_coord(plot_dev$dt_subset, "development", "Development Regions")
  ggsave("development_cell_type_polar.pdf",
         plot = polar_plots_dev$polar_plot,
         width = 8, height = 8)
  
  if (!is.null(polar_plots_dev$polar_plot_no_other)) {
    ggsave("development_cell_type_polar_no_other.pdf",
           plot = polar_plots_dev$polar_plot_no_other,
           width = 8, height = 8)
  }
}


#Extended  Data Fig8f -----
 p0 <- SpatialDimPlot(sp,group.by = 'new_cc_15',images = df[df$file=='KT1_2','image'],stroke = NA,pt.size.factor = 1.4,image.alpha = 0)+
    scale_fill_manual(values = c(color_mapping15,'other'='lightgrey'),name='')+
    theme_minimal() +
    coord_cartesian(xlim = common_xlim, ylim = common_xlim) +
    scale_x_continuous(breaks = common_xbreaks) +
    scale_y_continuous(breaks = common_xbreaks) +
    theme(
      aspect.ratio = 1,legend.position = 'top',legend.direction = 'horizontal')
  
  
  cols_mimer
  colnames(sp@meta.data)
  p1 <- SpatialDimPlot(sp,group.by = 'new_mimer',images = df[df$file=='KT1_2','image'],stroke = NA,pt.size.factor = 1.4,image.alpha = 0)+
    scale_fill_manual(values = cols_mimer,name='')+
    theme_minimal() +
    coord_cartesian(xlim = common_xlim, ylim = common_xlim) +
    scale_x_continuous(breaks = common_xbreaks) +
    scale_y_continuous(breaks = common_xbreaks) +
    theme(
      aspect.ratio = 1,legend.position = 'top',legend.direction = 'horizontal')
  p1


#Extended  Data Fig8g -----
load("KL.RData")
load("class.RData")
brain_cell_class_new
brain_sgraph_KL_mst_cons

require(purrr)
brain_sgraph_KL_mst_cons = reduce(brain_sgraph_KL_mst_cons, `+`) / length(brain_sgraph_KL_mst_cons)

library(tidyverse)
brain_cell_class_new<-brain_cell_class_new %>% reduce(inner_join, by = "id") 
brain_cell_class_new$freq.Freq<-rowMeans(brain_cell_class_new[,2:46])

long <- reshape2::melt(brain_sgraph_KL_mst_cons,id.vars= rownames(brain_sgraph_KL_mst_cons))
long
mst_cons_am <- brain_sgraph_KL_mst_cons
mst_cons_node <- data.frame(id=rownames(mst_cons_am), label=rownames(mst_cons_am))
directed = FALSE
if (!directed) mst_cons_am[upper.tri(mst_cons_am, diag = T)] <- NA

mst_cons_am <- data.frame(id=rownames(mst_cons_am), mst_cons_am, check.names=F)
mst_cons_edge <- reshape2::melt(mst_cons_am) %>% na.omit() %>% magrittr::set_colnames(c('from', 'to', 'value'))
mst_cons_edge

nodes <- data.frame(name = unique(union(mst_cons_edge$from, mst_cons_edge$to)))

nodes$number = brain_cell_class_new[brain_cell_class_new$id %in% nodes$name,"freq.Freq"]
nodes
class(nodes)

rownames(mst_cons_edge) <- 1:nrow(mst_cons_edge)
edges <- mst_cons_edge
colnames(edges) <- c("from","to","weighted")
class(edges)
edges
edges = edges[!edges$weighted == 0.00,]
library(tidygraph)
igraph::graph_from_data_frame(edges, vertices = nodes) %>% as_tbl_graph() -> g
gr1_layout2 <- create_layout(g, layout = "kk")
gr1_layout2
gr1_layout2[1,1:2] <- c(0.9,0.1)
gr1_layout2[2,1:2] <- c(-0.8,0.3)
gr1_layout2[3,1:2] <- c(0.15,0)
gr1_layout2[4,1:2] <- c(0.2,1)
gr1_layout2[5,1:2] <- c(-0.13,-0.7)

use_color=cell_color[cell_color$ct %in% gr1_layout2$name,"color"]

p<-ggraph(gr1_layout2) +
  scale_color_manual(values = use_color) +
  geom_edge_hive(aes(width = weighted),colour ="#6ea6cd")+
  geom_node_point(aes(size = number,colour = name))+ 
  scale_size(range = c(5,20),breaks=seq(50,450,100),limits = c(50,450))+
  geom_node_text(aes(label=name),size=3) +
  scale_edge_colour_gradientn(
    colours =  c("#4575b4","#f1f9d8","#f0663f")
  )+   guides(size= guide_legend(ncol = 2))+
  theme_graph() + expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5))


#Extended  Data Fig8h 8i -----
#leftz
ligand_target_matrix = readRDS("ligand_target_matrix.rds")
lr_network = readRDS("lr_network.rds")
weighted_networks = readRDS("weighted_networks.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% dplyr::inner_join(lr_network %>% dplyr::distinct(from,to), by = c("from","to"))
# filter out ligands not in ligand_target_matrix
ligands = lr_network$from %>% unique()
ligands = intersect(ligands, colnames(ligand_target_matrix))
receptors = lr_network$to %>% unique()
lr_network <- lr_network %>% filter(from %in% ligands & to %in% receptors) 
## receiver
xx<-subset(sc,tissue=='Tu');Idents(xx)<-xx$L4_C
receiver = "FB_C3_COL1A1"
expressed_genes_receiver = get_expressed_genes(receiver, xx, pct = 0.10)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
## sender
sender_celltypes = c("Mac_C2_SPP1")
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, xx, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
# DEG
my<-subset(xx,L3_C=='Fibroblast')
DE_table_receiver = FindMarkers(object = my, ident.1 = "FB_C3_COL1A1", min.pct = 0.10,only.pos = TRUE) %>% rownames_to_column("gene")
geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
# define potential ligands
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
# ligand activity
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities
best_upstream_ligands = ligand_activities %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
# activate target gene
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
save(vis_ligand_target,file = "Mac_CAF_NICHENET.RData")
new<-vis_ligand_target[c('SPP1','MMP9','IL1RN','TGFB1'),1:33]
new %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + 
  theme(axis.text = element_text(size=14,colour = "black")) + scale_fill_gradientn(colors =rev(getPalette(10)))

#right
## receiver
xx<-subset(sc,tissue=='Tu');Idents(xx)<-xx$L4_C
receiver = "Mac_C2_SPP1"
expressed_genes_receiver = get_expressed_genes(receiver, xx, pct = 0.10)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
## sender
sender_celltypes = c("FB_C3_COL1A1")
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, xx, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
# DEG
my<-subset(xx,L3_C=='Macrophage')
DE_table_receiver = FindMarkers(object = my, ident.1 = "Mac_C2_SPP1", min.pct = 0.10,only.pos = TRUE) %>% rownames_to_column("gene")
geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
# define potential ligands
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
# ligand activity
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities
best_upstream_ligands = ligand_activities %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
# activate target gene
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

new<-vis_ligand_target[c('TGFB3','COL5A3','COL1A1','CCL11','IGF2'),4:36]
new %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + 
  theme(axis.text = element_text(size=14,colour = "black")) + scale_fill_gradientn(colors =rev(getPalette(10)),breaks=c(0,0.002,0.004))


# E8j 8k
# 
library(ggplot2)
library(scales)
library(data.table)
library(dplyr)
library(tidyr)

analysis_type <- "development"

setwd('C:/Users/Ron/Desktop/Figs/E8/')
mimer_cor_development_ME <- fread('C:/Users/Ron/Desktop/Figs/E8/E8_ij_mimer_cor_development_ME.csv')
setDT(mimer_cor_development_ME)

if (analysis_type == "ME") {
  type_col <- "ME"
  filtered_types <- c("Nor-ME", "Hyp-ME", "MiD-ME", "MoD-ME", "SD&CA-ME", "ICA-ME")
} else {
  type_col <- "development"
  filtered_types <- c("Nor", "Hyp", "MiD", "MoD", "SD&CA", "ICA")
}

thresholds <- c("top5_perc_mimer", "top10_perc_mimer", "top20_perc_mimer")
threshold_names <- c("Top 5%", "Top 10%", "Top 20%")
threshold_colors <- c("Top 5%" = "#FF7F0E", "Top 10%" = "#3fa02c", "Top 20%" = "#1F77B4")

patient_stats_list <- list()

for (patient_id in unique(mimer_cor_development_ME$patient)) {
 patient_data <- mimer_cor_development_ME[patient == patient_id & get(type_col) %in% filtered_types]
  
  if (nrow(patient_data) > 0) {
    for (current_type in filtered_types) {
     type_data <- patient_data[get(type_col) == current_type]
      
      if (nrow(type_data) > 0) {
        for (i in seq_along(thresholds)) {
          threshold <- thresholds[i]
          threshold_name <- threshold_names[i]
          
          mimer_cells <- sum(type_data[[threshold]] == "MIMER", na.rm = TRUE)
          total_cells <- nrow(type_data)
          prop <- ifelse(total_cells > 0, mimer_cells / total_cells, 0)
          
          patient_stats_list[[paste(patient_id, current_type, threshold_name, sep = "_")]] <- 
            data.frame(
              patient = patient_id,
              type = current_type,
              threshold = threshold_name,
              prop = prop,
              n_cells = total_cells
            )
        }
      }
    }
  }
}

patient_stats <- rbindlist(patient_stats_list)

patient_stats$type <- factor(patient_stats$type, levels = filtered_types)
patient_stats$threshold <- factor(patient_stats$threshold, levels = threshold_names)

type_counts <- mimer_cor_development_ME[
  get(type_col) %in% filtered_types, 
  .(n_cells = .N), 
  by = .(patient, get(type_col))
]
setnames(type_counts, "get", "type")
type_counts$type <- factor(type_counts$type, levels = filtered_types)

type_patient_counts <- type_counts[, .(n_patients = .N), by = type]
type_patient_counts <- type_patient_counts[order(type)]

x_labels <- sapply(filtered_types, function(t) {
  n_patients <- type_patient_counts[type == t, n_patients]
  if (length(n_patients) > 0 && n_patients > 0) {
    return(paste0(t, "\n(n=", n_patients, ")"))
  } else {
    return(t)
  }
})

p <- ggplot(patient_stats, aes(x = type, y = prop, group = threshold, color = threshold)) +
 geom_line(aes(group = interaction(patient, threshold)), alpha = 0.3, linewidth = 0.5) +
  stat_summary(
    aes(group = threshold), 
    fun = mean, 
    geom = "line", 
    linewidth = 1.5
  ) +  stat_summary(
    aes(group = threshold), 
    fun = mean, 
    geom = "point", 
    size = 3
  ) + facet_wrap(~ patient, scales = 'free_y') + scale_color_manual(values = threshold_colors) +
    scale_x_discrete(limits = filtered_types, labels = x_labels) + scale_y_continuous(
    labels = percent_format(accuracy = 1)) +
  labs(
    title = ifelse(analysis_type == "ME", 
                   "MIMER+ Spot Proportion Across ME Types by Patient", 
                   "MIMER+ Spot Proportion Across Epi/Cancer Types by Patient"),
    x = ifelse(analysis_type == "ME", "ME Type", "Epi/Cancer Type"),
    y = "Proportion of MIMER+ Spots",
    color = "Threshold"
  ) + theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 10, face = "bold", margin = margin(r = 10)),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 10)),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_line(color = "gray95", linewidth = 0.3),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "gray90", color = "gray80")
  )

options(repr.plot.width=16,repr.plot.height=10)
print(p)
