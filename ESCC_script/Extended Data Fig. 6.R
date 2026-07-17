#---------------------------------------
#--------------Extended data.6a
#----------------------------------------

library(tidyverse)
library(RColorBrewer)
library(tidyr)
library(rstatix)
library(pheatmap)
library(Seurat)
library(gridExtra)

sp <- readRDS("/lustre/home/xhzh/scSpatial/xiaodu/codes/zxh/ALL_HALLMARK.rds")
meta <- read.table("/lustre/home/xhzh/scSpatial/xiaodu/codes/zxh/HALLMARK_custompath_TLS_meta.txt", sep = "\t", header=T)
tls_order <- readLines("/lustre/home/xhzh/scSpatial/xiaodu/codes/zxh/TLS/tls_order.txt")

meta_df <- meta[meta$tissue=="AF" | meta$tissue=="ANF" | meta$tissue=="NF" | meta$tissue=="TF", ]
TAN_subset <- subset(sp, subset=tissue %in% c("AF", "ANF", "NF", "TF"))

tls <- data.frame(key = tls_order, value = c(rep("TLS_G1", times=7), rep("TLS_G2", times=42), rep("TLS_G3", times=31)), value1 = c(rep("Peri_G1", times=7), rep("Peri_G2", times=42), rep("Peri_G3", times=31)))

load("/lustre/home/xhzh/scSpatial/xiaodu/codes/zxh/TLS/rctd_full_matrix.Rdata")
df <- as.data.frame(rctd_full_matrix)
df$celltype <- rownames(df)
max_rows <- apply(df, 2, function(col) rownames(df)[which.max(col)])
max_values <- apply(df, 2, max)
epi_cancer_spots <- names(max_values)[max_values > 0.5 & (max_rows == "cancer" | max_rows == "Epithelium")]

cytotoxic = c("GNLY", "IFNG", "NKG7", "PRF1", "GZMA", "GZMB", "GZMH", "GZMK")
icb = c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "BTLA", "CD274", "CD276")
TLS_Maturation_antitumor = c("STAT1", "ISG15", "IFI27", "IFI30", "MX1", "IFITM3", "GBP1", "BST2", "TAP1", "CXCL9", "CXCL10", "CXCL11", "CXCL13", "CX3CL1", "CCL28", "C7", "CFB1", "HLA-C")
Antiviral_mucosal_immune_activity = c("MUC5B", "BPIFB1", "BPIFB2", "TFF3", "DEFB1", "SCGB3A1", "PIGR", "JCHAIN", "WFDC2", "SAA1", "AZGP1", "LCN2", "PRR4", "ZG16B", "AGR2", "LTF", "SLPI", "SAA2", "CRISP3", "CXCL17")
Chemokine = c("CCL21","CCL19","CXCL12","CXCL14","CX3CL1","CCL17","CCL3","CXCL8","CCL2","XCL1","XCL2")
Chemokine_receptors = c("CCR2","CCR6","CCR5","CXCR4","CCR7","XCR1","CCR6","CCR10","CXCR3","CXCR4","CCR1","CCR3")
TLRs_and_Recptors = c("MYD88","TICAM1","TLR1","TLR2","TLR3","TLR4","TLR5","TLR6","TLR7","TLR8","TLR9","TLR11")
TLS_B_maturation_slgA_transport = c("IGHA1", "IGHD", "IGHM1", "IGLC1", "IGLC3", "IGLC71", "JCHAIN", "PIGR")

marker_list = list(cytotoxic, icb, TLS_Maturation_antitumor, Antiviral_mucosal_immune_activity, Chemokine, Chemokine_receptors, TLRs_and_Recptors, TLS_B_maturation_slgA_transport)
names(marker_list) = c("Cytotoxic gene", "ICB gene", "TLS Maturation antitumor", "Antiviral mucosal immune activity", "Chemokine", "Chemokine Receptors", "TLRs and Receptors", "TLS B maturation slgA transport")

TAN_subset@meta.data$tls_group_noepicaner <- "NA"
for (tls_label in tls_order) {
  core_cells <- meta_df[meta_df$TLS_cluster_label == tls_label, ]
  peripheral_cells <- meta_df[meta_df$TLS_cluster_label == "Non-TLS" & meta_df$peripheral_label_formatted == tls_label,]
  TAN_subset@meta.data$tls_group_noepicaner[TAN_subset@meta.data$cellID %in% rownames(core_cells) & !TAN_subset@meta.data$cellID %in% epi_cancer_spots] <- tls[tls$key==tls_label,]$value
  TAN_subset@meta.data$tls_group_noepicaner[TAN_subset@meta.data$cellID %in% rownames(peripheral_cells) & !TAN_subset@meta.data$cellID %in% epi_cancer_spots] <- tls[tls$key==tls_label,]$value1
}

table(TAN_subset$tls_group_noepicaner)
meta_df1 <- TAN_subset@meta.data

available_genes <- rownames(GetAssayData(TAN_subset, slot = "data", assay="SCT"))
plots <- vector("list", 8)
plots1 <- vector("list", 8)

for (i in c(1:8)) {
  genelist <- marker_list[[i]]
  matched_genes <- intersect(genelist, available_genes)
  missing_genes <- setdiff(genelist, available_genes)
  gene_lab <- names(marker_list)[i]
  print(paste0(gene_lab, ":", missing_genes))
  expr_cy <- GetAssayData(TAN_subset, slot = "data", assay="SCT")[matched_genes, ]
  exp <- data.frame()
  exp1 <- data.frame()
  for (group in c("TLS_G1","TLS_G2", "TLS_G3") ) {
    cells <- meta_df1[meta_df1$tls_group_noepicaner == group, ]
    for ( g in matched_genes ) {
      exp[group, g] <- mean(expr_cy[g,rownames(cells)], na.rm = TRUE)
    }
  }
  texp <- t(exp)
  print(gene_lab)
  p <- pheatmap(texp, border_color = "black",angle_col = "45",main=gene_lab, legend = TRUE, cluster_rows = F, cluster_cols = F, fontsize_row = 12, fontsize_col = 12, silent = TRUE)
  gt <- p$gtable
  plots[[i]] <- gt

  for (group in c("Peri_G1","Peri_G2", "Peri_G3") ) {
    cells <- meta_df1[meta_df1$tls_group_noepicaner == group, ]
    for ( g in matched_genes ) {
      exp1[group, g] <- mean(expr_cy[g,rownames(cells)], na.rm = TRUE)
    }
  }
  texp1 <- t(exp1)
  p1 <- pheatmap(texp1, border_color = "black",angle_col = "45",main=gene_lab, legend = TRUE, cluster_rows = F, cluster_cols = F, fontsize_row = 12, fontsize_col = 12, silent = TRUE)
  gt1 <- p1$gtable
  plots1[[i]] <- gt1
}
Cairo::CairoPDF(file = "geneset_exp_TLS_Peri.pdf", width = 8, height = 8) 
layout <- grid.arrange(plots[[1]], plots[[2]], ncol = 2)
layout <- grid.arrange(plots1[[1]], plots1[[2]], ncol = 2)
layout <- grid.arrange(plots[[3]], plots[[4]], ncol = 2)
layout <- grid.arrange(plots1[[3]], plots1[[4]], ncol = 2)
layout <- grid.arrange(plots[[5]], plots[[6]], ncol = 2)
layout <- grid.arrange(plots1[[5]], plots1[[6]], ncol = 2)
layout <- grid.arrange(plots[[7]], plots[[8]], ncol = 2)
layout <- grid.arrange(plots1[[7]], plots1[[8]], ncol = 2)
dev.off()

#### E6B
library(tidyverse)
library(ggpubr)
library(data.table)
library(ggsignif)

setwd('C:/Users/Ron/Desktop/Figs/E6/')
combined_data <- fread('C:/Users/Ron/Desktop/Figs/E6/E6_b_combined_data_celltype_boxplot.csv')

cluster_colors <- c('#E41A1C', '#377EB8', '#4DAF4A') 
options(repr.plot.width = 16, repr.plot.height = 12)

core_data <- combined_data %>% filter(Region == "Core")
peri_data <- combined_data %>% filter(Region == "Peri")

comparisons <- list(
  c("TLS_G1", "TLS_G2"),
  c("TLS_G1", "TLS_G3"),
  c("TLS_G2", "TLS_G3")
)


create_abundance_plot <- function(data, region_name) {
 
  y_positions <- data %>%
    group_by(CellType) %>%
    summarise(
      max_val = max(Abundance, na.rm = TRUE),
      q3 = quantile(Abundance, 0.75, na.rm = TRUE)
    ) %>%
    mutate(
      base_y = max_val * 1.05,step1 = base_y * 1.10, step2 = base_y * 1.20,step3 = base_y * 1.30  )
  
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
  
  valid_celltypes <- data %>%
    group_by(CellType) %>%
    summarise(groups = n_distinct(TLS_cluster)) %>%
    filter(groups == 3) %>%
    pull(CellType)
  
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

core_plot <- create_abundance_plot(core_data, "Core")
peri_plot <- create_abundance_plot(peri_data, "Peri")

print(core_plot)
print(peri_plot)

# ggsave("TLS_core_celltype_abundance.pdf", core_plot,
#        width = 12, height = 8)
# ggsave("TLS_peri_celltype_abundance.pdf", peri_plot,
#        width = 12, height = 12)


#--------------------------------
#-----Extended_data6.e
#----------------------------------

#------

db_spt<-fread(file = './VJC_cdr3_all_TRUST4.csv')  %>%  as.data.frame() %>%
  mutate(cell_id = gsub("\\.", "_", cell_id)) %>% filter(.,Omic == 'SPT')

db_spt$tissue <- sub("_SPT$", "", db_spt$Tissue)
db_spt <- db_spt[,c('cell_id','tissue','library','Patient','vjc_cdr3')]

db_spt$tissue[db_spt$library == 'A_JNP'] <- 'A'

library(dplyr)
library(readr)
library(stringr)

setwd("./spt_barcode_tissue_type_cor_by_ZWM")

file_list <- list.files(pattern = "\\.csv$", full.names = TRUE)

result_list <- list()

for (file in file_list) {
  filename <- basename(file) %>% str_remove("\\.csv$")
  parts <- unlist(strsplit(filename, "-"))
  
  sample_id <- parts[1]
  tissue_type <- paste(parts[-1], collapse = "-")
  
  df <- read_csv(file, show_col_types = FALSE) %>%
    rename(Barcode = 1) %>% mutate(
      cellID = paste(sample_id, Barcode, sep = "_"),  tissue = tissue_type,library = sample_id ) %>%
    select(cellID, tissue)
  
  result_list[[file]] <- df
}

final_df <- bind_rows(result_list)

library(dplyr)
merged_db <- left_join(
  db_spt, 
  final_df, 
  by = c("cell_id" = "cellID"))

db_spt_new <- merged_db %>%
  mutate(
    tissue = ifelse(
      !is.na(tissue.y),tissue.y,tissue.x )
  ) %>% select(-tissue.x, -tissue.y)

metadata <- fread(file = './Metadata_only_TAN_sample_TLS_mark_df_labeled.csv')

options(repr.plot.width = 12,repr.plot.height = 8)

selected_spots<-metadata$cellID[which(metadata$TLS_mark==1|!metadata$peripheral_label_formatted == '')]

for (locus in c("IGK", 'IGL', 'IGH')) {
  
  vdata <- db_spt_new[grep(locus, db_spt_new$vjc_cdr3), ]
  vdata <- vdata[vdata$cell_id %in% selected_spots, ]
  
  tissue_T <- unique(vdata$vjc_cdr3[vdata$tissue == "T"])
  tissue_N <- unique(vdata$vjc_cdr3[vdata$tissue == "N"])
  tissue_A <- unique(vdata$vjc_cdr3[vdata$tissue == "A"])
  
  tissue_list <- list(
    T = tissue_T,
    N = tissue_N,
    A = tissue_A
  )
  fit <- eulerr::euler(tissue_list)
  g<- plot(fit,
           quantities = TRUE,
           fills = c("#E69F00", "#56B4E9", "#009E73"),
           main = paste0('Spatial ', locus," Venn Diagram")
  )
  pdf(file=paste0("D:/Postdoc/Projects/Single_cell/ZWM/SPT_", locus, "_tissue_venn.pdf"),4,4)
  print(g)
  dev.off()
}

#------2.

data_scBCR = fread(file = './Add_huaxi_scBCR_db.new.tsv')
data_scBCR$cell_id<-sub("_contig_[^_]*$", "", data_scBCR$sequence_id) %>% sub("ScBCR","ScRNA",.)

cell_type_mapping_table<-data.table::fread('./meta_data.csv')

cell_type_mapping_table_base<-cell_type_mapping_table[cell_type_mapping_table$Project %in% c('Our'),c('orig.ident','Celltype_L1','Celltype_L3_add_cnv','Tissue')] %>% as.data.frame()
cell_type_mapping_table_base$cell_id=gsub("_(?!.*_)", ".", cell_type_mapping_table_base$orig.ident,perl=TRUE)

scBCR_all<-merge(data_scBCR,cell_type_mapping_table_base,by="cell_id")

scBCR_all$vjc_cdr3 <- paste(
  scBCR_all$v_call_10x,
  scBCR_all$j_call_10x,
  scBCR_all$c_call,
  scBCR_all$junction_10x,
  sep = "_"
)
scBCR_Bcells <- scBCR_all %>% filter(.,Celltype_L1=='B/Antibody-Secreting')
scBCR_Bcells <- scBCR_Bcells[,c('Tissue','vjc_cdr3','locus')] %>% as.data.frame()

library(VennDiagram)

for (locus in c("IGK", 'IGL', 'IGH')) {
 vdata <- scBCR_Bcells[grep(locus, scBCR_Bcells$locus), ]
  
  tissue_T <- unique(vdata$vjc_cdr3[vdata$Tissue == "T"])
  tissue_N <- unique(vdata$vjc_cdr3[vdata$Tissue == "N"])
  tissue_A <- unique(vdata$vjc_cdr3[vdata$Tissue == "A"])
  
  tissue_list <- list(
    T = tissue_T,
    N = tissue_N,
    A = tissue_A
  )
  fit <- eulerr::euler(tissue_list)
  g<- plot(fit, 
           quantities = TRUE, 
           fills = c("#E69F00", "#56B4E9", "#009E73"),
           main = paste0('ScBCR', locus," Venn Diagram")
  )
  pdf(file=paste0("./Singlecell_", locus, "_tissue_venn.pdf"),4,4)
  print(g)
  dev.off()
}