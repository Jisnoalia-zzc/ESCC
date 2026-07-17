#Fig s1 a&b&c -----
library(dplyr)
library(purrr)
TCR_info <- sample_info_excel[,c(1, 55:68)]
TCR_info <- TCR_info[-1, ]
colnames(TCR_info) <- TCR_info[1, ]
TCR_info <- TCR_info[-1, ]
TCR_info <- TCR_info[1:16,]

TCR_info[is.na(TCR_info)] <- 0
TCR_info <- as.data.frame(TCR_info)
rownames(TCR_info) <- TCR_info[, 1]
TCR_info <- TCR_info[,-1]
rownames(TCR_info) <- paste0('p', rownames(TCR_info))

TCR_info_df <- as.data.frame(TCR_info)
TCR_info_df$sample <- rownames(TCR_info) 
TCR_long <- pivot_longer(TCR_info_df, 
                             cols = -sample,        
                             names_to = "patient",   
                             values_to = "scTCR")    
TCR_long$sampleID <- gsub("-", ".", TCR_long$sampleID)

TCR_long <- TCR_long[3:4]

ATAC_info <- sample_info_excel[,c(1, 42:54)]
ATAC_info <- ATAC_info[-1, ]
colnames(ATAC_info) <- ATAC_info[1, ]
ATAC_info <- ATAC_info[-1, ]
ATAC_info <- ATAC_info[1:16,]

ATAC_info <- as.data.frame(ATAC_info)
rownames(ATAC_info) <- ATAC_info[, 1]
ATAC_info <- ATAC_info[,-1]
ATAC_info[is.na(ATAC_info)] <- 0

rownames(ATAC_info) <- paste0('p', rownames(ATAC_info))

ATAC_info_df <- as.data.frame(ATAC_info)
ATAC_info_df$sample <- rownames(ATAC_info) 
ATAC_long <- pivot_longer(ATAC_info_df, 
                         cols = -sample,        
                         names_to = "patient",   
                         values_to = "scATAC")    
ATAC_long$sampleID <- paste0(ATAC_long$sample,'_', ATAC_long$patient)
ATAC_long <- ATAC_long[3:4]

scRNA_long <- plot_df_use[1:2]

scRNA_data <- subset(scRNA_long, scRNA == 1)
scATAC_data <- subset(ATAC_long, scATAC == 1)
spatial_data <- subset(spatial_long, spatial == 1)
scTCR_data <- subset(TCR_long, scTCR == 1)

scATAC_data$scATAC <- as.numeric(scATAC_data$scATAC)
spatial_data$spatial <- as.numeric(spatial_data$spatial)
scTCR_data$scTCR <- as.numeric(scTCR_data$scTCR)

sort(unique(c(scRNA_data$sampleID, scATAC_data$sampleID, spatial_data$sampleID, scTCR_data$sampleID)))

merged_df <- list(scRNA_data, scATAC_data, spatial_data, scTCR_data) %>%
  reduce(full_join, by = "sampleID")
merged_df <- as.data.frame(merged_df)
merged_df[is.na(merged_df)] <- 0

rownames(merged_df) <- merged_df$sampleID
merged_df <- merged_df %>% separate(sampleID, into = c("patient", "tissue"), sep = "_")
merged_df$sampleID <- rownames(merged_df)
merged_df$rowsum <- rowSums(sapply(merged_df[, 3:6], as.numeric))
saveRDS(merged_df, file = "merged_df.rds")
merged_df <- merged_df[, c("sampleID", "scRNA", "scTCR", "scATAC", "spatial", "patient", "tissue", "rowsum")]
merged_df[, 2:5] <- lapply(merged_df[, 2:5], as.numeric)

merged_df$patient <- factor(merged_df$patient, levels = c("p1022", "p0901", "p1019", "p1104",
                                                          "p1119", "p1123", "p1201", "p1204",
                                                          "p1209", "p1210", "p1229", "p0118",
                                                          "p0128", "p0308", "p0330", "p0826"))
merged_df <- merged_df[order(merged_df$patient), ]


tt = merged_df[2:5]
tt$scTCR = as.numeric(gsub(1,2,tt$scTCR))
tt$scATAC = as.numeric(gsub(1,3,tt$scATAC))
tt$spatial = as.numeric(gsub(1,4,tt$spatial))
rownames(tt)  = merged_df$sampleID

anno = merged_df[,c(7,6)]
anno2 = left_join(anno,sample_info,by="patient")
rownames(anno2) = rownames(anno)
anno2[anno2=='-'] = NA
# rownames(anno) = plot_df_use$sampleID
cols_20 <- c("#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", 
             "#59A14F", "#8CD17D", "#B6992D", "#F1CE63", "#499894", 
             "#86BCB6", "#E15759", "#FF9D9A", "#79706E", "#BAB0AC", 
             "#D37295", "#FABFD2", "#B07AA1", "#D4A6C8", "#9D7660", 
             "#D7B5A6")
tissue_col = cols_20[1:length(unique(anno$tissue))]
names(tissue_col) = sort( unique(anno$tissue))
patient_col = cols_20[1:length(unique(anno$patient))]
names(patient_col) = unique(anno$patient)

stage_TNM = sort(unique(anno2$stage_TNM)[1:6])
cols_stage  = RColorBrewer::brewer.pal(6,"YlOrRd")
names(cols_stage) = stage_TNM
ann_colors=list(
  patient = patient_col,
  tissue=tissue_col,
  sex = c("F" = "#add18e", "M" ="#33c0ed"),
  pathology = c("squamous cell carcinoma" = "#1f639a",'precancerous lesions' = "#d38e91"),
  # age_group = c("0-45"="#d6e1c2","45-60"="#b6c9b8",">60" ="#9db7a5"),
  stage_TNM = cols_stage
)
anno2 = anno2[,c(5,3,4,1,2)]
color_list <- c("#f2f2f2", "#e7ca93", "#45a47a", "#e7959c", "#bacce9")  # 每个离散值对应一种颜色
breaks_list <- c(-1, 0, 1, 2, 3, 4)


pheatmap::pheatmap(t(tt),cluster_cols = F,cluster_rows = F,
                   cellwidth = 5,cellheight = 40,
                   legend = FALSE,
                   annotation_col = anno2,
                   # color = c("white","#bacce9"),
                   annotation_colors = ann_colors,
                   color = color_list, 
                   breaks = breaks_list, 
                   show_colnames = F,
                   filename = "~/meta_new_board.pdf",
                   width = 20,height = 20,
                   border_color = 'black'
)

dev.off()

library(ggvenn)
tt$sample = str_split(rownames(tt),'\\_',simplify = T)[,1]

# tt_sub = subset(tt,sample==sa)
scrna = ifelse(tt$scRNA==1,rownames(tt),0)
sctcr = ifelse(tt$scTCR==2,rownames(tt),0)
scatac =  ifelse(tt$scATAC==3,rownames(tt),0)
spatial =  ifelse(tt$spatial==4,rownames(tt),0)

venn_list = list(
  scRNA = scrna[scrna!=0],
  scTCR = sctcr[sctcr!=0],
  scATAC = scatac[scatac!=0],
  spatial = spatial[spatial!=0]
)

result_list = lapply(names(venn_list), function(pro){
  cc = unique(str_split( venn_list[[pro]],"\\_",simplify = T)[,1])
  return(cc)
})
names(result_list) = names(venn_list)
ggvenn(
  data = result_list,columns = NULL,show_elements = F,label_sep = "\n",show_percentage = T, digits = 1,
  fill_color = c("#e7ca93", "#45a47a", "#e7959c", "#bacce9"), fill_alpha = 0.7, stroke_color = "white",stroke_alpha = 0.7, stroke_size = 0.5,
  stroke_linetype = "solid", set_name_color = "black", set_name_size = 6,  text_color = "black", text_size = 4 )+
  theme(plot.title = element_text(hjust = 0.5,size = 20))
ggsave(paste("~/2-patient_level_ggvenn_new.pdf"),width = 6,height = 6)


ggvenn(data = venn_list, columns = NULL, show_elements = F, label_sep = "\n", show_percentage = T, digits = 1,
fill_color = c("#e7ca93", "#45a47a", "#e7959c", "#bacce9"), fill_alpha = 0.7, stroke_color = "white", stroke_alpha = 0.7,
  stroke_size = 0.5, stroke_linetype = "solid", set_name_color = "black",set_name_size = 6, text_color = "black", text_size = 4)+
  theme(plot.title = element_text(hjust = 0.5,size = 20))
ggsave(paste("~/2-sample_level_ggvenn_new.pdf"),width = 6,height = 6)


#Fig s1 d -----
load("/workspace/wangxiliang/esca/sc_rna/scRNA_t.RData")
my$L4_C<-Idents(my)
my$L4_C<-factor(my$L4_C,levels = c('T_C2_MKI67','T_C1_HSPA1A','CD8','CD4','NK_C2_XCL1','NK_C1_CD16'))
T_plotG <- c('SPON2','FCGR3A','CX3CR1','CD160','AREG','XCL1','FCER1G','XCL2',
             'MAL','CD4','ICOS','IL6ST','CD8A','CCL4L2','CD8B','GZMK','HSPA1B','HSPA1A','DNAJB1','JUN','STMN1',
             'MKI67','RRM2','TOP2A')
celltype <- c('NK_C1_CD16',  'NK_C2_XCL1','CD4','CD8', 'T_C1_HSPA1A',  'T_C2_MKI67')
ctidx <- c(4,4,4,4,4,4)
getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))
dotplot(my,T_plotG,celltype,ctidx,color = rev(getPalette(10)))


#Fig s1 e -----
FetchData(my,c("UMAP_1",'UMAP_2','subC'))  %>% 
  ggplot() +
  ggrastr::geom_point_rast(aes(UMAP_1,UMAP_2,color = subC),size =.4/.pt) +
  scale_color_manual(values =c('#EC593B','#F7F499','#1E90FF','#20B2AA','#DF9A89','#B0DD8A')) +
  theme_classic() +
  coord_fixed(ratio = 1) + 
  theme(panel.background = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 7.5),
        axis.title =  element_text(size = 7.5),
        axis.text = element_text(size = 7.5)) + 
  guides(colour = guide_legend(ncol = 1,
                               override.aes=list(shape=19, size=5, linetype=0)))


#Fig s1 f -----
T_plotG <- c('SPON2','FCGR3A','CX3CR1','CD160','AREG','XCL1','FCER1G','XCL2',
             'MAL','CD4','ICOS','IL6ST','CD8A','CCL4L2','CD8B','GZMK','HSPA1B','HSPA1A','DNAJB1','JUN','STMN1',
             'MKI67','RRM2','TOP2A')
plot_df <- read.csv("/home/users/zhangzhichao/workspace/project/ESCC/NKT_marker_gene_mean_z.csv",header = T)
plot_df$subC = factor(plot_df$subC,
                       levels = c("NK_C1_CD16","NK_C2_XCL1","CD4","CD8","T_C2_MKI67","T_C1_HSPA1A"))
plot_df$variable = factor(plot_df$variable ,
                          levels = T_plotG)

celltype <- c('NK_C1_CD16',  'NK_C2_XCL1','CD4','CD8', 'T_C1_HSPA1A',  'T_C2_MKI67')
ctidx <- c(4,4,4,4,4,4)
getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))
dotplot(plot_df,T_plotG,celltype,ctidx,color = rev(getPalette(10)))


#Fig s1 g -----
Idents(sc)<-sc$tissue
t<-as.matrix(table(Idents(sc),sc@meta.data$L1_C))
t<-as.data.frame(t)
total<-as.data.frame(table(Idents(sc)))
colnames(total)<-c('Var1','total')
t<-merge(t,total,all = TRUE)
colnames(t)<-c('Tissue','Cell_Subcluster','Cell_Number','Total')
t$ratio<-t$Cell_Number/t$Total
t$Tissue<-factor(t$Tissue,levels = c('PBMC','nLN','pLN','Nor','Adj','Tu'))
ggplot(t, aes(Tissue,ratio, fill=Cell_Subcluster)) +geom_col(position = 'stack', width = 0.6)+
  theme(axis.text.x = element_text(size=14,angle=60,hjust=1,colour = "black"),panel.background = element_blank(),
        axis.text.y = element_text(size=14,colour = "black"),panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.7),legend.text = element_text(size = 14),legend.title = element_blank(),
        legend.background = element_blank(),legend.key = element_blank(),title = element_text(size = 14))+
  xlab('Tissue')+ylab('Proportion')+scale_fill_manual(values=pal_igv("default")(20))

Idents(sc)<-sc$patient
t<-as.matrix(table(Idents(sc),sc@meta.data$L1_C))
t<-as.data.frame(t)
total<-as.data.frame(table(Idents(sc)))
colnames(total)<-c('Var1','total')
t<-merge(t,total,all = TRUE)
colnames(t)<-c('Patient','Cell_Subcluster','Cell_Number','Total')
t$ratio<-t$Cell_Number/t$Total
ggplot(t, aes(Patient,ratio, fill=Cell_Subcluster)) +geom_col(position = 'stack', width = 0.6)+
  theme(axis.text.x = element_text(size=14,angle=60,hjust=1,colour = "black"),panel.background = element_blank(),
        axis.text.y = element_text(size=14,colour = "black"),panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.7),legend.text = element_text(size = 14),legend.title = element_blank(),
        legend.background = element_blank(),legend.key = element_blank(),title = element_text(size = 14))+
  xlab('Patient')+ylab('Proportion')+scale_fill_manual(values =pal_igv("default")(20))
