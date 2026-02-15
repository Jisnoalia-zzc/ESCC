#Extended Data Fig 3 -----
#Extended Data Fig 3b -----
library(Seurat)
codex = qread(file = '~/workspace/proj_ESCC_STW_ZWM_2022_01/zzc/final_analysis/final_codex.qs')
codex$Celltype.y = factor(codex$Celltype.y,
                          levels = c("Tcells","Bcells","Myeloids",
                                     "Endothelium","Fibroblast","Epithelium"))

codex@meta.data %>% dplyr::select(Celltype.y,subCelltype) %>%
  distinct(subCelltype,.keep_all = T) %>%
  arrange(Celltype.y) -> rank
rank$subCelltype = gsub("_","-",
                        rank$subCelltype)
codex$subCelltype = gsub('_','-',
                         codex$subCelltype)
codex$subCelltype = factor(codex$subCelltype,
                           levels = rank$subCelltype)
qsave(codex,file = '~/workspace/proj_ESCC_STW_ZWM_2022_01/liuliqiu//final_analysis/final_codex_20260215.qs')
DimPlot(codex,group.by = 'subCelltype',raster = T)
#Extended Data Fig 3c -----
Idents(codex) = 'subCelltype'

codex <- NormalizeData(object = codex, normalization.method = "CLR", margin = 2)
codex  = ScaleData(codex)
avg = AverageExpression(codex,layer='data')
avg = avg$Akoya
rank = rank %>% filter(subCelltype!='Epithelium')
avg = avg[,match(rank$subCelltype,colnames(avg))]

palette_length <- 100
my_color <- colorRampPalette(c("#2696f2","#58a3e8", "white","#ff8f6b","#d9100b"))(palette_length)
my_breaks <- c(seq(-4, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05,4, length.out=floor(palette_length/2)))
anno_cols = rank
names(major) = unique(anno_cols$Celltype.y)
annotation_colors = list(
  subCelltype = annotation_colors$subCelltype[names(annotation_colors$subCelltype) !="Epithelium"],
  Celltype.y = annotation_colors$Celltype.y
  
)
rownames(anno_cols)=anno_cols$subCelltype
anno_cols = anno_cols[,c(2,1)]
save(cell_color,major,file = '~/ESCC_codex/codex_color.Rdata')
pdf("~/ESCC_codex/CODEX_Average_Scaled_Celltype_Marker_Expression.pdf",width = 10,height = 15)
pheatmap::pheatmap(

  t(scale(t(avg))),
  scale = 'none',
  cluster_cols = F,
  color = my_color,
  annotation_col = anno_cols,
  annotation_colors = annotation_colors,
  # display_numbers = T,
  breaks = my_breaks,
  cellwidth = 8,
  cellheight = 8,
  fontsize_row = 9,
  fontsize_col = 9,
  main = "CODEX Average Scaled Celltype Marker Expression",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean")
dev.off()
#Extended Data Fig 3d -----
load("/BGFS1/home/zhangzc/ESCC_codex/cor/scRNA_avg.Rdata")
p2r = openxlsx::read.xlsx('/BGFS1/home/zhangzc/ESCC_codex/cor/rna_protein_df_edit.xlsx')

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
scRNA_avg
scRNA_avg %>%
  mutate(RNA = rownames(.)) %>%
  left_join(p2r,by = "RNA") %>%
  filter(RNA !='KRT18') %>%
  column_to_rownames("Protein") %>%
  dplyr::select(-RNA) -> scRNA_avg2
all(rownames(scRNA_avg2)==rownames(codex_avg))
setdiff(rownames(scRNA_avg2),rownames(codex_avg))
setdiff(rownames(codex_avg),rownames(scRNA_avg2))
scRNA_avg2 = scRNA_avg2[match(rownames(codex_avg),rownames(scRNA_avg2)),]
scRNA_avg3 = t(scale(t(scRNA_avg2)))
codex_avg2 = t(scale(t(codex_avg)))
num_cols_scRNA <- ncol(scRNA_avg3)
num_cols_codex <- ncol(codex_avg2)
cross_correlation_matrix <- matrix(NA, nrow = num_cols_scRNA, ncol = num_cols_codex)
rownames(cross_correlation_matrix) <- colnames(scRNA_avg3)
colnames(cross_correlation_matrix) <- colnames(codex_avg2)

for (i in 1:num_cols_scRNA) {
  for (j in 1:num_cols_codex) {
    # 计算 matrix_A 的第 i 列和 matrix_B 的第 j 列之间的相关性
    cross_correlation_matrix[i, j] <- cor(scRNA_avg3[, i], codex_avg2[, j],method = 'pearson')
  }
}

scRNA_cols = rep('grey',num_cols_scRNA)
names(scRNA_cols) =  colnames(scRNA_avg3)
scRNA_cols[names(scRNA_cols) %in% c("CD4_C7_OX40","CD8_C6_CD39",'FB_C3_COL1A1','Mac_C2_SPP1','Endo_C3_RGCC')] = 'brown'
codex_cols = rep('grey',num_cols_codex)
Mimer = c("SPP1+Mac",'ISG15+Mac','Ki-67+Mac',
          'CD4T-aTreg','Tcells.cc',
          'PDCD1+Tex','PDCD1+Tex2','LAG3+Tex',
          'COL1A1+CAF','POSTN+CAF',
          'PLVAP+Endo','Ki-67+Endo')
names(codex_cols) =  colnames(codex_avg2)
codex_cols[names(codex_cols) %in% Mimer] = 'brown'
anno_cor_colors = list(
  codex_mimer = codex_cols,
  scRNA_mimer = scRNA_cols
)
anno_cols_codex = data.frame(
  row.names = names(codex_cols) ,
  codex_mimer = names(codex_cols)
)
anno_cols_scRNA = data.frame(
  row.names = names(scRNA_cols) ,
  scRNA_mimer = names(scRNA_cols)
)

palette_length <- 100
my_color <- colorRampPalette(c("#2696f2","#58a3e8", "white","#ff8f6b","#d9100b"))(palette_length)
my_breaks <- c(seq(-1, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05,1, length.out=floor(palette_length/2)))

pdf('scRNA_codex_cor_pheatmap.pdf',width = 8,height = 10 )
pheatmap::pheatmap(cross_correlation_matrix,scale = 'none',
                   clustering_distance_rows = "euclidean",
                   clustering_distance_cols = "euclidean",
                   color = my_color,
                   annotation_col = anno_cols_codex,
                   annotation_row = anno_cols_scRNA,
                   annotation_colors = anno_cor_colors,
                   breaks = my_breaks,
                   fontsize_row = 6,fontsize_col = 6,
                   cellwidth  = 6,cellheight  = 6)
dev.off()
#Extended Data Fig 3e -----
