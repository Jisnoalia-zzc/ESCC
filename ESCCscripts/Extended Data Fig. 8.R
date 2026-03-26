#Extended  Data Fig8 -----
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

# E8C1
library(ggplot2)
library(dplyr)
library(data.table)

setwd('C:/Users/Ron/Desktop/Figs/E8/')
file_stats <- fread('C:/Users/Ron/Desktop/Figs/E8/E8_c1_file_stats.csv')
# 定义组织类型颜色映射
tissue_colors <- c(
  "TF" = "#8B0000",
  "LpF" = "#CD5C5C", 
  "AF" = "#034272",
  "NF" = "#0d6cc0",
  "LnF" = "#87CEEB"
)

# 过滤实际存在的组织类型
present_tissues <- unique(file_stats$tissue)
available_colors <- tissue_colors[names(tissue_colors) %in% present_tissues]

# 处理p=0的情况
file_stats <- file_stats %>%
  mutate(
    p_adj = ifelse(p_value == 0, .Machine$double.xmin, p_value),
    neg_log_p = -log10(p_adj)
  )

# 计算y轴范围
finite_p <- file_stats$neg_log_p[is.finite(file_stats$neg_log_p)]
y_max <- if(length(finite_p) > 0) max(finite_p, na.rm = TRUE) else 10
y_max <- if(is.infinite(y_max) | is.na(y_max)) 10 else y_max
y_limit <- y_max * 1.2

# 处理Inf值
file_stats <- file_stats %>%
  mutate(neg_log_p_plot = ifelse(is.infinite(neg_log_p), y_limit * 0.9, neg_log_p))

# 创建火山图
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

# 添加显著标签（可选）
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
# 2. 整合坐标
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

# 定义颜色映射
tissue_colors <- c(
"TF" = "#8B0000",       # 深红色
"LpF" = "#CD5C5C",      # 中等红色
"AF" = "#034272",       # 深蓝色
"NF" = "#0d6cc0",       # 蓝色
"LnF" = "#87CEEB",      # 浅蓝色
"mixed" = "#808080"     # 灰色
)

# 过滤出在数据中实际存在的组织类型
present_tissues <- unique(plot_data$tissue)
available_colors <- tissue_colors[names(tissue_colors) %in% present_tissues]

# 创建基础图
p <- ggplot(plot_data, aes(x = diff_obs_perm, y = neg_log10_p)) +
# 背景线
geom_hline(yintercept = -log10(p_threshold), 
            linetype = "dashed", color = "red", alpha = 0.5) +
geom_vline(xintercept = 0, 
            linetype = "dashed", color = "gray", alpha = 0.5) +
# 点 - 按组织类型着色，统一形状
geom_point(aes(color = tissue), 
            shape = 19,  # 统一使用圆形
            size = point_size, alpha = alpha) +
# 颜色设置
scale_color_manual(values = available_colors) +
# 标签和主题
labs(
    title = "Moran's I: Observed vs Permutation",
    x = "Difference (Observed - Permutation Mean)",
    y = "-log10(p-value)",
    color = "Tissue Type"
) +
theme_minimal() +
theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
)

# 添加统计信息 - 每个组织类型的Yes/No放在同一行
stats_df <- plot_data %>%
group_by(tissue) %>%
summarise(
    total = n(),
    significant = sum(Significant == "Yes"),
    nonsignificant = sum(Significant == "No"),
    .groups = "drop"
) %>%
mutate(
    label = sprintf("%s: Yes=%d, No=%d (Total=%d)", 
                    tissue, significant, nonsignificant, total)
)

p


#Extended  Data Fig8c -----
sub = subset(sp,development == 'unknown',invert=T)
table(sub$new_mimer,sub$file) %>% as.data.frame.array() %>% as.data.frame() -> plot_df


melted_data <- plot_df %>% 
  rownames_to_column("rownames") %>%
  melt()
melted_data %>% 
  group_by(variable) %>%
  mutate(sum = sum(value),
         ratio = value/sum) -> spot_type_ratio 
spot_type_ratio %>% filter(rownames=="mimer") %>% arrange(desc(ratio)) -> rank_mimer
spot_type_ratio$variable = factor(spot_type_ratio$variable,levels = rank_mimer$variable)
spot_type_ratio

spot_type_ratio$rownames = factor(spot_type_ratio$rownames,levels = c("mimer","CD8.C6.PDCD1",'CD4.C7.OX40','FB.C3.COL1A1','Endo.C3.RGCC',"Mac.C2.SPP1","other"))
cols_mimer = c("mimer" = 'brown','other'='grey','CD4.C7.OX40' = "#1fb1aa",
               "Mac.C2.SPP1" = "#c2bc8d","CD8.C6.PDCD1" = "#569c9d",
               "Endo.C3.RGCC" = "#775095","FB.C3.COL1A1" = "#b96c35")

theme_niwot <- function(){
  theme(
    legend.key=element_blank(),   # 图例键为空
    legend.text = element_text(color="black",size=10), # 定义图例文本
    legend.spacing.x=unit(0.1,'cm'), # 定义文本书平距离
    legend.key.width=unit(0.5,'cm'), # 定义图例水平大小
    legend.key.height=unit(0.5,'cm'), # 定义图例垂直大小
    legend.background=element_blank()) 
}
spot_type_ratio %>% 
  # filter(rownames!="other") %>%
  ggplot(.,aes(x=variable,y = ratio))+
  geom_col(aes(fill = rownames),size=2,width = 0.7#,theta = "x",start=0
  ) +
  scale_fill_manual(values = cols_mimer)+
  # coord_polar()+
  # theme_minimal()+
  theme(
    axis.title = element_blank(),
    # axis.ticks = element_blank(),
    axis.text.y = element_text(size=14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.direction = 'vertical',
    panel.grid.major=element_blank(), 
    panel.grid.minor=element_blank(), 
    panel.background = element_blank()
  )+labs(title ="Epithelium/cancer dominted area",fill = 'spot type') +
  theme_niwot()-> bar
bar

sub@meta.data %>% 
  mutate(new_dev = paste0(file,'_',development)) %>%
  distinct(new_dev,.keep_all = TRUE) %>% 
  dplyr::select(file,development,new_dev) -> m1 
m1$file = factor(m1$file,
            levels = rank_mimer$variable)  
m1$development = factor(m1$development,
                        levels = c('Nor','Hyp','MiD','MoD','SD&CA','ICA','MCA'))
m1 %>%
  mutate(group = "development") %>%
  ggplot(aes(file, group, fill = development)) +
  geom_tile(aes(width = 0.7, height = 0.2), position = position_dodge(width = 0.5)) +
  scale_y_discrete(expand = c(0, 0), position = "right") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = domin_cols)+
  # theme_void() +
  theme(axis.title = element_blank(),
        # axis.ticks = element_blank(),
        axis.text.y = element_text(size=14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.direction = 'vertical',
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background = element_blank()) +
  theme_niwot() -> group3

library(aplot)
bar %>% insert_bottom(group3,height = 0.02) 
ggsave("spot_ratio_in_slide_dominted.pdf",width = 9,height = 4)

{
  table(sub$new_mimer,sub$development) %>% as.data.frame.array() %>% as.data.frame() -> plot_df

  melted_data <-plot_df[,1:7] %>% 
    rownames_to_column("rownames") %>%
    melt()
  melted_data %>% 
    # filter(rownames != 'other') %>%
    group_by(variable) %>%
    mutate(sum = sum(value),
           ratio = value/sum) -> spot_type_ratio 
  spot_type_ratio$rownames = factor(spot_type_ratio$rownames,levels = c("mimer","CD8.C6.PDCD1",'CD4.C7.OX40','FB.C3.COL1A1','Endo.C3.RGCC',"Mac.C2.SPP1","other"))

  spot_type_ratio %>% 
    # filter(rownames!="other") %>%
    ggplot(.,aes(x=variable,y = ratio))+
    geom_hline(
      aes(yintercept = y), 
      data.frame(y = c(0:4) * 0.25),
      color = "lightgrey"
    ) +
    geom_col(aes(fill = rownames),size=2,width = 0.7#,theta = "x",start=0
    ) +
    geom_segment(
      aes(
        x = reorder(variable, ratio),
        y = 0,
        xend = reorder(variable, ratio),
        yend = 1
      ),
      linetype = "dashed",
      color = "gray12"
    )+
    scale_fill_manual(values =cols_mimer)+
    coord_polar(start = 0)+
    annotate(
      "text",
      x = rep(0.5,5),
      y = seq(0,1,0.25),
      label = seq(0,1,0.25))+
    # coord_radial(start =0,end=1.8*pi ,r_axis_inside=TRUE)+
    theme_minimal()+
    theme(
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(color = "gray12", size = 12),
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom"
    )+labs(title = "Epithelium/cancer\ndominted area",fill = 'spot type')
  ggsave("Dominted_spot_ratio_bar_with_other.pdf",width = 5,height = 5)

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

# 定义细胞类型映射（去掉_above_median后缀）
cell_type_mapping <- c(
  "CD8_C6_CD39_above_median" = "CD8_C6_CD39",
  "CD4_C7_OX40_above_median" = "CD4_C7_OX40", 
  "Mac_C2_SPP1_above_median" = "Mac_C2_SPP1",
  "FB_C3_COL1A1_above_median" = "FB_C3_COL1A1",
  "Endo_C3_RGCC_above_median" = "Endo_C3_RGCC"
)

# 设置tissue映射和顺序
tissue_map <- c(A = "Adj", N = "N", T = "T", Ln = "LnN", Lp = "LnP")
tissue_levels <- c("N", "Adj", "T", "LnP", "LnN")
tissue_colors <- c(
  "N" = "#6888B9",
  "Adj" = "#9AB292", 
  "T" = "#C09DB0",
  "LnP" = "#CC834F",
  "LnN" = "#DCC13A"
)

# 定义细胞类型颜色和顺序（去掉_above_median后缀）
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

# 定义stage颜色
stage_colors <- c(
  "Nor" = "#1B9E77",           # 浅蓝
  "Hyp" = "#D95F02",           # 中蓝
  "MiD" = "#666666",           # 深蓝
  "MoD" = "#E7298A",           # 深蓝黑
  "SD&CA" = "#66A61E",         # 浅黄
  "ICA" = "#E6AB02",           # 橙色
  "MCA" = "#A6761D",            # 深橙
  "Nor-ME" = "#1B9E77",           # 浅蓝
  "Hyp-ME" = "#D95F02",           # 中蓝
  "MiD-ME" = "#666666",           # 深蓝
  "MoD-ME" = "#E7298A",           # 深蓝黑
  "SD&CA-ME" = "#66A61E",         # 浅黄
  "ICA-ME" = "#E6AB02",           # 橙色
  "MCA-ME" = "#A6761D"            # 深橙
)

# 函数：为细胞定义类型
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

# 函数：分析并绘图
analyze_and_plot <- function(dt_subset, title, filename_prefix, group_var) {
  # 添加tissue列
  dt_subset[, tissue := tissue_map[tissue_SPT]]
  dt_subset <- dt_subset[tissue %in% tissue_levels]
  
  # 定义细胞类型
  cell_types <- define_cell_type(dt_subset)
  dt_subset <- merge(dt_subset, cell_types, by = c("cellID", "file", "patient", "tissue_SPT"))
  dt_subset[, tissue := factor(tissue, levels = tissue_levels)]
  dt_subset[, cell_type := factor(cell_type, levels = cell_type_levels)]
  
  # 计算样本组成
  sample_comp <- dt_subset[, .(count = .N), by = .(file, tissue, patient, cell_type)]
  sample_totals <- sample_comp[, .(total = sum(count)), by = .(file, tissue, patient)]
  sample_comp <- merge(sample_comp, sample_totals, by = c("file", "tissue", "patient"))
  sample_comp[, percentage := count / total * 100]
  
  # 创建样本ID并排序
  sample_comp[, sample_id := paste(patient, file, tissue, sep = "_")]
  mimer_percent <- sample_comp[cell_type == "MIMER", .(sample_id, mimer_percent = percentage)]
  sample_comp <- merge(sample_comp, mimer_percent, by = "sample_id", all.x = TRUE)
  sample_comp[is.na(mimer_percent), mimer_percent := 0]
  
  # 计算stage/ME组成
  stage_comp <- dt_subset[, .(count = .N), by = .(file, tissue, patient, get(group_var))]
  setnames(stage_comp, "get", "stage")
  stage_totals <- stage_comp[, .(total = sum(count)), by = .(file, tissue, patient)]
  stage_comp <- merge(stage_comp, stage_totals, by = c("file", "tissue", "patient"))
  stage_comp[, percentage := count / total * 100]
  stage_comp[, sample_id := paste(patient, file, tissue, sep = "_")]
  
  # 合并排序信息
  stage_comp <- merge(stage_comp, unique(sample_comp[, .(sample_id, mimer_percent)]), by = "sample_id")
  
  # 设置因子顺序
  stage_levels <- names(stage_colors)[names(stage_colors) %in% unique(stage_comp$stage)]
  if(length(stage_levels) > 0) {
    stage_comp[, stage := factor(stage, levels = stage_levels)]
  }
  
  # 确定样本顺序
  sample_order <- unique(sample_comp[order(-mimer_percent)]$sample_id)
  sample_comp[, sample_id := factor(sample_id, levels = sample_order)]
  stage_comp[, sample_id := factor(sample_id, levels = sample_order)]
  
  sample_comp <- sample_comp[!is.na(sample_id)]
  stage_comp <- stage_comp[!is.na(sample_id)]
  
  # 创建样本信息用于tissue颜色条
  sample_info <- unique(sample_comp[, .(sample_id, tissue)])
  sample_info[, tissue_color := tissue_colors[as.character(tissue)]]
  
  # 绘制细胞类型柱状图（第一排）
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
  
  # 绘制stage/ME柱状图（第二排）
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
  
  # 绘制tissue颜色条（第三排）
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
  
  # 组合图形，按照8:1:1的比例
  combined_plot <- p_cell / p_stage / p_tissue + plot_layout(heights = c(8, 0.53, 0.7))
  
  print(combined_plot)
  # 保存
  ggsave(paste0(filename_prefix, "_", group_var, ".pdf"),
         plot = combined_plot, 
         width = 16, height = 6)
  
  # 保存数据
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

# 函数：绘制极坐标图
plot_polar_coord <- function(dt_subset, group_var, title_suffix) {
  # 计算每个group_var中各种细胞类型的数量
  stage_counts <- dt_subset[, .(count = .N), by = .(get(group_var), cell_type)]
  setnames(stage_counts, "get", "stage")
  
  # 计算每个stage的总数和比例
  stage_totals <- stage_counts[, .(total = sum(count)), by = stage]
  stage_counts <- merge(stage_counts, stage_totals, by = "stage")
  stage_counts[, ratio := count / total * 100]
  
  # 设置阶段顺序
  stage_levels <- c('Nor', 'Hyp', 'MiD', 'MoD', 'SD&CA', 'ICA', 'MCA')
  if (group_var == "ME") {
    stage_levels <- paste0(stage_levels, "-ME")
  }
  
  # 只保留存在的stage
  stage_levels <- stage_levels[stage_levels %in% unique(stage_counts$stage)]
  stage_counts[, stage := factor(stage, levels = stage_levels)]
  
  # 设置细胞类型顺序
  stage_counts[, cell_type := factor(cell_type, levels = cell_type_levels)]
  
  # 定义极坐标图的颜色
  polar_cols <- cell_type_colors
  
  # 绘制极坐标图
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
  
  # 绘制不包含Other的极坐标图
  stage_counts_no_other <- stage_counts[cell_type != "Other"]
  
  if (nrow(stage_counts_no_other) > 0) {
    # 计算最大比例以设置y轴
    # max_ratio <- max(stage_counts_no_other$ratio, na.rm = TRUE)
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

# 主分析
dt <- mimer_median_enrich_development_ME

# 1. ME分析（只使用ME != "unknown"的数据）
dt_me <- dt[ME != "unknown"]
plot_me <- analyze_and_plot(
  dt_me, 
  "Cell Type Composition in Microenvironment",
  "ME_cell_type_composition",
  "ME"  # 指定使用ME列
)

# 绘制ME区域的极坐标图
if (!is.null(plot_me$dt_subset)) {
  polar_plots_me <- plot_polar_coord(plot_me$dt_subset, "ME", "ME Regions")
  
  # 保存极坐标图
  ggsave("ME_cell_type_polar.pdf", 
         plot = polar_plots_me$polar_plot, 
         width = 8, height = 8)
  
  if (!is.null(polar_plots_me$polar_plot_no_other)) {
    ggsave("ME_cell_type_polar_no_other.pdf", 
           plot = polar_plots_me$polar_plot_no_other, 
           width = 8, height = 8)
  }
}


# 2. Development分析（只使用development != "unknown"的数据）
dt_dev <- dt[development != "unknown"]
plot_dev <- analyze_and_plot(
  dt_dev,
  "Cell Type Composition in Epi/Tumor",
  "development_cell_type_composition",
  "development"  # 指定使用development列
)

# 绘制Development区域的极坐标图
if (!is.null(plot_dev$dt_subset)) {
  polar_plots_dev <- plot_polar_coord(plot_dev$dt_subset, "development", "Development Regions")
  # 保存极坐标图
  ggsave("development_cell_type_polar.pdf",
         plot = polar_plots_dev$polar_plot,
         width = 8, height = 8)
  
  if (!is.null(polar_plots_dev$polar_plot_no_other)) {
    ggsave("development_cell_type_polar_no_other.pdf",
           plot = polar_plots_dev$polar_plot_no_other,
           width = 8, height = 8)
  }
}


#Extended  Data Fig8d -----
sub = subset(sp,ME == 'unknown',invert=T)

for (i in 1:nrow(df)) {
  SpatialDimPlot(sub,group.by = "new_mimer",
                 images = df$image[i],crop=F,image.alpha = 0.7,
                 # pt.size.factor = 1.5,
                 stroke = NA) + labs(title = df$file[i]) +
    scale_fill_manual(values = cols_mimer ) &
    guides(fill = guide_legend(override.aes = list(size=5)))
  ggsave(filename = paste0("/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhangzhichao/Figure/Mimer/mimer_spot_vs_only/ME/",df$file[i],"_spot_single.pdf"),width = 4,height = 5)
}

sub = subset(sp,development == 'unknown',invert=T)
table(sub$new_mimer,sub$development) %>% as.data.frame.array() %>% as.data.frame() -> plot_df


melted_data <-plot_df[,1:7] %>% 
  rownames_to_column("rownames") %>%
  melt()
melted_data %>% 
  group_by(variable) %>%
  mutate(sum = sum(value),
         ratio = value/sum) -> spot_type_ratio 
spot_type_ratio$rownames = factor(spot_type_ratio$rownames,levels = c("mimer","PDCD1_only",'OX40_only','COL1A1_only','RGCC_only',"SPP1_only","other"))
spot_type_ratio %>% 
  filter(rownames!="other") %>%
  ggplot(.,aes(x=variable,y = ratio))+
  geom_hline(
    aes(yintercept = y), 
    data.frame(y = c(0:3) * 0.1),
    color = "lightgrey"
  ) +
  geom_col(aes(fill = rownames),size=2,width = 0.7#,theta = "x",start=0
  ) +
  geom_segment(
    aes(
      x = reorder(variable, ratio),
      y = 0,
      xend = reorder(variable, ratio),
      yend = 0.4
    ),
    linetype = "dashed",
    color = "gray12"
  )+
  coord_polar(start = 0)+
  annotate(
    "text",
    x = rep(0.5,5),
    y = seq(0,0.4,0.1),
    label = seq(0,0.4,0.1))+
  # coord_radial(start =0,end=1.8*pi ,r_axis_inside=TRUE)+
  theme_minimal()+
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(color = "gray12", size = 12),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )+labs(title = "Epithelium/cancer\ndominted area",fill = 'spot type')
ggsave("/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhangzhichao/Figure/Mimer/mimer_spot_vs_only/dominated_spot_ratio_bar.pdf",width = 5,height = 5)
spot_type_ratio %>% 
  # filter(rownames!="other") %>%
  ggplot(.,aes(x=variable,y = ratio))+
  geom_hline(
    aes(yintercept = y), 
    data.frame(y = c(0:4) * 0.25),
    color = "lightgrey"
  ) +
  geom_col(aes(fill = rownames),size=2,width = 0.7#,theta = "x",start=0
  ) +
  geom_segment(
    aes(
      x = reorder(variable, ratio),
      y = 0,
      xend = reorder(variable, ratio),
      yend = 1
    ),
    linetype = "dashed",
    color = "gray12"
  )+
  scale_fill_manual(values =cols_mimer)+
  coord_polar(start = 0)+
  annotate(
    "text",
    x = rep(0.5,5),
    y = seq(0,1,0.25),
    label = seq(0,1,0.25))+
  # coord_radial(start =0,end=1.8*pi ,r_axis_inside=TRUE)+
  theme_minimal()+
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(color = "gray12", size = 12),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )+labs(title = "Epithelium/cancer\ndominted area",fill = 'spot type')
ggsave("/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhangzhichao/Figure/Mimer/mimer_spot_vs_only/dominated_spot_ratio_bar_with_other.pdf",width = 5,height = 5)
#Extended  Data Fig8e -----
 p = subset(sp,development=='unknown',invert=T)
  p$new_mimer_region = ifelse(p$new_mimer %in% 'mimer','mimer','other')
  table(p$new_mimer_region,p$file) %>% as.data.frame.array() ->plot_df
  plot_df[plot_df!=0] = 1
  plot_df
  meta_plot = p@meta.data %>% dplyr::distinct(file,.keep_all = T) %>% 
    select(file,patient) %>% arrange(patient)
  rownames(meta_plot) = meta_plot$file
  plot_df = plot_df[,match(rownames(meta_plot),colnames(plot_df))]
  meta_plot$file = NULL
  anno_cols = c("#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", 
                "#59A14F", "#8CD17D", "#B6992D", "#F1CE63", "#499894", 
                "#86BCB6", "#E15759", "#FF9D9A", "#79706E", "#BAB0AC", 
                "#D37295", "#FABFD2", "#B07AA1", "#D4A6C8", "#9D7660", 
                "#D7B5A6")[1:15]
  names(anno_cols)=unique(meta_plot$patient)
  pdf("~/workspace/esca/picture/pheatmap_number_of_MIMER_spot_in_slide_development.pdf",
      width = 20 ,height = 10)
  pheatmap::pheatmap(plot_df,cluster_rows = F,cluster_cols = F,
                     cellwidth = 20,cellheight = 20,color = c("white","brown"),
                     annotation_col = meta_plot,
                     gaps_col = c(4,8,12,16,20,24,29,34,37,41,45,49,53,58),
                     annotation_colors = list(patient = anno_cols ),
                     breaks = c(0, 0.5, 1))
  dev.off()
  
  
  p = subset(sp,ME=='unknown',invert=T)
  p$new_mimer_region = ifelse(p$new_mimer %in% 'mimer','mimer','other')
  table(p$new_mimer_region,p$file) %>% as.data.frame.array() ->plot_df
  plot_df[plot_df!=0] = 1
  plot_df
  meta_plot = p@meta.data %>% dplyr::distinct(file,.keep_all = T) %>% 
    select(file,patient) %>% arrange(patient)
  rownames(meta_plot) = meta_plot$file
  plot_df = plot_df[,match(rownames(meta_plot),colnames(plot_df))]
  meta_plot$file = NULL
  anno_cols = c("#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", 
                "#59A14F", "#8CD17D", "#B6992D", "#F1CE63", "#499894", 
                "#86BCB6", "#E15759", "#FF9D9A", "#79706E", "#BAB0AC", 
                "#D37295", "#FABFD2", "#B07AA1", "#D4A6C8", "#9D7660", 
                "#D7B5A6")[1:14]
  names(anno_cols)=unique(meta_plot$patient)
  pdf("~/workspace/esca/picture/pheatmap_number_of_MIMER_spot_in_slide_ME.pdf",
      width = 20 ,height = 10)
  pheatmap::pheatmap(plot_df,cluster_rows = F,cluster_cols = F,
                     cellwidth = 20,cellheight = 20,color = c("white","brown"),
                     annotation_col = meta_plot,
                     gaps_col = c(4,8,12,16,20,24,28,33,36,40,44,48,52),
                     annotation_colors = list(patient = anno_cols ),
                     breaks = c(0, 0.5, 1))
  dev.off()
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
# 计算均值
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
#节点数据
nodes <- data.frame(name = unique(union(mst_cons_edge$from, mst_cons_edge$to)))

nodes$number = brain_cell_class_new[brain_cell_class_new$id %in% nodes$name,"freq.Freq"]
nodes
class(nodes)
#边数据
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
#Extended  Data Fig8h -----
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


# E8ij
# 加载必要的包
library(ggplot2)
library(scales)
library(data.table)
library(dplyr)
library(tidyr)

# 设置分析类型：'development' 或 'ME'
analysis_type <- "development"  # 可更改为 "development" 或 "ME"

# 设置数据表为data.table格式
setwd('C:/Users/Ron/Desktop/Figs/E8/')
mimer_cor_development_ME <- fread('C:/Users/Ron/Desktop/Figs/E8/E8_ij_mimer_cor_development_ME.csv')
setDT(mimer_cor_development_ME)

# 根据分析类型处理数据
if (analysis_type == "ME") {
  # 使用ME列
  type_col <- "ME"
  filtered_types <- c("Nor-ME", "Hyp-ME", "MiD-ME", "MoD-ME", "SD&CA-ME", "ICA-ME")
} else {
  # 使用development列
  type_col <- "development"
  filtered_types <- c("Nor", "Hyp", "MiD", "MoD", "SD&CA", "ICA")
}

# 阈值定义
thresholds <- c("top5_perc_mimer", "top10_perc_mimer", "top20_perc_mimer")
threshold_names <- c("Top 5%", "Top 10%", "Top 20%")
threshold_colors <- c("Top 5%" = "#FF7F0E", "Top 10%" = "#3fa02c", "Top 20%" = "#1F77B4")

# 计算每个病人每个类型每个阈值的MIMER比例
patient_stats_list <- list()

for (patient_id in unique(mimer_cor_development_ME$patient)) {
  # 获取该病人的数据
  patient_data <- mimer_cor_development_ME[patient == patient_id & get(type_col) %in% filtered_types]
  
  if (nrow(patient_data) > 0) {
    for (current_type in filtered_types) {
      # 获取该病人该类型的所有细胞
      type_data <- patient_data[get(type_col) == current_type]
      
      if (nrow(type_data) > 0) {
        for (i in seq_along(thresholds)) {
          threshold <- thresholds[i]
          threshold_name <- threshold_names[i]
          
          # 计算MIMER+细胞数
          mimer_cells <- sum(type_data[[threshold]] == "MIMER", na.rm = TRUE)
          total_cells <- nrow(type_data)
          prop <- ifelse(total_cells > 0, mimer_cells / total_cells, 0)
          
          # 存储结果
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

# 合并所有结果
patient_stats <- rbindlist(patient_stats_list)

# 设置因子顺序
patient_stats$type <- factor(patient_stats$type, levels = filtered_types)
patient_stats$threshold <- factor(patient_stats$threshold, levels = threshold_names)

# 计算每个类型在每个病人中的样本数
type_counts <- mimer_cor_development_ME[
  get(type_col) %in% filtered_types, 
  .(n_cells = .N), 
  by = .(patient, get(type_col))
]
setnames(type_counts, "get", "type")
type_counts$type <- factor(type_counts$type, levels = filtered_types)

# 计算每个类型出现在多少个病人中
type_patient_counts <- type_counts[, .(n_patients = .N), by = type]
type_patient_counts <- type_patient_counts[order(type)]

# 创建x轴标签（包含病人数）
x_labels <- sapply(filtered_types, function(t) {
  n_patients <- type_patient_counts[type == t, n_patients]
  if (length(n_patients) > 0 && n_patients > 0) {
    return(paste0(t, "\n(n=", n_patients, ")"))
  } else {
    return(t)
  }
})

# 绘制折线图
p <- ggplot(patient_stats, aes(x = type, y = prop, group = threshold, color = threshold)) +
  # 为每个病人绘制折线
  geom_line(aes(group = interaction(patient, threshold)), alpha = 0.3, linewidth = 0.5) +
  # 计算并绘制所有病人的均值线
  stat_summary(
    aes(group = threshold), 
    fun = mean, 
    geom = "line", 
    linewidth = 1.5
  ) +
  # 添加均值点
  stat_summary(
    aes(group = threshold), 
    fun = mean, 
    geom = "point", 
    size = 3
  ) +
  # 分面显示
  facet_wrap(~ patient, scales = 'free_y') +
  # 颜色设置
  scale_color_manual(values = threshold_colors) +
  # 设置x轴标签
  scale_x_discrete(limits = filtered_types, labels = x_labels) +
  # 设置y轴
  scale_y_continuous(
    labels = percent_format(accuracy = 1)#,
    #limits = c(0, 1)
  ) +
  # 标签和标题
  labs(
    title = ifelse(analysis_type == "ME", 
                   "MIMER+ Spot Proportion Across ME Types by Patient", 
                   "MIMER+ Spot Proportion Across Epi/Cancer Types by Patient"),
    x = ifelse(analysis_type == "ME", "ME Type", "Epi/Cancer Type"),
    y = "Proportion of MIMER+ Spots",
    color = "Threshold"
  ) +
  # 主题设置
  theme_minimal(base_size = 12) +
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

# 显示图形
options(repr.plot.width=16,repr.plot.height=10)
print(p)

# # 获取病人数
# n_patients <- length(unique(patient_stats$patient))
# cat(sprintf("共绘制了 %d 个病人的折线图\n", n_patients))
# cat("每个病人包含的类型:\n")
# print(table(type_counts[, .(patient, type)]))

# # 保存PDF
# pdf_file <- ifelse(analysis_type == "ME", 
#                    "ME_mimer_perc_by_patient_lines.pdf", 
#                    "development_mimer_perc_by_patient_lines.pdf")
# pdf(pdf_file, width = 12, height = 8)  # 调整宽度以容纳多个病人
# print(p)
# dev.off()
# cat(sprintf("PDF已保存到: %s\n", pdf_file))

#Extended  Data Fig8i -----
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

#Extended  Data Fig8j -----
sp@meta.data$cc_15_combine = sp@meta.data$cc_15
sp@meta.data$cc_15_combine = ifelse(sp$cc_15_combine %in% c("CC4","CC7","CC1","CC6","CC5","CC14"),
                                    "Mimer_resident_niches",sp$cc_15_combine)

DimPlot(sp,group.by='cc_15_combine',raster=T,reduction = "umap.ischia15",label = T,repel = T)
development = subset(sp,development =='unknown',invert=T)
SpatialDimPlot(development,group.by  = c("development"),images = df[df$file=='QP','image'],stroke = NA,crop = F,pt.size.factor = 1.3)+
   scale_fill_manual(values = c("SD&CA"="#66A61E"))
ggsave(filename = "development_QP_SpatialDimPlot.pdf",width = 6,height = 6)


development@meta.data %>%
  group_by(file,development,cc_15_combine) %>% 
  summarise(number = n()) %>%
  ungroup() %>% 
  group_by(file,development) %>%
  mutate(sum = sum(number),
         ratio = number/sum) -> dev_df

head(dev_df)
library(ggrastr)
my_levels<-c('Nor','Hyp','MiD','MoD','SD&CA','ICA','MCA')
# Idents(p)<-factor(Idents(p),levels = my_levels)
mycolor<-c("#1B9E77","#D95F02","#666666","#E7298A","#66A61E","#E6AB02","#A6761D")
my_comparison<-list(c('Nor','Hyp'),c('Hyp','MiD'),c('MiD','MoD'),c('MoD','SD&CA'),c('SD&CA','ICA'),c('ICA','MCA'))
dev_df$development = factor(dev_df$development,
                            levels = my_levels)

dev_df %>% 
  filter(file != 'N3NT') %>% 
  filter(cc_15_combine =="Mimer_resident_niches") %>%
  ggplot(.,aes(x=development,y=ratio,color=development))+
  geom_boxplot()+
  geom_jitter(width=0.2)+ 
  scale_fill_manual(values = mycolor)+
  scale_color_manual(values = mycolor)+
  stat_compare_means(comparisons =my_comparison )+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0,size = 14), 
        legend.position= "none",
        strip.background = element_rect(color="black",fill =  "#76B7B2", linetype = "solid",size=1),
        strip.text = element_text( color = "black",hjust = 0.5, size = 18,),
        plot.title=element_text(size="20", color="brown",hjust = 0.5),
        axis.text.y = element_text(size = 14) ,
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color="black", size=1, linetype="solid")
  )+labs(x='',y="The proportion of Mimer resident niches in slide")

ggsave(filename = 'Mimer_resident_niches_ratio_in_slide.pdf',width = 6,height = 5)


ME = subset(sp,ME=='unknown',invert=T)
ME@meta.data %>%
  group_by(file,ME,cc_15_combine) %>% 
  summarise(number = n()) %>%
  ungroup() %>% 
  group_by(file,ME) %>%
  mutate(sum = sum(number),
         ratio = number/sum) -> dev_df

head(dev_df)
library(ggrastr)
my_levels<-c('Nor-ME','Hyp-ME','MiD-ME','MoD-ME','SD&CA-ME','ICA-ME','MCA-ME')
# Idents(p)<-factor(Idents(p),levels = my_levels)
mycolor<-c("#1B9E77","#D95F02","#666666","#E7298A","#66A61E","#E6AB02","#A6761D")
my_comparison<-list(c('Nor-ME','Hyp-ME'),c('Hyp-ME','MiD-ME'),c('MiD-ME','MoD-ME'),c('MoD-ME','SD&CA-ME'),c('SD&CA-ME','ICA-ME'),c('ICA-ME','MCA-ME'))
dev_df$ME = factor(dev_df$ME,
                            levels = my_levels)
dev_df %>%
  filter(file != 'N3NT') %>% 
  filter(cc_15_combine =="Mimer_resident_niches") %>%
  ggplot(.,aes(x=ME,y=ratio,color=ME))+
  geom_boxplot(aes(fill=ME))+geom_point()+
  scale_fill_manual(values = mycolor)+
  scale_color_manual(values = mycolor)+
  stat_compare_means(comparisons =my_comparison )

dev_df %>% 
  filter(file != 'N3NT') %>% 
  filter(cc_15_combine =="Mimer_resident_niches") %>%
  ggplot(.,aes(x=ME,y=ratio,color=ME))+
  geom_boxplot()+
  geom_jitter(width=0.2)+ 
  scale_fill_manual(values = mycolor)+
  scale_color_manual(values = mycolor)+
  stat_compare_means(comparisons =my_comparison )+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,size = 14,vjust = 0.5,hjust = 1), 
        legend.position= "none",
        strip.background = element_rect(color="black",fill =  "#76B7B2", linetype = "solid",size=1),
        strip.text = element_text( color = "black",hjust = 0.5, size = 18,),
        plot.title=element_text(size="20", color="brown",hjust = 0.5),
        axis.text.y = element_text(size = 14) ,
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color="black", size=1, linetype="solid")
  )+labs(x='',y="The proportion of Mimer resident niches in slide")


ggsave(filename = 'Mimer_resident_niches_ratio_in_slide_ME.pdf',width = 6,height = 5)

####E8k
# E8_k
library(ggplot2)
library(scales)
library(data.table)
library(dplyr)
library(tidyr)
library(gridExtra)  # 用于多图布局

setwd('C:/Users/Ron/Desktop/Figs/E8/')
Module_score_add_me_dev_tissue <- fread('C:/Users/Ron/Desktop/Figs/E8/E8_k_Module_score_add_me_dev_tissue.csv')
# 设置数据表为data.table格式
setDT(Module_score_add_me_dev_tissue)

# 定义要分析的模块
modules <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7")
module_names <- c("Module 1", "Module 2", "Module 3", "Module 4", "Module 5", "Module 6", "Module 7")

# 设置分析类型：同时分析 'development' 和 'ME'
analysis_types <- c("development", "ME")

# 1. 检查原始development和ME列的取值
cat("原始development列的取值分布:\n")
print(table(Module_score_add_me_dev_tissue$development))

cat("\n原始ME列的取值分布:\n")
print(table(Module_score_add_me_dev_tissue$ME))

# 定义分析类型对应的列和类型列表
type_defs <- list(
  development = list(
    type_col = "development",
    all_types = c("Nor", "Hyp", "MiD", "MoD", "SD&CA", "ICA", "MCA"),
    nor_type = "Nor"
  ),
  ME = list(
    type_col = "ME",
    all_types = c("Nor-ME", "Hyp-ME", "MiD-ME", "MoD-ME", "SD&CA-ME", "ICA-ME", "MCA-ME"),
    nor_type = "Nor-ME"
  )
)

# 过滤掉unknown和MCA类型
for (type in names(type_defs)) {
  if (type == "ME") {
    type_defs[[type]]$filtered_types <- type_defs[[type]]$all_types[!type_defs[[type]]$all_types %in% c("unknown", "MCA-ME")]
  } else {
    type_defs[[type]]$filtered_types <- type_defs[[type]]$all_types[!type_defs[[type]]$all_types %in% c("unknown", "MCA")]
  }
}

# 定义共同类型（去掉-ME后缀）
common_types <- c("Nor", "Hyp", "MiD", "MoD", "SD&CA", "ICA")

# 定义阈值
threshold_percentile <- 0.95  # 使用top5%阈值
threshold_name <- "Top 5%"

# 初始化存储所有模块结果的列表
all_modules_results <- list()

# 对每个模块进行分析
for (module_idx in seq_along(modules)) {
  module <- modules[module_idx]
  module_name <- module_names[module_idx]
  
  cat(sprintf("\n=== 分析模块: %s (%s) ===\n", module, module_name))
  
  # 复制数据，避免修改原始数据
  data_copy <- copy(Module_score_add_me_dev_tissue)
  
  # 计算该模块的top5%阈值
  module_scores <- data_copy[[module]]
  threshold_value <- quantile(module_scores, threshold_percentile, na.rm = TRUE)
  
  cat(sprintf("模块 %s 的top5%%阈值: %.6f\n", module, threshold_value))
  cat(sprintf("模块 %s 分数范围: %.6f 到 %.6f\n", module, 
              min(module_scores, na.rm = TRUE), 
              max(module_scores, na.rm = TRUE)))
  
  # 创建top5标记列
  data_copy[, top5_status := ifelse(get(module) >= threshold_value, "Top5", "Non-Top5")]
  
  # 初始化存储所有分析类型结果的列表
  module_results <- list()
  
  # 对每种分析类型进行处理
  for (analysis_type in analysis_types) {
    cat(sprintf("\n分析类型: %s\n", analysis_type))
    
    type_info <- type_defs[[analysis_type]]
    type_col <- type_info$type_col
    filtered_types <- type_info$filtered_types
    nor_type <- type_info$nor_type
    
    # 过滤掉unknown和MCA类型的数据
    base_data <- data_copy[get(type_col) %in% filtered_types]
    
    # 为ME类型创建通用类型名（去掉-ME后缀）
    if (analysis_type == "ME") {
      base_data[, common_type := gsub("-ME", "", get(type_col))]
    } else {
      base_data[, common_type := get(type_col)]
    }
    
    # 计算每个通用类型实际出现在哪些文件中
    type_files <- base_data[, 
                            .(files = unique(file)), 
                            by = common_type
    ]
    setnames(type_files, "common_type", "type")
    
    # 计算每个通用类型的文件数（样本数）
    type_n_files <- type_files[, .(n_files = .N), by = type]
    type_n_files$type <- factor(type_n_files$type, levels = common_types)
    type_n_files <- type_n_files[order(type)]
    
    # 初始化存储结果
    sample_stats_list <- list()
    
    for (current_common_type in common_types) {
      # 获取包含当前通用类型的所有文件
      files_with_type <- unique(base_data[common_type == current_common_type, file])
      
      if (length(files_with_type) > 0) {
        # 对于每个包含当前通用类型的文件，计算Top5比例
        for (current_file in files_with_type) {
          # 获取该文件当前通用类型的所有细胞
          type_cells <- base_data[file == current_file & common_type == current_common_type]
          
          # 计算Top5细胞数
          top5_cells <- sum(type_cells$top5_status == "Top5", na.rm = TRUE)
          
          # 计算总细胞数
          total_cells <- nrow(type_cells)
          
          # 计算比例
          prop <- ifelse(total_cells > 0, top5_cells / total_cells, 0)
          
          # 保存结果
          sample_stats_list[[paste(current_common_type, current_file, sep="_")]] <- data.frame(
            file = current_file,
            common_type = current_common_type,
            top5_cells = top5_cells,
            total_cells = total_cells,
            prop = prop,
            threshold = threshold_name,
            module = module,
            analysis_type = analysis_type
          )
        }
      } else {
        cat(sprintf("  警告: 通用类型 %s 在%s数据中未找到\n", current_common_type, analysis_type))
      }
    }
    
    # 合并所有结果
    if (length(sample_stats_list) > 0) {
      sample_stats <- rbindlist(sample_stats_list)
    } else {
      sample_stats <- data.table(
        file = character(),
        common_type = factor(levels = common_types),
        top5_cells = integer(),
        total_cells = integer(),
        prop = numeric(),
        threshold = character(),
        module = character(),
        analysis_type = character()
      )
    }
    
    # 设置通用类型为因子
    sample_stats$common_type <- factor(sample_stats$common_type, levels = common_types)
    
    # 计算每个通用类型的统计量
    mean_stats <- sample_stats[, 
                               .(mean_prop = mean(prop, na.rm = TRUE),
                                 se_prop = sd(prop, na.rm = TRUE) / sqrt(.N),
                                 median_prop = median(prop, na.rm = TRUE),
                                 sd_prop = sd(prop, na.rm = TRUE),
                                 n_samples = .N,
                                 n_nonzero = sum(prop > 0)), 
                               by = common_type
    ]
    
    # 对于没有数据的通用类型，添加0行
    missing_types <- setdiff(common_types, mean_stats$common_type)
    if (length(missing_types) > 0) {
      for (mt in missing_types) {
        mean_stats <- rbindlist(list(
          mean_stats,
          data.table(
            common_type = mt,
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
    mean_stats[, module := module]  
    mean_stats[, analysis_type := analysis_type]
    
    # 确保通用类型顺序
    mean_stats$common_type <- factor(mean_stats$common_type, levels = common_types)
    mean_stats <- mean_stats[order(common_type)]
    
    # 进行Wilcoxon检验 - 每个类型与Nor比较
    # 获取Nor类型的数据
    nor_data <- sample_stats[common_type == "Nor", prop]
    n_nor <- length(nor_data)
    
    # 初始化存储比较结果的列表
    comparison_results <- list()
    
    # 对每个非Nor类型与Nor进行比较
    for (current_common_type in common_types) {
      if (current_common_type != "Nor") {
        # 获取当前类型的数据
        current_data <- sample_stats[common_type == current_common_type, prop]
        n_current <- length(current_data)
        
        # 进行Wilcoxon检验
        if (n_nor >= 2 && n_current >= 2) {
          test_result <- wilcox.test(current_data, nor_data, exact = FALSE)
          p_value <- test_result$p.value
        } else {
          p_value <- NA
        }
        
        # 计算显著性标记
        significance <- ifelse(is.na(p_value), "NA",
                               ifelse(p_value < 0.001, "***",
                                      ifelse(p_value < 0.01, "**",
                                             ifelse(p_value < 0.05, "*", "ns"))))
        
        # 存储结果
        comparison_results[[current_common_type]] <- list(
          p_value = p_value,
          significance = significance,
          n_nor = n_nor,
          n_current = n_current
        )
      }
    }
    
    # 存储结果
    module_results[[analysis_type]] <- list(
      sample_stats = sample_stats,
      mean_stats = mean_stats,
      comparison_results = comparison_results,
      n_nor = n_nor
    )
  }
  
  # 存储模块结果
  all_modules_results[[module]] <- list(
    development = module_results[["development"]],
    ME = module_results[["ME"]],
    module_name = module_name
  )
}

# 准备绘图
# 创建存储所有图的列表
plot_list <- list()

# 设置图表标题
base_title <- "Top5% Spot Proportion Across Types"
x_label <- "Type"

# 对每个模块创建单独的图
for (module_idx in seq_along(modules)) {
  module <- modules[module_idx]
  module_name <- module_names[module_idx]
  
  results <- all_modules_results[[module]]
  
  if (!is.null(results)) {
    # 合并两种分析类型的mean_stats
    dev_mean_stats <- results$development$mean_stats
    me_mean_stats <- results$ME$mean_stats
    
    # 在mean_stats中添加样本数信息
    dev_mean_stats$type_with_n <- paste0(dev_mean_stats$common_type, " (n=", dev_mean_stats$n_samples, ")")
    me_mean_stats$type_with_n <- paste0(me_mean_stats$common_type, " (n=", me_mean_stats$n_samples, ")")
    
    # 合并数据
    combined_data <- rbindlist(list(dev_mean_stats, me_mean_stats), fill = TRUE)
    combined_data$analysis_type <- factor(combined_data$analysis_type, levels = c("development", "ME"))
    
    # 确保通用类型顺序
    combined_data$common_type <- factor(combined_data$common_type, levels = common_types)
    
    # 为每个通用类型创建标签（包含两种分析类型的样本数）
    type_labels <- sapply(common_types, function(t) {
      dev_n <- if (t %in% dev_mean_stats$common_type) {
        dev_mean_stats[common_type == t, n_samples]
      } else {
        0
      }
      
      me_n <- if (t %in% me_mean_stats$common_type) {
        me_mean_stats[common_type == t, n_samples]
      } else {
        0
      }
      
      return(paste0(t, "\nDev(n=", dev_n, "), ME(n=", me_n, ")"))
    })
    
    # 准备标注数据 - 每个类型与Nor比较的结果
    annotation_data <- data.frame()
    
    # 处理development类型的显著性标记
    for (current_type in common_types) {
      if (current_type != "Nor") {
        comp_result <- results$development$comparison_results[[current_type]]
        
        if (!is.null(comp_result) && !is.na(comp_result$p_value)) {
          # 获取当前类型的位置
          type_pos <- which(common_types == current_type)
          
          # 获取development类型的统计量
          type_mean <- dev_mean_stats[common_type == current_type, mean_prop]
          type_se <- dev_mean_stats[common_type == current_type, se_prop]
          
          if (length(type_mean) > 0 && length(type_se) > 0) {
            y_pos <- type_mean + type_se
            
            # 创建标注数据
            annotation_data <- rbind(
              annotation_data,
              data.frame(
                type = current_type,
                x_pos = type_pos,
                y_pos = y_pos,
                p_value = comp_result$p_value,
                significance = comp_result$significance,
                n_current = comp_result$n_current,
                analysis_type = "development"
              )
            )
          }
        }
      }
    }
    
    # 处理ME类型的显著性标记
    for (current_type in common_types) {
      if (current_type != "Nor") {
        comp_result <- results$ME$comparison_results[[current_type]]
        
        if (!is.null(comp_result) && !is.na(comp_result$p_value)) {
          # 获取当前类型的位置
          type_pos <- which(common_types == current_type)
          
          # 获取ME类型的统计量
          type_mean <- me_mean_stats[common_type == current_type, mean_prop]
          type_se <- me_mean_stats[common_type == current_type, se_prop]
          
          if (length(type_mean) > 0 && length(type_se) > 0) {
            y_pos <- type_mean + type_se
            
            # 创建标注数据
            annotation_data <- rbind(
              annotation_data,
              data.frame(
                type = current_type,
                x_pos = type_pos,
                y_pos = y_pos,
                p_value = comp_result$p_value,
                significance = comp_result$significance,
                n_current = comp_result$n_current,
                analysis_type = "ME"
              )
            )
          }
        }
      }
    }
    
    # 计算y轴上限
    y_max_dev <- max(dev_mean_stats$mean_prop + dev_mean_stats$se_prop, na.rm = TRUE, default = 0)
    y_max_me <- max(me_mean_stats$mean_prop + me_mean_stats$se_prop, na.rm = TRUE, default = 0)
    y_max <- max(y_max_dev, y_max_me, annotation_data$y_pos, na.rm = TRUE, default = 0)
    y_limit <- ifelse(y_max > 0, y_max * 1.15, 0.1)  # 增加上限以容纳显著性标记
    
    # 创建带有轻微错位的x轴位置
    # 将因子类型转换为数值，然后根据分析类型进行微调
    combined_data$x_pos <- as.numeric(combined_data$common_type)
    
    # 为development类型向左偏移0.1，为ME类型向右偏移0.1
    combined_data$x_pos_adjusted <- combined_data$x_pos
    combined_data$x_pos_adjusted[combined_data$analysis_type == "development"] <- 
      combined_data$x_pos[combined_data$analysis_type == "development"] - 0.1
    combined_data$x_pos_adjusted[combined_data$analysis_type == "ME"] <- 
      combined_data$x_pos[combined_data$analysis_type == "ME"] + 0.1
    
    # 为显著性标记也做类似的错位调整
    if (nrow(annotation_data) > 0) {
      annotation_data$x_pos_adjusted <- annotation_data$x_pos
      annotation_data$x_pos_adjusted[annotation_data$analysis_type == "development"] <- 
        annotation_data$x_pos[annotation_data$analysis_type == "development"] - 0.1
      annotation_data$x_pos_adjusted[annotation_data$analysis_type == "ME"] <- 
        annotation_data$x_pos[annotation_data$analysis_type == "ME"] + 0.1
    }
    
    # 创建折线图
    p <- ggplot(combined_data, aes(x = x_pos_adjusted, y = mean_prop, group = analysis_type, color = analysis_type)) +
      # 添加均值点
      geom_point(size = 3) +
      # 添加连接线
      geom_line(linewidth = 1) +
      # 添加误差棒
      geom_errorbar(
        aes(ymin = mean_prop - se_prop, ymax = mean_prop + se_prop),
        width = 0.2,
        linewidth = 0.8
      )
    
    # 只在有显著结果时添加显著性标记
    if (nrow(annotation_data) > 0) {
      # 只保留显著的结果（去掉"ns"和"NA"）
      sig_annotations <- annotation_data[annotation_data$significance %in% c("*", "**", "***"), ]
      
      if (nrow(sig_annotations) > 0) {
        # 为显著性标记添加额外的高度偏移
        y_offset <- max(combined_data$mean_prop + combined_data$se_prop, na.rm = TRUE) * 0.1
        sig_annotations$y_pos <- sig_annotations$y_pos + y_offset
        
        # 只显示显著性符号，不显示样本数
        p <- p + 
          geom_text(
            data = sig_annotations,
            aes(x = x_pos_adjusted, y = y_pos, label = significance, color = analysis_type),
            size = 5,
            fontface = "bold",
            vjust = 0.5,
            show.legend = FALSE
          )
      }
    }
    
    # 完成图表设置
    p <- p +
      # 设置x轴标签
      scale_x_continuous(
        breaks = 1:length(common_types),
        labels = type_labels,
        limits = c(0.5, length(common_types) + 0.5)
      ) +
      # 设置y轴
      scale_y_continuous(
        labels = percent_format(accuracy = 1),
        limits = c(0, y_limit),
        expand = expansion(mult = c(0.05, 0.2))  # 为显著性标记留出更多空间
      ) +
      # 设置颜色
      scale_color_discrete(
        name = "Analysis Type",
        labels = c("development", "ME")
      ) +
      # 标签和标题
      labs(
        title = paste(module_name, "\n", base_title),
        x = x_label, 
        y = "Proportion of Top5% Spots"
      ) +
      # 主题设置
      theme_minimal(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 10, face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(size = 10, face = "bold", margin = margin(r = 10)),
        plot.title = element_text(
          hjust = 0.5, 
          face = "bold", 
          size = 12,
          margin = margin(b = 5)
        ),
        legend.title = element_text(face = "bold"),
        legend.position = "top",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_line(color = "gray95", linewidth = 0.3),
        panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
      )
    
    plot_list[[module]] <- p
  }
}

# 组合所有图到一个大图中
if (length(plot_list) > 0) {
  # 计算合适的布局（3列，3行）
  n_cols <- 3
  n_rows <- ceiling(length(plot_list) / n_cols)
  
  # 创建组合图
  combined_plot <- grid.arrange(
    grobs = plot_list,
    ncol = n_cols,
    nrow = n_rows,
    top = "Combined Analysis - All Modules (Top 5% Threshold)"
  )
  
  # 显示组合图
  print(combined_plot)
}

# 输出所有模块的统计检验结果
cat("\n========== 所有模块的Wilcoxon检验结果 ==========\n")

for (module_idx in seq_along(modules)) {
  module <- modules[module_idx]
  module_name <- module_names[module_idx]
  
  cat(sprintf("\n模块: %s (%s)\n", module, module_name))
  
  # 输出development类型的结果
  cat("=== development分析 ===\n")
  if (!is.null(all_modules_results[[module]]$development)) {
    results <- all_modules_results[[module]]$development
    
    cat(sprintf("  Nor 样本数: %d\n", results$n_nor))
    
    for (current_type in common_types) {
      if (current_type != "Nor") {
        comp_result <- results$comparison_results[[current_type]]
        
        if (!is.null(comp_result)) {
          cat(sprintf("  %s vs Nor:\n", current_type))
          cat(sprintf("    %s样本数: %d\n", current_type, comp_result$n_current))
          cat(sprintf("    p值: %.4e\n", comp_result$p_value))
          cat(sprintf("    显著性: %s\n\n", comp_result$significance))
        }
      }
    }
  }
  
  # 输出ME类型的结果
  cat("=== ME分析 ===\n")
  if (!is.null(all_modules_results[[module]]$ME)) {
    results <- all_modules_results[[module]]$ME
    
    cat(sprintf("  Nor 样本数: %d\n", results$n_nor))
    
    for (current_type in common_types) {
      if (current_type != "Nor") {
        comp_result <- results$comparison_results[[current_type]]
        
        if (!is.null(comp_result)) {
          cat(sprintf("  %s vs Nor:\n", current_type))
          cat(sprintf("    %s样本数: %d\n", current_type, comp_result$n_current))
          cat(sprintf("    p值: %.4e\n", comp_result$p_value))
          cat(sprintf("    显著性: %s\n\n", comp_result$significance))
        }
      }
    }
  }
}

# 保存PDF
tryCatch({
  pdf_file <- "combined_modules_top5_perc_vs_Nor.pdf"
  
  # 设置PDF大小，根据布局调整
  pdf_width <- 15
  pdf_height <- 5 * n_rows
  
  pdf(pdf_file, width = pdf_width, height = pdf_height)
  
  # 重新绘制组合图
  grid.arrange(
    grobs = plot_list,
    ncol = n_cols,
    nrow = n_rows,
    top = "Combined Analysis - All Modules (Top 5% Threshold)"
  )
  
  dev.off()
  cat(sprintf("\nPDF已成功保存到: %s\n", pdf_file))
  cat(sprintf("PDF尺寸: %.1f x %.1f 英寸\n", pdf_width, pdf_height))
}, error = function(e) {
  cat(sprintf("保存PDF时出错: %s\n", e$message))
  cat("当前工作目录:", getwd(), "\n")
  
  # 尝试在临时目录中保存
  pdf_file <- file.path(tempdir(), "combined_modules_top5_perc_vs_Nor.pdf")
  
  pdf_width <- 15
  pdf_height <- 5 * n_rows
  
  pdf(pdf_file, width = pdf_width, height = pdf_height)
  
  grid.arrange(
    grobs = plot_list,
    ncol = n_cols,
    nrow = n_rows,
    top = "Combined Analysis - All Modules (Top 5% Threshold)"
  )
  
  dev.off()
  cat(sprintf("PDF已保存到临时目录: %s\n", pdf_file))
})
