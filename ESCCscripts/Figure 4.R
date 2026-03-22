#Figure 4a -----
load("./permutation_All.RData")
all_long <- tidyr::gather(all_sum_mat, key = "sample", value = "closeness", -"cell")


all_long$sample <- lapply(strsplit(all_long$sample, '[.]'),function(x){res<-x[1];return(res);})%>%unlist()
mycell<-"CD4_C7_OX40;CD8_C6_CD39;Mac_C2_SPP1;Endo_C3_RGCC;FB_C3_COL1A1"
all_long$cell_type = all_long$cell
all_long$cell_type[which(all_long$cell!=mycell)]<-"other"
all_long$cell_type[which(all_long$cell==mycell)]<-"mimer"

all = rbind(all_long,all_long_other)
all$cell_type  = factor(all$cell_type,
                             levels = c("mimer",'other','CD4','CD8','Mac','Endo','FB'))
mycolor<-c("#1B9E77","grey","#1e90ff","#20b2aa",'#ba6338',"#ce3d32","#f0e685")
my_comparison<-list(c('mimer','other'))


all_cell_percentage<-as.data.frame(all_cell_percentage)
rownames(all_cell_percentage)<-lapply(strsplit(rownames(all_cell_percentage), '[.]'),function(x){res<-x[1];return(res);})%>%unlist()
ggplot(as.data.frame(all_cell_percentage), aes(x = mimer_percentage))+ geom_density(color = "black", fill = "gray")
drop<-rownames(all_cell_percentage)[all_cell_percentage$mimer_percentage<0.05]
all_sub<-all[all$sample%in%drop==F,]


p<-ggplot(all_sub, aes(x=cell_type, y=closeness, fill = cell_type)) + 
  geom_boxplot(width=0.5, outlier.shape = NA, alpha = 0.8)+ 
  scale_fill_manual(values = mycolor)+ stat_compare_means(paired = FALSE, comparisons = my_comparison) + 
  theme(axis.text = element_text(color = "black",size = 12),
        axis.title = element_text(size = 14)) + 
  theme_classic() + 
  NoLegend()

p+coord_fixed(ratio = 0.5)

# Figure 4b -----
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


# Figure 4c -----
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
    
# Figure 4f -----
df<-data.frame('Score Type'=c('Immunosuppression','Angiogenesis','Invasion and metastasis','Tumor suppression','Phagocytosis','Inflammatory macrophage','Anti-inflammatory macrophage','Exhaustion','Toxicity','MHCI score'))
df$gene<-'';df[1,]$gene<-paste('CD274','PDCD1LG2','ICOSLG','CD276','VTCN1','VSIR','IDO1','TGFB1','IL10','CCL5','CCL17','CCL22','CXCL8',
                               'CCL16','CCL18','IL17A','CCL24','IL1R1','IL4','IL13','ARG1','CD40LG','TPH1','XBP1','SOCS1',sep = ",")
df[2,]$gene<-paste('VEGFA','CXCL8','CXCL1','CXCL2','TEK','ANGPT2','CXCL12','FLT1','CXCR4','CTSB','CTSS','SEMA4D','PROK2',
                   'CXCR1','S1PR1','PDGFA','FGF2','MMP12','VEGFC','DNTTIP2','PGF','HIF1A','NFE2L2','CCND2','CCNE1','CD44',
                   'E2F3','EDN1','EZH2','FGF18','FGFR1','FYN','HEY1','ITGAV','JAG1','JAG2','NOTCH1','PTK2','SPP1','STC1',
                   'TNFAIP6','TYMP','VAV2',sep = ",")
df[3,]$gene<-paste('MMP7','MMP2','MMP3','MMP9','CCL2','CCL3','CCL22','TGFB1','EGF','CCR2','FLT1','CTSS','CTSB','CXCR3',
                   'WNT7B','ALOX5','CDH1','CCL18','CCL24','S1PR1','CXCL16','MAPK7','CSF1','SPARC','TLR4','VCAM1','CCL20',sep = ",")
df[4,]$gene<-paste('IL1B','IL6','IL12A','IL23A','TNF','CXCL9','CXCL10','TLR2','TLR4','CXCL11','IFNG','CD40','FCGR2A',
                   'ITGAX','IFNGR1','HLA-DPB1','HLA-DPA1','HLA-DRA','HLA-DRB1','HLA-DQA1','CD74','HLA-DRB5','IRF5',sep = ",")
df[5,]$gene<-paste('MRC1','CD163','MERTK','C1QB','FCRLA','CD5L','CD81','GPNMB','CD36',sep = ",")
df[6,]$gene<-paste('IL23A','TNF','CD86','IL1A','IL1B','IL6','CCL5','IRF5','IRF1','CD40','CD74','HLA-DPA1','HLA-DPB1',
                   'HLA-DRB5','HLA-DRB1','CD83','CD68','CD80','S100A8','S100A9',sep = ",")
df[7,]$gene<-paste('CCL20','CCL18','CCL22','LYVE1','VEGFA','CTSB','CTSD','TGFB1','MMP19','MMP9','CLEC7A','WNT7B','FASLG',
                   'TNFSF12','TNFSF8','CD276','MSR1','FN1','IRF4','MRC1','CD163','ARG1','IL10','MARCO',sep = ",")
df[8,]$gene<-paste("LAYN","TIGIT","CTLA4","PDCD1","HAVCR2","ENTPD1","EOMES","TOX",sep = ",")
df[9,]$gene<-paste("PRF1","IFNG","GNLY","NKG7","GZMB","GZMH","FGFBP2","CCL4",sep = ",")
df[10,]$gene<-paste("HLA-A","HLA-B","HLA-C","TAP1","TAP2","B2M",sep = ",")

head(df)

marker_list = list(
  exhaustion = c("LAYN","TIGIT","CTLA4","PDCD1","HAVCR2","ENTPD1","EOMES","TOX"),
  pro_fibrotic_signature = c( "CCL22", "CSF1", "CHIT1", "FCMR", "SIGLEC15", "CTSK", "COL22A1", "CHI3L1",
                              "SPP1", "SPARC", "MMP9", "MMP7", "GPC4", "TIMP3", "MATK", "LIPA", "PALLD", "FDX1", "LPL"),
  ECM = ECM ,
  Cytolytics_effector = c("PRF1","IFNG","GNLY","NKG7","GZMB","GZMH","FGFBP2","CCL4"),
  
  Phagocytosis = c('MRC1','CD163','MERTK','C1QB','FCRLA','CD5L','CD81','GPNMB','CD36'),
  Anti_inflammatory = c('CCL20','CCL18','CCL22','LYVE1','VEGFA','CTSB','CTSD','TGFB1','MMP19','MMP9','CLEC7A','WNT7B','FASLG',
                        'TNFSF12','TNFSF8','CD276','MSR1','FN1','IRF4','MRC1','CD163','ARG1','IL10','MARCO'),
  
  Immunosuppression = c('CD274','PDCD1LG2','ICOSLG','CD276','VTCN1','VSIR','IDO1','TGFB1','IL10','CCL5','CCL17','CCL22','CXCL8',
                        'CCL16','CCL18','IL17A','CCL24','IL1R1','IL4','IL13','ARG1','CD40LG','TPH1','XBP1','SOCS1'),
  Angiogenesis = c('VEGFA','CXCL8','CXCL1','CXCL2','TEK','ANGPT2','CXCL12','FLT1','CXCR4','CTSB','CTSS','SEMA4D','PROK2',
                   'CXCR1','S1PR1','PDGFA','FGF2','MMP12','VEGFC','DNTTIP2','PGF','HIF1A','NFE2L2','CCND2','CCNE1','CD44',
                   'E2F3','EDN1','EZH2','FGF18','FGFR1','FYN','HEY1','ITGAV','JAG1','JAG2','NOTCH1','PTK2','SPP1','STC1',
                   'TNFAIP6','TYMP','VAV2'),
  Inflammatory = c('IL23A','TNF','CD86','IL1A','IL1B','IL6','CCL5','IRF5','IRF1','CD40','CD74','HLA-DPA1','HLA-DPB1',
                   'HLA-DRB5','HLA-DRB1','CD83','CD68','CD80','S100A8','S100A9'),
  Tumorsuppression = c('IL1B','IL6','IL12A','IL23A','TNF','CXCL9','CXCL10','TLR2','TLR4','CXCL11','IFNG','CD40','FCGR2A',
                       'ITGAX','IFNGR1','HLA-DPB1','HLA-DPA1','HLA-DRA','HLA-DRB1','HLA-DQA1','CD74','HLA-DRB5','IRF5'),
  pro_metastasis = c('MMP7','MMP2','MMP3','MMP9','CCL2','CCL3','CCL22','TGFB1','EGF','CCR2','FLT1','CTSS','CTSB','CXCR3',
                     'WNT7B','ALOX5','CDH1','CCL18','CCL24','S1PR1','CXCL16','MAPK7','CSF1','SPARC','TLR4','VCAM1','CCL20')
)
sp = AddModuleScore(sp,features = marker_list)
colnames(sp@meta.data)

sp@meta.data = dplyr::rename(sp@meta.data,exhaustion=Cluster1, pro_fibrotic_signature=Cluster2, ECM=Cluster3,
                             Cytolytics_effector=Cluster4, Phagocytosis=Cluster5,Anti_inflammatory=Cluster6 ,
                             Immunosuppression = Cluster7,Angiogenesis = Cluster8,Inflammatory = Cluster9,
                             Tumorsuppression = Cluster10,pro_metastasis = Cluster11)
fea <-names(marker_list)

fea = c("exhaustion", "pro_fibrotic_signature" ,"ECM" , "Immunosuppression" , "Angiogenesis", "Tumorsuppression", "pro_metastasis" )
dominated = subset(sp,development%in%"unknown",invert=T)
dominated = subset(dominated,new_mimer=="other",invert=T)
myd=as.data.frame(scale(dominated@meta.data[,fea]));
myd=cbind(myd,dominated@meta.data[,c('tissue','new_mimer')]);
head(myd)
# myd=myd[,-6];
myd=dplyr::group_by(myd, new_mimer) %>% dplyr::summarise(across(all_of(fea), mean))
myd=as.data.frame(myd);myd=tibble::column_to_rownames(myd,var='new_mimer'); myd=t(as.matrix(myd))
library(pheatmap)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
colors = c(colorRampPalette(colors = rev(getPalette(10))[1:5])(length(bk)/2),
           colorRampPalette(colors = rev(getPalette(10))[6:10])(length(bk)/2))
pdf("mimer_function_in_dominated.pdf",width = 8,height = 8)
# myd[myd< -0.5] = -2
pheatmap(myd, cellwidth=30, cellheight=30, cluster_cols=TRUE, clustering_method='ward.D',
         cluster_rows=F, treeheight_row=0, kmeans_k=NA, border_color="white", 
         scale="row",
         drop_levels=T,
         show_rownames=TRUE, show_colnames=TRUE, 
         color = colors,
         # color=rev(getPalette(10)),
         legend_breaks=seq(-2,2,0.5),breaks = bk,
         fontsize=12, fontsize_row=15, legend=T, fontsize_col=15 #, display_numbers=ifelse(myd > 0.37,"*","")
         )
dev.off()

ME = subset(sp,ME%in%"unknown",invert=T)
ME = subset(ME,new_mimer=="other",invert=T)
myd=as.data.frame(scale(ME@meta.data[,fea]));
myd=cbind(myd,ME@meta.data[,c('tissue','new_mimer')]);
head(myd)
# myd=myd[,-6];
myd=dplyr::group_by(myd, new_mimer) %>% dplyr::summarise(across(all_of(fea), mean))
myd=as.data.frame(myd);myd=tibble::column_to_rownames(myd,var='new_mimer'); myd=t(as.matrix(myd))
library(pheatmap)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))
pdf("mimer_function_in_ME.pdf",width = 8,height = 8)
pheatmap(myd, cellwidth=30, cellheight=30, cluster_cols=TRUE, 
         clustering_method='ward.D',
         cluster_rows=F, treeheight_row=0, kmeans_k=NA, border_color="white",
         scale="row",
         drop_levels=T,
         show_rownames=TRUE, show_colnames=TRUE, 
         color = colors,
         legend_breaks=seq(-2,2,0.5),breaks = bk,
         # color=rev(getPalette(10)),
         fontsize=12, fontsize_row=15, legend=T, fontsize_col=15 #, display_numbers=ifelse(myd > 0.37,"*","")
)
dev.off()
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

# Figure 4h -----
#up
load('sp_correct.RData')
df<-read.table("sp_table.txt", sep='\t',header = T)
gene<-c("CXCL13","ACP5","LAG3","PHLDA1","HAVCR2","RGS2","FOXP3","TIGIT","BATF","TNFRSF18","TNFRSF4","TNFRSF9","COL1A1","COL3A1","COL1A2","SPARC","FN1","POSTN","PLVAP","COL4A1","COL4A2","HSPG2","VWF","IGFBP7","SPP1","APOC1","MMP12","MMP9","FBP1","APOE")
Idents(sp)<-sp$development;
p<-subset(sp,development=='unknown',invert=T)
Idents(p)<-p$development
p<-AddModuleScore(p,features = list(gene),name = 'MIMER')
length(colnames(p@meta.data))
colnames(p@meta.data)[84]<-'MIMER'
my_levels<-c('Nor','Hyp','MiD','MoD','SD&CA','ICA','MCA')
Idents(p)<-factor(Idents(p),levels = my_levels)
mycolor<-c("#1B9E77","#D95F02","#666666","#E7298A","#66A61E","#E6AB02","#A6761D")
my_comparison<-list(c('Nor','Hyp'),c('Hyp','MiD'),c('MiD','MoD'),c('MoD','SD&CA'),c('SD&CA','ICA'),c('ICA','MCA'))

VlnPlot(p,features = c('MIMER'),pt.size = 0,assay = 'spatial',cols = mycolor)+
  geom_signif(comparisons = my_comparison,y_position = 1.4,textsize=4,step_increase=0)+NoLegend()+
  theme(axis.title.x = element_blank())+ggtitle("")+ylab('MIMER expression')
#bottom
Idents(sp)<-sp$development;
p<-subset(sp,development=='unknown',invert=T)
Idents(p)<-p$development
p<-AddModuleScore(p,features = list(gene),name = 'MIMER')
length(colnames(p@meta.data))
colnames(p@meta.data)[84]<-'MIMER'
my_levels<-c('Nor','Hyp','MiD','MoD','SD&CA','ICA','MCA')
Idents(p)<-factor(Idents(p),levels = my_levels)
mycolor<-c("#1B9E77","#D95F02","#666666","#E7298A","#66A61E","#E6AB02","#A6761D")
my_comparison<-list(c('Nor','Hyp'),c('Hyp','MiD'),c('MiD','MoD'),c('MoD','SD&CA'),c('SD&CA','ICA'),c('ICA','MCA'))

VlnPlot(p,features = c('MIMER'),pt.size = 0,assay = 'spatial',cols = mycolor)+
  geom_signif(comparisons = my_comparison,y_position = 1.4,textsize=4,step_increase=0)+NoLegend()+
  theme(axis.title.x = element_blank())+ggtitle("")+ylab('MIMER expression')
# Figure 4i
library(survival)
library(survminer)
pdf(file="hazard_ratio.pdf", width=8, height=6)
data <- read.table("data1.txt",header=T)

model <- coxph(Surv(Survival_time, Survival_status) ~ Age_2 + Grade_2 + AJCC_2 + MIMER, data = data)
ggforest(model, data=data, main='Hazard ratio', refLabel = '1', noDigits =3)

fit <- survfit(Surv(Survival_time, Survival_status) ~MIMER, data=data)
ggsurvplot(fit, data = data, conf.int = TRUE, pval = TRUE, fun = "pct", xlim = c(0, 80), risk.table = TRUE, size = 1, linetype = "strata", palette = c("#E64B35CC", "#4DBBD5CC"), legend = "top",legend.labs = c("MIMER high", "MIMER low"), legend.title = "")


# Figure 4j-----
total_le = rbind(hyp_LE,SDCA_LE)
total_me = rbind(hyp_ME,SDCA_ME)
total_me$distance = gsub("far","Distal",
                         gsub("near","Proximal",
                              total_me$distance))

ggplot(total_me,aes(x=development,y=MIMER,
                          color=distance,
                          fill=distance))+
  geom_violin()+
  stat_compare_means(mapping = aes(condition='distance'),method = "wilcox.test")+NoLegend()+theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size=16,color = "black"),
        axis.title.y  = element_text(size=18,color = "black"),
        legend.text  = element_text(size = 12),
        legend.title =  element_text(size = 14))+ggtitle("")+
  ylab('ME\nMIMER expression')

total_major$distance = gsub("far","Distal",
                            gsub("near","Proximal",
                                 total_major$distance))
ggplot(total_major,aes(x=development,y=MIMER,
                             color=distance,
                             fill=distance))+
  geom_violin()+
  stat_compare_means(mapping = aes(condition='distance'),method = "wilcox.test")+NoLegend()+theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size=16,color = "black"),
        axis.title.y  = element_text(size=18,color = "black"),
        legend.text  = element_text(size = 12),
        legend.title =  element_text(size = 14))+ggtitle("")
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
