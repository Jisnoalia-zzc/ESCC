# E7 c-d
library(ggplot2)
library(scales)
library(data.table)
library(dplyr)
library(tidyr)
library(patchwork)
library(cowplot)
setwd('C:/Users/Ron/Desktop/Figs/E7/')
mimer_median_enrich_development_ME <- fread('C:/Users/Ron/Desktop/Figs/E7/E7_c-g_mimer_median_enrich_development_ME.csv')

# 设置目标细胞类型列（使用原始比例值）
cell_prop_columns <- c("CD8.C6.CD39", "CD4.C7.OX40", "Mac.C2.SPP1", 
                       "FB.C3.COL1A1", "Endo.C3.RGCC")

# 设置数据表格式
setDT(mimer_median_enrich_development_ME)

# 将相关列转换为字符类型以避免因子比较问题
mimer_median_enrich_development_ME[, ME := as.character(ME)]
mimer_median_enrich_development_ME[, development := as.character(development)]

# 定义阶段顺序
stage_order <- c("Nor", "Hyp", "MiD", "MoD", "SD&CA", "ICA")
# 定义完整的ME和Development阶段（包括后缀）
me_levels <- paste0(stage_order, "-ME")
dev_levels <- stage_order

# 定义分类颜色
class_colors <- c("ME" = "#E41A1C", "Development" = "#377EB8")

# 存储单个图形的列表
plot_list <- list()

# --- Helper Function for Statistics and Comparisons (修改为处理单个细胞类型连续值) ---
calculate_stats_and_comparisons_cell_prop <- function(data_subset, stage_col, levels_to_use, column_name, class_name) {
  # Filter data based on the specified column and levels
  filtered_data <- data_subset[get(stage_col) %in% levels_to_use]
  filtered_data[, type := get(stage_col)]
  filtered_data$type <- factor(filtered_data$type, levels = levels_to_use)
  
  # 检查目标列是否存在
  if(!column_name %in% names(filtered_data)) {
    stop(paste("Column", column_name, "not found in data"))
  }
  
  # 按文件和阶段计算目标细胞的比例均值（使用连续值）
  file_means <- filtered_data[
    , .(
      mean_value = mean(get(column_name), na.rm = TRUE)
    ), 
    by = .(file, type)
  ]
  
  # 计算每个阶段的汇总统计量
  mean_stats <- file_means[
    , .(
      mean_prop = mean(mean_value, na.rm = TRUE),
      se_prop = ifelse(.N > 1, sd(mean_value, na.rm = TRUE) / sqrt(.N), 0),
      n_samples = .N
    ), 
    by = type
  ]
  
  mean_stats$class <- class_name
  mean_stats$column <- column_name
  
  # 获取Nor阶段的数据用于统计检验
  if(class_name == "ME") {
    nor_stage <- "Nor-ME"
  } else {
    nor_stage <- "Nor"
  }
  
  # 提取Nor阶段的均值数据
  nor_data <- file_means[type == nor_stage, mean_value]
  
  # 统计检验结果 - 每个阶段与Nor比较
  comp_results <- list()
  
  # 对除了Nor以外的每个阶段进行与Nor的比较
  for(stage in levels_to_use) {
    if(stage != nor_stage) {
      # 提取当前阶段的均值数据
      stage_data <- file_means[type == stage, mean_value]
      
      n_nor <- length(nor_data)
      n_stage <- length(stage_data)
      
      # 只有当两个组都有至少2个样本时才进行检验
      if (n_nor >= 2 && n_stage >= 2) {
        comp_test <- wilcox.test(nor_data, stage_data, exact = FALSE)
        p_val <- comp_test$p.value
      } else {
        p_val <- NA
      }
      
      comp_results[[stage]] <- list(
        p.value = p_val, 
        n_nor = n_nor, 
        n_stage = n_stage,
        stage_name = stage,
        base_stage = ifelse(class_name == "ME", gsub("-ME$", "", stage), stage)
      )
    }
  }
  
  return(list(
    stats = mean_stats,
    comparisons = comp_results,
    file_means = file_means  # 保存文件级均值数据
  ))
}

# --- Helper Function for creating x-axis labels ---
create_x_labels <- function(me_stats, dev_stats, stage_order, me_levels, dev_levels) {
  x_labels <- character(length(stage_order))
  
  for (i in seq_along(stage_order)) {
    stage <- stage_order[i]
    
    # 获取ME分类的样本数
    me_stage_name <- paste0(stage, "-ME")
    me_n <- if (me_stage_name %in% me_stats$type) {
      me_stats[type == me_stage_name, n_samples]
    } else {
      0
    }
    
    # 获取Development分类的样本数
    dev_n <- if (stage %in% dev_stats$type) {
      dev_stats[type == stage, n_samples]
    } else {
      0
    }
    
    # 创建标签，确保me_n和dev_n是单个数值
    if (length(me_n) == 1 && length(dev_n) == 1) {
      x_labels[i] <- sprintf("%s\nME:n=%d\nDev:n=%d", stage, me_n, dev_n)
    } else {
      x_labels[i] <- stage  # 如果获取不到样本数，只显示阶段名称
    }
  }
  
  return(x_labels)
}

# --- Main Loop over Columns ---
for (column_name in cell_prop_columns) {
  # 创建友好的列名用于图形标题
  friendly_name <- column_name
  cat(sprintf("\nProcessing column: %s\n", friendly_name))
  
  # Process ME data
  me_result <- calculate_stats_and_comparisons_cell_prop(
    mimer_median_enrich_development_ME, "ME", me_levels, column_name, "ME"
  )
  me_stats <- me_result$stats
  me_comps <- me_result$comparisons
  me_file_means <- me_result$file_means
  
  # Process Development data
  dev_result <- calculate_stats_and_comparisons_cell_prop(
    mimer_median_enrich_development_ME, "development", dev_levels, column_name, "Development"
  )
  dev_stats <- dev_result$stats
  dev_comps <- dev_result$comparisons
  dev_file_means <- dev_result$file_means
  
  # Combine stats for plotting
  combined_stats <- rbindlist(list(me_stats, dev_stats))
  # Map back to base stage names for x-axis consistency
  combined_stats[, base_stage := gsub("-ME$", "", type)] # Remove -ME suffix
  combined_stats$base_stage <- factor(combined_stats$base_stage, levels = stage_order)
  
  # 创建包含样本数的x轴标签
  x_labels <- create_x_labels(me_stats, dev_stats, stage_order, me_levels, dev_levels)
  names(x_labels) <- stage_order
  
  # Prepare annotation data for p-values
  annotation_list <- list()
  
  # Add ME annotations (每个阶段与Nor比较)
  for(stage_name in names(me_comps)) {
    comp <- me_comps[[stage_name]]
    if(!is.na(comp$p.value)) {
      # 找到该阶段在stage_order中的位置
      stage_pos_in_order <- which(stage_order == comp$base_stage)
      # 找到该阶段ME分类的y最大值
      y_max_vals <- combined_stats[class == "ME" & base_stage == comp$base_stage, mean_prop]
      y_pos <- if(length(y_max_vals) > 0) max(y_max_vals, na.rm = TRUE) * 1.1 else 0.01
      
      annotation_list[[paste0("ME_", stage_name)]] <- data.table(
        x_pos_index = stage_pos_in_order,
        y_pos = y_pos,
        p_value = comp$p.value,
        class = "ME",
        comparison = paste0("vs Nor (n=", comp$n_nor, " vs ", comp$n_stage, ")"),
        stage = comp$base_stage
      )
    }
  }
  
  # Add Development annotations (每个阶段与Nor比较)
  for(stage_name in names(dev_comps)) {
    comp <- dev_comps[[stage_name]]
    if(!is.na(comp$p.value)) {
      # 找到该阶段在stage_order中的位置
      stage_pos_in_order <- which(stage_order == comp$base_stage)
      # 找到该阶段Development分类的y最大值
      y_max_vals <- combined_stats[class == "Development" & base_stage == comp$base_stage, mean_prop]
      y_pos <- if(length(y_max_vals) > 0) max(y_max_vals, na.rm = TRUE) * 1.1 else 0.01
      
      annotation_list[[paste0("Dev_", stage_name)]] <- data.table(
        x_pos_index = stage_pos_in_order,
        y_pos = y_pos,
        p_value = comp$p.value,
        class = "Development",
        comparison = paste0("vs Nor (n=", comp$n_nor, " vs ", comp$n_stage, ")"),
        stage = comp$base_stage
      )
    }
  }
  
  # 检查是否有注释数据
  if (length(annotation_list) > 0) {
    annotation_data <- rbindlist(annotation_list)
    annotation_data[, significance := ifelse(p_value < 0.001, "***",
                                             ifelse(p_value < 0.01, "**",
                                                    ifelse(p_value < 0.05, "*", "ns")))]
    annotation_data[, formatted_p := sapply(p_value, function(p) {
      if (is.na(p)) return("NA")
      if (p < 0.0001) formatC(p, format = "e", digits = 2)
      else if (p < 0.001) format(round(p, 4), nsmall = 4, scientific = FALSE)
      else if (p < 0.01) format(round(p, 3), nsmall = 3, scientific = FALSE)
      else if (p < 0.1) format(round(p, 2), nsmall = 2, scientific = FALSE)
      else format(round(p, 2), nsmall = 2, scientific = FALSE)
    })]
    annotation_data[, label := paste0("p=", formatted_p, " ", significance)]
  } else {
    annotation_data <- data.table()
  }
  
  # Create the plot for this cell type
  p <- ggplot(combined_stats, aes(x = base_stage, y = mean_prop, group = class, color = class)) +
    geom_line(linewidth = 1, aes(linetype = class)) +
    geom_point(size = 2.5, position = position_dodge(width = 0.1)) +
    geom_errorbar(aes(ymin = mean_prop - se_prop, ymax = mean_prop + se_prop),
                  width = 0.1, linewidth = 0.8, position = position_dodge(width = 0.1)) +
    scale_color_manual(values = class_colors, name = "Classification") +
    scale_linetype_manual(values = c("ME" = "solid", "Development" = "dashed"), 
                          name = "Classification") +
    scale_x_discrete(limits = stage_order, labels = x_labels) +
    labs(
      title = friendly_name,
      x = "Disease Stage", 
      y = "Mean Proportion"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, lineheight = 0.8),
      axis.title = element_text(face = "bold", size = 11),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.position = "top",
      legend.key.width = unit(1.5, "lines"),
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
    )
  
  # 如果有p值注释数据，添加到图形中
  if (nrow(annotation_data) > 0) {
    # 为ME和Development分别设置不同的x偏移量，避免重叠
    annotation_data[, x_offset := ifelse(class == "ME", -0.2, 0.2)]
    annotation_data[, x_pos := x_pos_index + x_offset]
    
    # 调整y位置，避免重叠
    annotation_data[, y_offset := ifelse(class == "ME", 0.00005, 0.0001)]
    annotation_data[, y_pos_adj := y_pos + y_offset]
    
    p <- p + geom_text(
      data = annotation_data,
      aes(x = x_pos, y = y_pos_adj, label = label, color = class),
      inherit.aes = FALSE,
      size = 2.8, fontface = "bold", vjust = -0.5, show.legend = FALSE
    ) +
      annotate("text", x = 1, y = max(combined_stats$mean_prop) * 1.2, 
               label = "Nor (Reference)", size = 3, color = "gray40")
  }
  
  # Store the plot
  plot_list[[friendly_name]] <- p
}

# 检查是否有图形生成
if (length(plot_list) > 0) {
  # 使用patchwork合并图形
  final_plot <- wrap_plots(plot_list, ncol = 3) +
    plot_annotation(
      title = "MIMER Cell Type Proportions: ME vs Development Classification (Mean ± SEM)",
      subtitle = paste(
        "Analysis based on cell type proportion values",
        "Solid lines: ME classification | Dashed lines: Development classification",
        "All stages compared to Nor (reference group)",
        "* p<0.05, ** p<0.01, *** p<0.001 (Colored by classification type)",
        sep = "\n"
      ),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40", lineheight = 1.2)
      )
    )
  
  # 调整图形大小以适应新的x轴标签
  options(repr.plot.width = 18, repr.plot.height = 12)
  
  # 打印最终合并的图形
  print(final_plot)
} else {
  cat("No plots generated. Please check your data and target columns.\n")
}

#####################
library(psych)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(cowplot)
library(Seurat)
library(dplyr)
require("DESeq2")
# Extended Data Fig 7 -----
# Extended Data Fig 7a -----
####### calculate cell num ratio and abundance
load('scRNA.RData');
colnames(sc@meta.data)
table(sc$tissue)
primary_all = subset(sc,tissue=='Tu')
# normalize each sample by all cells in that sample, 
# in that dataset and condition
deseq2_norm_single = function(data){
  
  a.mtx = data %>% 
    dplyr::select(sample_id, sub, sub_num) %>%
    pivot_wider(.,names_from = sub, values_from = sub_num) %>%
    replace(., is.na(.),0) %>%
    column_to_rownames(var='sample_id')     %>%
    as.matrix() %>%
    t()
  
  sf <- estimateSizeFactorsForMatrix(a.mtx, type ="poscounts")
  a.norm <- log2(sweep(a.mtx,2,sf,"/")+1)
  
  a.normdf = a.norm %>% 
    t() %>%
    as.data.frame() %>%
    as_tibble(.,rownames='sample_id') %>%
    pivot_longer(!sample_id, names_to = 'sub',values_to = 'abundance')%>%
    mutate(`sf.DESeq` = sf[sample_id])
  return(a.normdf)
} 

load('scRNA_cd4.RData');
my@meta.data %>% select(c(subC)) -> CD4
load('scRNA_cd8.RData');
my@meta.data %>% select(c(subC)) -> CD8
load('scRNA_b.RData');
my@meta.data %>% select(c(subC)) -> Bcells
load('scRNA_endo.RData');
my@meta.data %>% select(c(subC)) -> Endo
load('scRNA_fb.RData');
my@meta.data %>% select(c(subC)) -> Fibro
load('scRNA_myeloid.RData');
my@meta.data %>% select(c(subC)) -> Myeloid
CD4$major = 'CD4'
CD8$major = 'CD8'
Bcells$major = "Bcell"
Endo$major = "Endo"
Fibro$major = "Fibro"
Myeloid$major = "Myeloid"
df = rbind(CD4,CD8,Bcells,Endo,Fibro,Myeloid)
df$cellID = rownames(df)
primary_all$cellID = rownames(primary_all@meta.data)
phe = primary_all@meta.data %>% select(c(cellID,loc,metastasis,tissue,patient))
df2 = df %>% left_join(phe,by = "cellID")
df2 = df2 %>% filter(tissue=="Tu")
rownames(df2) = df2$cellID
my$tissue = stringr::str_split(my$orig.ident,"_",simplify = T)[,2]
sub$cellID = rownames(sub@meta.data)
df2$newC = df2$subC

df2 = df2 %>% filter(newC!="Smooth muscle cell")

dat<-as.data.frame(table(df2$patient,df2$newC));
colnames(dat)<-c('sample_id','sub','sub_num');
dat<-deseq2_norm_single(dat)

# ratio -----
t<-as.data.frame(table(df2$patient,df2$newC));colnames(t)<-c('sample_id','sub','sub_num')
d<-as.data.frame(table(df2$newC,df2$major));colnames(d)<-c('sub','major','freq');
d<-d[d$freq>0,];t<-merge(t,d[,1:2],by='sub')
total<-as.data.frame(table(df2$patient,df2$major));
colnames(total)<-c('sample_id','major','major_num');t<-merge(t,total,by=c('sample_id','major'))
t$ratio<-t$sub_num/t$major_num
dat<-merge(dat,t[,c('sample_id','sub','sub_num','ratio')],by=c('sample_id','sub'))
############### Imm ALL R+NR+U ######################################
samples = dat %>%
    mutate(sub = factor(sub, levels=dat$sub %>% unique()))
samples[is.na(samples$ratio),'ratio']<-0

abun = samples %>% 
    dplyr::select(sub,sample_id,ratio  ) %>%
    pivot_wider(names_from = sub, values_from =ratio ) %>%
    column_to_rownames(var='sample_id')
# cor -----
res = corr.test(abun)
r1 = res$p 
r1[lower.tri(r1,diag = T)] = 0
r1 = r1 + t(r1)
res0 = r1 %>%
    as_tibble() %>%
    mutate_all(function(x)  case_when(
        (x>=0.01 & x<0.05)~'*', 
        (x>=0.001 & x<0.01)~'**',
        x<0.001~'***',
        x>0.05 ~' '))

n_CM = 7
res1 = pheatmap(res$r,clustering_distance_rows = 'correlation',
                clustering_distance_cols = 'correlation',
                clustering_method = 'ward.D2',
                color = colorRampPalette(rev(RColorBrewer::brewer.pal(11,'RdBu')))(100),
                breaks=seq(-1,1,length.out=101),silent = F,border_color = NA,
                # display_numbers =res0,
                # display_numbers = T,
                number_color = 'black',cutree_cols = n_CM,
                cutree_rows =n_CM )
res1

 

modules = data.frame('sub'=names(cutree(res1$tree_row,k=n_CM)),
                     'module'= cutree(res1$tree_row,k=n_CM)) 
modules = modules[res1$tree_row$labels[res1$tree_row$order],]
modules = modules %>%mutate(module=paste0('M',match(modules$module,unique(modules$module))))

rownames(modules)<-NULL;annotation_col = modules %>%
  distinct(sub,module) %>%
  column_to_rownames(var='sub') %>%
  mutate(module = factor(module))
ann_colors = list(
  module = c(M1="#FF6666",M2="#FFB266",M3="#FFCC99",M4="#66B266",M5="#66CCB2",M6="#6699CC",M7="#B266CC"))

pheatmap(res$r,clustering_distance_rows = 'correlation',annotation_col = annotation_col,annotation_colors = ann_colors,
         clustering_distance_cols = 'correlation',clustering_method = 'ward.D2',show_colnames = F,
         color = colorRampPalette(rev(RColorBrewer::brewer.pal(11,'RdBu')))(100),
         breaks=seq(-1,1,length.out=101),silent = F,border_color = NA,
         # display_numbers =res0,
         number_color = 'black',cutree_cols = n_CM)

# Extended Data Fig 7b -----
library("FactoMineR")
library("factoextra")
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(cluster)
library(ggbreak)
set.seed(123)
produc_pca <- PCA(scale(rall_m), ncp = 3, graph = FALSE)
produc_hcpc <- HCPC(produc_pca, graph = FALSE,method = "ward",min = 7)
fviz_cluster(produc_hcpc,
             repel = TRUE, 
             show.clust.cent = TRUE, labelsize=7,
             palette = "jco",  ggtheme = theme_minimal(),
             main = "")+ 
  theme_classic()+
  scale_x_break(c(4,6.5),
                space = 0.2,
                scales = 0.2)
# Extended Data Fig 7c -----
gene1<-c('CXCL13','ACP5','LAG3','PHLDA1','HAVCR2','RGS2','PLPP1','RHOB','SNX9','CCL5','CD8A','CD3D') # CD39+CD8
gene2<-c('FOXP3','TIGIT','BATF','TNFRSF18','TNFRSF4','TNFRSF9','IL32','CD4','IL10','IL2RA') #Treg
gene3<-c('SPP1','APOC1','MMP12','MMP9','FBP1','APOE','CTSB','CD68','CCL3','TYROBP') # SPP1+MAC
gene4<-c('PLVAP','COL4A1','COL4A2','HSPG2','VWF','IGFBP7','PECAM1','SERPINE1','SPARC','INSR') # Endo
gene5<-c('COL1A1','COL3A1','COL1A2','SPARC','FN1','POSTN','CST1','MMP11','CTHRC1','COL6A3') # CAF
gene6<-c('RNASE1','CCL18','C1QA','C1QB','C1QC','SELENOP','F13A1','PLTP','LGMN','LYVE1','CD68') # LYVE1+Mac
Idents(sp)<-sp$development;p<-subset(sp,ident=c('Nor','SD&CA','ICA'))
my_levels<-c('Nor','SD&CA','ICA')
Idents(p)<-factor(Idents(p),levels = my_levels)
p<-AddModuleScore(p,features = list(gene3),name = 'SPP1.TAM');p<-AddModuleScore(p,features = list(gene1),name = 'CD39.CD8')
p<-AddModuleScore(p,features = list(gene2),name = 'OX40.Treg');p<-AddModuleScore(p,features = list(gene4),name = 'PLVAP.Endo')
p<-AddModuleScore(p,features = list(gene5),name = 'POSTN.CAF');p<-AddModuleScore(p,features = list(gene6),name = 'LYVE1.Mac')
fea<-c('SPP1.TAM1','PDCD1.CD81','OX40.Treg1','PLVAP.Endo1','POSTN.CAF1','LYVE1.Mac1')
myd=as.data.frame(scale(p@meta.data[,fea]))
myd=cbind(myd,p@meta.data[,c('file','development','tissue','metastasis')])
colnames(myd)<-c('SPP1+TAM_score','PDCD1+CD8_score','OX40+Treg_score','PLVAP+Endo_score','POSTN+CAF_score','LYVE1+Mac_score','file','Development','tissue','metastasis')
s1<-myd[myd$Development=='Nor',];s1<-aggregate(s1[,1:6],by=list(Group=s1$file),mean);colnames(s1)[1]<-'file';s1$Development<-'Nor'
s2<-myd[myd$Development=='SD&CA',];s2<-aggregate(s2[,1:6],by=list(Group=s2$file),mean);colnames(s2)[1]<-'file';s2$Development<-'SD&CA'
s3<-myd[myd$Development=='ICA',];s3<-aggregate(s3[,1:6],by=list(Group=s3$file),mean);colnames(s3)[1]<-'file';s3$Development<-'ICA'
new<-rbind(s1,s2);new<-rbind(new,s3)
mycolor<-c("#1B9E77","#66A61E","#E6AB02")
my_comparison<-list(c('Nor','SD&CA'),c('Nor','ICA'));ggboxplot(new,x="Development",y="SPP1+TAM_score",color = "Development",palette = mycolor,add = "jitter",add.params = list(size=0.4))+
  stat_compare_means(comparisons = my_comparison)+NoLegend()