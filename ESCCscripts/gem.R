setwd("/workspace/wangxiliang/project/escc_baseline_zhanQiMin/")

library(dplyr,lib.loc='/workspace/zhangzhichao/miniconda3/envs/scRNA/lib/R/library')
library(ggplot2,lib.loc='/workspace/zhangzhichao/miniconda3/envs/scRNA/lib/R/library')
library(ggpubr,lib.loc='/workspace/zhangzhichao/miniconda3/envs/scRNA/lib/R/library')
library(qs,lib.loc='/workspace/zhangzhichao/miniconda3/envs/scRNA/lib/R/library')
library(Matrix,lib.loc='/workspace/zhangzhichao/miniconda3/envs/scRNA/lib/R/library'); library(doParallel)
library(Seurat,lib.loc='/workspace/zhangzhichao/miniconda3/envs/scRNA/lib/R/library')
library(patchwork,lib.loc='/workspace/zhangzhichao/miniconda3/envs/scRNA/lib/R/library')

#### GEM in PD1+ Tex, >= 10 spots
sp = qs::qread('input/sp_correct_with_rctd_20250523.qs'); DefaultAssay(sp) = 'SCT'
mys = c('B_JT1','D_JT2','IT2','IT3','KT1_1','KT1_2','KT2','MLN_N','MLN_P','MT3','PT3LN_P_2','ST1','TT2_2','WT1',
        'XT2_1','XT2_2','YLN_3LN_P','YT1','ZT3')
res = data.frame()
for(i in 1:length(mys)){
  mc = read.csv('input/cellID_mimer_each_celltype_enrich_mark_for_wxl.csv',header=T)
  mc = mc[mc$CD8_C6_CD39_above_median %in% 'TRUE' & mc$file %in% mys[i],]
  sc = subset(sp, cells=mc$cellID); sc@meta.data$class = 'With'
  sc@meta.data[ rownames(sc@meta.data) %in% mc[mc$mimer %in% 'Others','cellID'],'class'] = 'WithOut'
  md = FetchData(sc, vars=c('GEM','class'))
  tmp = aggregate(x=md$GEM, by=list(md$class), FUN=mean); tmp$file = mys[i]; res = rbind(res,tmp) }
colnames(res)[c(1,2)] = c('class','val'); write.table(res,file='gem_Tex.txt',quote=F,sep='\t',row.names=F,col.names=T)

ggpaired(data=res, x="class", y="val", color="class", line.color="gray", line.size=0.4, palette="jco", id='file',
         xlab=F, point.size=3, ylab='GEM') + theme(legend.position='none') +
  stat_compare_means(comparisons=list(c('With','WithOut')),method='wilcox.test',label="p.format",size=5,paired=TRUE)
ggsave(filename='gem_Tex.pdf', device='pdf', dpi=300, width=6, height=6)

#### SPP1+ TAM, >= 10 spots
sp = qs::qread('input/sp_correct_with_rctd_20250523.qs'); DefaultAssay(sp) = 'SCT'
mys = c('B_JT1','D_JT2','E3_C','IT1','IT2','IT3','KT1_2','KT2','MLN_N','MLN_P','MT3','MT3_2','NT2','NT2_2','PT1',
        'PT3LN_P_2','QP','SLN_1','ST1','TT2','TT2_2','WT1','XT2_1','XT2_2','YLN_3LN_P','YT1','YT2','ZT2','ZT3','ZLN_P')
res = data.frame()
for(i in 1:length(mys)){
  mc = read.csv('input/cellID_mimer_each_celltype_enrich_mark_for_wxl.csv',header=T)
  mc = mc[mc$Mac_C2_SPP1_above_median %in% 'TRUE' & mc$file %in% mys[i],]
  sc = subset(sp, cells=mc$cellID); sc@meta.data$class = 'With'
  sc@meta.data[ rownames(sc@meta.data) %in% mc[mc$mimer %in% 'Others','cellID'],'class'] = 'WithOut'
  md = FetchData(sc, vars=c('CALR','APP','IL15','IL10','IGF1','class'))
  tmp = aggregate(x=md[,c(1:5)], by=list(md$class), FUN=mean); tmp$file = mys[i]; res = rbind(res,tmp) }
colnames(res)[1]='class'; write.table(res,file='genes_SPP1.txt',quote=F,sep='\t',row.names=F,col.names=T)

for(i in c('CALR','APP','IL15','IL10','IGF1')){
  ggpaired(data=res, x="class", y=i, color="class", line.color="gray", line.size=0.4, palette="jco", id='file',
           xlab=F, point.size=3, ylab=i) + theme(legend.position='none') +
    stat_compare_means(comparisons=list(c('With','WithOut')), method='wilcox.test', label="p.format", size=5, paired=TRUE)
  ggsave(filename=paste0(i,'_SPP1.pdf'), device='pdf', dpi=300, width=6, height=6) }

#### RGCC+ End, >= 10 spots
sp = qs::qread('input/sp_correct_with_rctd_20250523.qs'); DefaultAssay(sp) = 'SCT'
mys = c('B_JT1','D_JT2','E2_D','E3_C','IT1','IT2','IT3','KNP','KT1_1','KT1_2','KT2','MLN_N','MLN_P','MT3','MT3_2',
        'NT2','NT2_2','PT1','PT1_2','PT3LN_P_2','QP','QT1','SLN_1','ST1','TT2','TT2_2','WLn_N','WT1','XLn_1Ln_2',
        'XT2_1','XT2_2','YLN_3LN_P','YP','YT1','YT2','ZLN_N','ZT2','ZT3','ZLN_P')
res = data.frame()
for(i in 1:length(mys)){
  mc = read.csv('input/cellID_mimer_each_celltype_enrich_mark_for_wxl.csv',header=T)
  mc = mc[mc$Endo_C3_RGCC_above_median %in% 'TRUE' & mc$file %in% mys[i],]
  sc = subset(sp, cells=mc$cellID); sc@meta.data$class = 'With'
  sc@meta.data[ rownames(sc@meta.data) %in% mc[mc$mimer %in% 'Others','cellID'],'class'] = 'WithOut'
  md = FetchData(sc, vars=c('CALR','APP','IL15','IL10','IGF1','class'))
  tmp = aggregate(x=md[,c(1:5)], by=list(md$class), FUN=mean); tmp$file = mys[i]; res = rbind(res,tmp) }
colnames(res)[1]='class'; write.table(res,file='genes_RGCC.txt',quote=F,sep='\t',row.names=F,col.names=T)

for(i in c('CALR','APP','IL15','IL10','IGF1')){
ggpaired(data=res, x="class", y=i, color="class", line.color="gray", line.size=0.4, palette="jco", id='file',
         xlab=F, point.size=3, ylab=i) + theme(legend.position='none') +
  stat_compare_means(comparisons=list(c('With','WithOut')), method='wilcox.test', label="p.format", size=5, paired=TRUE)
ggsave(filename=paste0(i,'_RGCC.pdf'), device='pdf', dpi=300, width=6, height=6) }

#### OX40+ Treg, >= 10 spots
sp = qs::qread('input/sp_correct_with_rctd_20250523.qs'); DefaultAssay(sp) = 'SCT'
mys = c('B_JT1','D_JT2','E3_C','IT1','IT2','IT3','KT1_1','KT1_2','KT2','MLN_N','MLN_P','MT3','MT3_2','NT2','NT2_2',
        'PT1','PT1_2','PT3LN_P_2','QP','SLN_1','ST1','TT2','TT2_2','WLn_N','WT1','XLn_1Ln_2','XT2_1','XT2_2',
        'YLN_3LN_P','YP','YT1','YT2','ZLN_N','ZT3','ZLN_P')
res = data.frame()
for(i in 1:length(mys)){
  mc = read.csv('input/cellID_mimer_each_celltype_enrich_mark_for_wxl.csv',header=T)
  mc = mc[mc$CD4_C7_OX40_above_median %in% 'TRUE' & mc$file %in% mys[i],]
  sc = subset(sp, cells=mc$cellID); sc@meta.data$class = 'With'
  sc@meta.data[ rownames(sc@meta.data) %in% mc[mc$mimer %in% 'Others','cellID'],'class'] = 'WithOut'
  md = FetchData(sc, vars=c('CALR','APP','IL15','IL10','IGF1','class'))
  tmp = aggregate(x=md[,c(1:5)], by=list(md$class), FUN=mean); tmp$file = mys[i]; res = rbind(res,tmp) }
colnames(res)[1]='class'; write.table(res,file='genes_OX40.txt',quote=F,sep='\t',row.names=F,col.names=T)

for(i in c('CALR','APP','IL15','IL10','IGF1')){
  ggpaired(data=res, x="class", y=i, color="class", line.color="gray", line.size=0.4, palette="jco", id='file',
           xlab=F, point.size=3, ylab=i) + theme(legend.position='none') +
    stat_compare_means(comparisons=list(c('With','WithOut')), method='wilcox.test', label="p.format", size=5, paired=TRUE)
  ggsave(filename=paste0(i,'_OX40.pdf'), device='pdf', dpi=300, width=6, height=6) }


#### COL1A1+ CAF, >= 10 spots
sp = qs::qread('input/sp_correct_with_rctd_20250523.qs'); DefaultAssay(sp) = 'SCT'
mys = c('B_JT1','C_JT2','D_JT2','E2_D','E3_C','IT1','IT2','IT3','KT1_1','KT1_2','KT2','MLN_P','MT3','MT3_2','NT2',
        'NT2_2','PT1','PT1_2','PT3LN_P_2','QP','QT1','QT1_2','SLN_1','ST1','TT2','TT2_2','WT1','XT2_1','XT2_2',
        'YLN_3LN_P','YT1','YT2','ZT2','ZT3','ZLN_P')
res = data.frame()
for(i in 1:length(mys)){
  mc = read.csv('input/cellID_mimer_each_celltype_enrich_mark_for_wxl.csv',header=T)
  mc = mc[mc$FB_C3_COL1A1_above_median %in% 'TRUE' & mc$file %in% mys[i],]
  sc = subset(sp, cells=mc$cellID); sc@meta.data$class = 'With'
  sc@meta.data[ rownames(sc@meta.data) %in% mc[mc$mimer %in% 'Others','cellID'],'class'] = 'WithOut'
  md = FetchData(sc, vars=c('CALR','APP','IL15','IL10','IGF1','class'))
  tmp = aggregate(x=md[,c(1:5)], by=list(md$class), FUN=mean); tmp$file = mys[i]; res = rbind(res,tmp) }
colnames(res)[1]='class'; write.table(res,file='genes_COL1A1.txt',quote=F,sep='\t',row.names=F,col.names=T)

for(i in c('CALR','APP','IL15','IL10','IGF1')){
  ggpaired(data=res, x="class", y=i, color="class", line.color="gray", line.size=0.4, palette="jco", id='file',
           xlab=F, point.size=3, ylab=i) + theme(legend.position='none') +
    stat_compare_means(comparisons=list(c('With','WithOut')), method='wilcox.test', label="p.format", size=5, paired=TRUE)
  ggsave(filename=paste0(i,'_COL1A1.pdf'), device='pdf', dpi=300, width=6, height=6) }
##