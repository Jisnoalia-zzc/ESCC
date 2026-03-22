# SupplementaryFig. 12 ------
# Supplementary Fig. 12 a -----
my<-subset(sc,L2_C=='T')
n<-subset(my,tissue=='nLN');p<-subset(my,tissue=='pLN')
t<-as.data.frame(table(n@meta.data$patient,n@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(n@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'nLN'
com<-t
t<-as.data.frame(table(p@meta.data$patient,p@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(p@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'pLN'
com<-rbind(com,t)
my_comparisons<-list(c('nLN','pLN')) 
sub<-com[com$Subcluster=='CD8_C6_PDCD1',]
colnames(sub)[5]<-"Ratio of CD8_C6_PDCD1 in T";sub$Tissue<-factor(sub$Tissue,levels = c('nLN','pLN'))
ggboxplot(sub,x="Tissue",y="Ratio of CD8_C6_PDCD1 in T",color = "Tissue",palette = c("#FC8D62","#8DA0CB"),add = "jitter",add.params = list(size=1.5))+
  stat_compare_means(method='anova',paired = FALSE)


my<-subset(sc,L2_C=='T')
n<-subset(my,tissue=='nLN');p<-subset(my,tissue=='pLN')
t<-as.data.frame(table(n@meta.data$patient,n@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(n@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'nLN'
com<-t
t<-as.data.frame(table(p@meta.data$patient,p@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(p@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'pLN'
com<-rbind(com,t)
my_comparisons<-list(c('nLN','pLN')) 
sub<-com[com$Subcluster=='CD4_C7_OX40',]
colnames(sub)[5]<-"Ratio of CD4_C7_OX40 in T";sub$Tissue<-factor(sub$Tissue,levels = c('nLN','pLN'))
ggboxplot(sub,x="Tissue",y="Ratio of CD4_C7_OX40 in T",color = "Tissue",palette = c("#FC8D62","#8DA0CB"),add = "jitter",add.params = list(size=1.5))+
  stat_compare_means(method='anova',paired = FALSE)

my<-subset(sc,L1_C=='Endothelium')
n<-subset(my,tissue=='nLN');p<-subset(my,tissue=='pLN')
t<-as.data.frame(table(n@meta.data$patient,n@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(n@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'nLN'
com<-t
t<-as.data.frame(table(p@meta.data$patient,p@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(p@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'pLN'
com<-rbind(com,t)
my_comparisons<-list(c('nLN','pLN')) 
sub<-com[com$Subcluster=='Endo_C3_RGCC',]
colnames(sub)[5]<-"Ratio of Endo_C3_RGCC in Endothelium";sub$Tissue<-factor(sub$Tissue,levels = c('nLN','pLN'))
ggboxplot(sub,x="Tissue",y="Ratio of Endo_C3_RGCC in Endothelium",color = "Tissue",palette = c("#FC8D62","#8DA0CB"),add = "jitter",add.params = list(size=1.5))+
  stat_compare_means(method='anova',paired = FALSE)

my = subset(sc, L3_C == 'Macrophage' | L3_C == 'Monocyte')
n<-subset(my,tissue=='nLN');p<-subset(my,tissue=='pLN')
t<-as.data.frame(table(n@meta.data$patient,n@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(n@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'nLN'
com<-t
t<-as.data.frame(table(p@meta.data$patient,p@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(p@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'pLN'
com<-rbind(com,t)
my_comparisons<-list(c('nLN','pLN')) 
sub<-com[com$Subcluster=='Mac_C2_SPP1',]
colnames(sub)[5]<-"Ratio of Mac_C2_SPP1 in Macrophage+Monocyte";sub$Tissue<-factor(sub$Tissue,levels = c('nLN','pLN'))
ggboxplot(sub,x="Tissue",y="Ratio of Mac_C2_SPP1 in Macrophage+Monocyte",color = "Tissue",palette = c("#FC8D62","#8DA0CB"),add = "jitter",add.params = list(size=1.5))+
  stat_compare_means(method='anova',paired = FALSE)

my<-subset(sc,L1_C=='Fibroblast')
n<-subset(my,tissue=='nLN');p<-subset(my,tissue=='pLN')
t<-as.data.frame(table(n@meta.data$patient,n@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(n@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'nLN'
com<-t
t<-as.data.frame(table(p@meta.data$patient,p@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(p@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'pLN'
com<-rbind(com,t)
my_comparisons<-list(c('nLN','pLN')) 
sub<-com[com$Subcluster=='FB_C3_COL1A1',]
colnames(sub)[5]<-"Ratio of FB_C3_COL1A1 in Fibroblast";sub$Tissue<-factor(sub$Tissue,levels = c('nLN','pLN'))
ggboxplot(sub,x="Tissue",y="Ratio of FB_C3_COL1A1 in Fibroblast",color = "Tissue",palette = c("#FC8D62","#8DA0CB"),add = "jitter",add.params = list(size=1.5))+
  stat_compare_means(method='anova',paired = FALSE)


# Supplementary Fig. 12 b -----
overlap = intersect(VDJ_lpf$clone.id,VDJ_tf$clone.id)
VDJ_lt$share = ifelse(VDJ_lt$clone.id %in% overlap,paste0("share",VDJ_lt$tissue),VDJ_lt$tissue )
tab.1 <- table(VDJ_lt$share,VDJ_lt$patient)
tab.1 <- tab.1[,c(3,4,6,7)]
tab.2 <- melt(tab.1)
colnames(tab.2) <- c("type","patient","number")
tab.2$tt <- tab.2$type
tab.2$tt <- gsub("share","",tab.2$tt)
tab.2 %>% group_by(type,patient)
tab.2<- ddply(tab.2,.(patient),transform,percent=number/sum(number)*100) 
tab.2<- ddply(tab.2,.(patient,tt),transform,percent2=number/sum(number)*100) 
tab.2$label = paste0(sprintf("%.1f", tab.2$percent), "%")
tab.2$label2 = paste0(sprintf("%.1f", tab.2$percent2), "%")
tab.2$sub_type = rep(pro,nrow(tab.2))
tab.2$percentage = tab.2$percent/100
tab.2$percentage2 = tab.2$percent2/100 
tab.2 %>%
    ggplot(aes(fill=type,x=patient,y=percentage,label = label))+
    geom_bar(position = position_stack(),stat = "identity")+
    geom_text(aes(label =label), 
              position = position_stack(vjust = 0.5), size = 6)+
    scale_y_continuous(labels = scales::percent)+
    # geom_text(aes(label=label),vjust=3,size=6,color="black")+
    theme_classic()+
    theme(legend.position = "top",
          legend.text = element_text(size=12),
          axis.text = element_text(size=12),
          axis.title = element_text(size=12,hjust = 0.5),
          plot.title  = element_text(hjust = 0.5,size = 16,face = "bold"),
          axis.text.x = element_text(size = 20,face = "bold"),
          axis.text.y = element_text(size = 20,face = "bold"),
          axis.title.x = element_text(size = 20,face = "bold"),
          axis.title.y = element_text(size = 20,face = "bold"),
          legend.title = element_text(size = 20,face = "bold")
    )+labs(title = pro)+scale_fill_d3()
  
  
# Supplementary Fig. 12 c -----
#Endothelium
my<-subset(sc,L1_C=='Endothelium' & metastasis == 'Y')
u<-subset(my,tissue=='Tu');p<-subset(my,tissue=='pLN')
t<-as.data.frame(table(u@meta.data$patient,u@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(u@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'Tu'
com<-t
t<-as.data.frame(table(p@meta.data$patient,p@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(p@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'pLN'
com<-rbind(com,t)
my_comparisons<-list(c('pLN','Tu')) 
sub<-com[com$Subcluster=='Endo_C3_RGCC',]
colnames(sub)[5]<-"Ratio of Endo_C3_RGCC in Endothelium";sub$Tissue<-factor(sub$Tissue,levels = c('pLN','Tu'))
ggboxplot(sub,x="Tissue",y="Ratio of Endo_C3_RGCC in Endothelium",color = "Tissue",palette = c("#8DA0CB","#c87137"),add = "jitter",add.params = list(size=1.5))+
  stat_compare_means(method='anova',paired = FALSE)

#Fibroblast
my = subset(sc, L1_C=='Fibroblast');my<-subset(my,metastasis == 'Y')
u<-subset(my,tissue=='Tu');p<-subset(my,tissue=='pLN')
t<-as.data.frame(table(u@meta.data$patient,u@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(u@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'Tu'
com<-t
t<-as.data.frame(table(p@meta.data$patient,p@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(p@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'pLN'
com<-rbind(com,t)
my_comparisons<-list(c('pLN','Tu')) 
sub<-com[com$Subcluster=='FB_C3_COL1A1',]
colnames(sub)[5]<-"Ratio of FB_C3_COL1A1 in Fibroblast";sub$Tissue<-factor(sub$Tissue,levels = c('pLN','Tu'))
ggboxplot(sub,x="Tissue",y="Ratio of FB_C3_COL1A1 in Fibroblast",color = "Tissue",palette = c("#8DA0CB","#c87137"),add = "jitter",add.params = list(size=1.5))+
  stat_compare_means(method='anova',paired = FALSE)

# Supplementary Fig. 12 d -----
my<-subset(sc,L1_C=='Endothelium')
n<-subset(my,tissue=='nLN');p<-subset(my,tissue=='pLN')
t<-as.data.frame(table(n@meta.data$patient,n@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(n@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'nLN'
com<-t
t<-as.data.frame(table(p@meta.data$patient,p@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(p@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'pLN'
com<-rbind(com,t)
my_comparisons<-list(c('nLN','pLN')) 
sub<-com[com$Subcluster=='Endo_C3_RGCC',]
colnames(sub)[5]<-"Ratio of Endo_C3_RGCC in Endothelium";sub$Tissue<-factor(sub$Tissue,levels = c('nLN','pLN'))
ggboxplot(sub,x="Tissue",y="Ratio of Endo_C3_RGCC in Endothelium",color = "Tissue",palette = c("#FC8D62","#8DA0CB"),add = "jitter",add.params = list(size=1.5))+
  stat_compare_means(method='anova',paired = FALSE)

my<-subset(sc,L1_C=='Fibroblast')
n<-subset(my,tissue=='nLN');p<-subset(my,tissue=='pLN')
t<-as.data.frame(table(n@meta.data$patient,n@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(n@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'nLN'
com<-t
t<-as.data.frame(table(p@meta.data$patient,p@meta.data$L4_C))
colnames(t)<-c('Patient','Subcluster','Cell_num')
total<-as.data.frame(table(p@meta.data$patient))
colnames(total)<-c('Patient','total_num');t<-merge(t,total,by='Patient')
t$ratio<-t$Cell_num/t$total_num*100
t<-t[!is.nan(t$ratio),]
t$Tissue<-'pLN'
com<-rbind(com,t)
my_comparisons<-list(c('nLN','pLN')) 
sub<-com[com$Subcluster=='FB_C3_COL1A1',]
colnames(sub)[5]<-"Ratio of FB_C3_COL1A1 in Fibroblast";sub$Tissue<-factor(sub$Tissue,levels = c('nLN','pLN'))
ggboxplot(sub,x="Tissue",y="Ratio of FB_C3_COL1A1 in Fibroblast",color = "Tissue",palette = c("#FC8D62","#8DA0CB"),add = "jitter",add.params = list(size=1.5))+
  stat_compare_means(method='anova',paired = FALSE)

# Supplementary Fig. 12 e -----
my = subset(sc, L3_C == 'Cancer' & metastasis == 'Y');Idents(my)<-my$tissue
deg<-FindMarkers(my,assay='RNA',logfc.threshold = 0,ident.1 = 'pLN',ident.2 = 'Tu',only.pos = F,max.cells.per.ident = 1000)
logFC_cutoff <- 0.3
deg$change = as.factor(
  ifelse(deg$p_val < 0.05 & abs(deg$avg_log2FC) > logFC_cutoff,
         ifelse(deg$avg_log2FC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
dat<-deg;dat$gene_name<-rownames(dat)
m2d = function(x){
  mean(abs(x))+2*sd(abs(x))
}  
gene<-c('CD24','RAC1','EGFR','S100A9','CXCL17','S100A8','SERPINE1','MDK','STMN1','TOP2A','TUBA1B','TUBB','TUBA1A','HLA-DRB1','HLA-DRA','HLA-A','CD74')
dat<-arrange(dat,desc(abs(avg_log2FC)))
dat<-dat[!(duplicated(dat$gene_name)),]
dat$gene <- as.factor(ifelse(dat$gene_name %in% gene, 'Y', 'N'))
dat$labels <- ''; for (i in gene) {dat[dat$gene_name == i, "labels"] <- i}
ggplot(data = dat, aes(x = avg_log2FC, y = -log10(p_val))) + 
  geom_point(alpha = 0.4, size = 3.5,aes(color = change)) + 
  geom_text_repel(aes(label=labels),size=5,max.overlaps = 10000,force=2,show.legend=F)+
  scale_color_manual(values = c("blue", "gray", "red")) + 
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), lty = 4, col = "black", lwd = 0.8) + 
  geom_hline(yintercept = -log10(0.05), lty = 4,col = "black", lwd = 0.8) + 
  theme_bw() + 
  labs(x = "log2FC", y = "-log10(P)") + 
  theme(plot.title = element_text(hjust = 0.5),panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),axis.text= element_text(colour = "black",size=14),axis.title = element_text(colour = "black",size=14),
        legend.text = element_text(colour = "black",size=14),legend.title = element_text(colour = "black",size=14))


# Supplementary Fig. 12 f -----
library(GseaVis)
library(clusterProfiler)
my = subset(sc, L3_C == 'Cancer' & metastasis == 'Y');Idents(my)<-my$tissue
deg<-FindMarkers(my,assay='RNA',ident.1 = 'pLN',ident.2='Tu',only.pos = T);deg$Gene<-rownames(deg)
rank<-deg[,c('Gene','avg_log2FC')]
geneList<-rank$avg_log2FC
names(geneList)=rank$Gene 
geneList=sort(geneList,decreasing = TRUE) 
human_H<-msigdbr(species = "human", category = "H")
human_C2<-msigdbr(species = "human", category = "C2")
human_C2<-human_C2[human_C2$gs_subcat%in%c('CP:BIOCARTA','CP:KEGG','CP:REACTOME'),]
human<-rbind(human_C2,human_H);gmt<-human %>% dplyr::select(gs_name,gene_symbol)
gseaRes <- GSEA(geneList = geneList,TERM2GENE = gmt, minGSSize = 0, maxGSSize = 1000,
                pvalueCutoff = 1, pAdjustMethod = "BH",verbose = FALSE)
geneSetID = c('REACTOME_DEVELOPMENTAL_BIOLOGY','REACTOME_ADAPTIVE_IMMUNE_SYSTEM')
gseaNb(object = gseaRes,geneSetID = geneSetID,addPoint = F, addPval = T,
       pCol = 'black', pHjust = 0,subPlot = 2,pvalX = 1,pvalY = 1.1,curveCol = brewer.pal(n = 2, name = "Set1"))

# S12g
setwd('/public/home/xlwang454/proj_paper_revise/ng/')
library(readxl);library(dplyr);library(tidyr);library(ggplot2);library(patchwork);library(openxlsx)

df <- readxl::read_xlsx('GEMCD8_TcelldistributionafterQC.xlsx')
endpoints <- c("%CD8A"="% CD8A Positive Cells", "%CD8+PD1+"="% CD8+ PD1+ Positive Cells",
               "%CD8+GEM+"="% CD8+GEM+ Positive Cells", "%CD8+PD1+GEM+"="% CD8+PD1+GEM+ Positive Cells",
               "CD8PD1GEM/CD8 (prop)"="prop_CD8PD1GEM_of_CD8", "CD8PD1GEM/CD8PD1 (prop)"="prop_CD8PD1GEM_of_CD8PD1")
df_pt <- df %>% dplyr::group_by(ID,tissue) %>% dplyr::summarise(across(all_of(unname(endpoints)), ~mean(as.numeric(.x),na.rm=TRUE)),.groups="drop")
paired_ids <- df_pt %>% dplyr::count(ID, name="n_tissue") %>% dplyr::filter(n_tissue==2) %>% dplyr::pull(ID)
df_pt <- df_pt %>% dplyr::mutate(is_paired = ID %in% paired_ids)

df_pt <- df_pt %>% dplyr::mutate(tissue = factor(tissue, levels=c("Tumor","Adjacent")))

stats_tbl <- lapply(names(endpoints), function(lbl) { col <- endpoints[[lbl]]
  x <- df_pt %>% dplyr::filter(tissue == "Tumor") %>% dplyr::pull(!!sym(col)) %>% as.numeric()
  y <- df_pt %>% dplyr::filter(tissue == "Adjacent") %>% dplyr::pull(!!sym(col)) %>% as.numeric()
  p_unpaired <- suppressWarnings(wilcox.test(x, y, paired = FALSE, exact = FALSE)$p.value)
  wide <- df_pt %>% dplyr::filter(ID %in% paired_ids) %>% dplyr::select(ID,tissue, !!sym(col)) %>%
    tidyr::pivot_wider(names_from = tissue, values_from = !!sym(col)) %>% tidyr::drop_na()
  p_paired <- if (nrow(wide) >= 3) {
    suppressWarnings(wilcox.test(wide$Tumor, wide$Adjacent, paired=TRUE, exact=FALSE)$p.value)} else {NA_real_}

  tibble(endpoint = lbl, p_unpaired_MW = p_unpaired, p_paired_Wilcoxon = p_paired, n_pairs_used = nrow(wide)) }) %>%
  dplyr::bind_rows() %>% dplyr::mutate(q_unpaired_MW=p.adjust(p_unpaired_MW,method="BH"), q_paired_Wilcoxon=p.adjust(p_paired_Wilcoxon, method="BH"))

transform_for_plot <- function(lbl, x) {
  if (startsWith(lbl, "%")) {  list(y = log10(x + 0.01), ylab = "log10(% + 0.01)")} else {list(y = x, ylab = "proportion")  }}
make_panel <- function(lbl, colname) {
  dd <- df_pt %>% dplyr::select(ID, tissue, is_paired, value = !!sym(colname)) %>% dplyr::mutate(value = as.numeric(value))
  tf <- transform_for_plot(lbl, dd$value); dd <- dd %>% dplyr::mutate(value_plot = tf$y)
  paired_long <- dd %>% dplyr::filter(ID %in% paired_ids) %>% dplyr::select(ID, tissue, value_plot) %>% tidyr::drop_na()
  row <- stats_tbl %>% dplyr::filter(endpoint == lbl)
  ann <- sprintf("paired q=%s\nunpaired q=%s\npaired n=%d", signif(row$q_paired_Wilcoxon, 3), signif(row$q_unpaired_MW, 3), row$n_pairs_used)
  
  set.seed(1)
  ggplot(dd, aes(x = tissue, y = value_plot)) + geom_violin(trim = FALSE, alpha = 0.25, color = "black", linewidth = 0.4) +
    geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.4, color = "black", linewidth = 0.4) +
    geom_point(data = dd, position = position_jitter(width = 0.08, height = 0), size = 1.6, alpha = 0.45) +
    geom_line(data = paired_long, aes(group = ID), linewidth = 0.35, alpha = 0.9) + geom_point(data = paired_long,size = 2.0, alpha=0.9) +
    scale_x_discrete(labels=c("Tumor" = "Tumor", "Adjacent" = "Adjacent")) + labs(title=lbl, x=NULL, y=tf$ylab) +
    annotate("text", x = -Inf, y = Inf, label = ann, hjust = -0.05, vjust = 1.15, size = 3) + theme_classic(base_size = 11) +
    theme(plot.title = element_text(size = 12, face = "plain"), axis.text.x = element_text(size = 10),axis.title.y = element_text(size = 10)) }

p_list <- mapply(make_panel, names(endpoints), unname(endpoints), SIMPLIFY = FALSE)
fig_main <- (p_list[[1]] | p_list[[2]] | p_list[[3]]) / (p_list[[4]] | p_list[[5]] | p_list[[6]]) +
  plot_annotation(title = "Tumor vs Adjacent (patient-level mean; violin+box with paired connections)")
ggsave("Figure1_Tumor_vs_Adjacent_violin_paired.pdf", fig_main, width = 12.5, height = 7.5, units = "in", device = cairo_pdf)

for (i in seq_along(p_list)) {
  fname <- paste0("Figure1_", gsub("[^A-Za-z0-9_+]+", "_", names(endpoints)[i]), ".pdf")
  ggsave(fname, p_list[[i]], width = 4.2, height = 3.6, units = "in", device = cairo_pdf) }
#