#Supplementary Fig. 13 a -----
co=subset(sc,subset =L4_C=='CD8_C6_CD39')
FeatureScatter(co, pt.size = 2, feature1 = 'GEM', feature2 = 'ENTPD1') +stat_cor(method = "pearson",size=6,vjust=0.5)+
  stat_regline_equation(aes(label=paste(..eq.label..)),size=6,vjust=1.6)+stat_smooth(method='lm', color="black", se=TRUE)+labs(title='')+guides(size=FALSE)+theme(legend.title = element_blank(),legend.text = element_text(size=15))+
  scale_color_manual(values ='#F37F4F')


#Supplementary Fig. 13 b -----
load('ESCC_RNA_FPKM_edit.RData')
head(meta)
library(ggrepel)
library(ggrastr)

exprSet %>% rownames_to_column("gene") %>% 
  filter(gene %in%c("GEM",'PDCD1','CD8A'))  %>%
  column_to_rownames("gene") %>% t()  %>%
  as.data.frame()->tmp

tmp$gem  = tmp$GEM/tmp$CD8A
tmp$pdcd1 = tmp$PDCD1/tmp$CD8A

ggplot(tmp, aes(x=log2(gem), y=log2(pdcd1))) +
  geom_point(color="#F37F4F",size=4)+ 
  stat_cor(method = "pearson",size=6,vjust=0.5,digits = 2)+
  stat_smooth(method='lm', 
              color="black", se=TRUE)+labs(title='')+guides(size=FALSE)+
  theme(legend.title = element_blank(),legend.text = element_text(size=15))+
  theme_classic()+
  theme(axis.text = element_text(size=12,colour = "black"),
        axis.title = element_text(size=15))+
  labs(x="log2(GEM/CD8A FPKM)", y="log2(PDCD1/CD8A FPKM")


#Supplementary Fig. 13 c -----
DefaultAssay(sp)<-'spatial'
p<-subset(sp,CD8A>0&ENTPD1>0&GEM>0&PDCD1>0)
DefaultAssay(p)<-'spatial';p$group<-'CD8+CD39+';Idents(p)<-p$group
FeatureScatter(p, pt.size = 2, feature1 = 'GEM', feature2 = 'PDCD1') +stat_cor(method = "pearson",size=6,vjust=0.5)+
  stat_regline_equation(aes(label=paste(..eq.label..)),size=6,vjust=1.6)+stat_smooth(method='lm', color="black", se=TRUE)+labs(title='')+guides(size=FALSE)+theme(legend.title = element_blank(),legend.text = element_text(size=15))+
  scale_color_manual(values ='#F37F4F')


# Supplementary Fig. 13 g
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