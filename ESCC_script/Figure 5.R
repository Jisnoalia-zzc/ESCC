######################
#--------- Figure 5b----------
#######################

library(ggplot2)
library(dplyr)
library(ggpubr)
library(data.table)

file_stats <- fread('C:/Users/Ron/Desktop/Figs/F5/file_stats.csv')

if (exists("file_stats")) {
  diff_data <- file_stats %>%
    filter(tissue %in% c("LnF", "LpF")) %>%
    select(file, tissue, diff) %>%
    mutate(tissue = factor(tissue, levels = c("LnF", "LpF")))
  
  LnF_diff <- diff_data %>% filter(tissue == "LnF") %>% pull(diff)
  LpF_diff <- diff_data %>% filter(tissue == "LpF") %>% pull(diff)
  
  if (length(LnF_diff) > 0 && length(LpF_diff) > 0) {
    wilcox_test <- wilcox.test(diff ~ tissue, data = diff_data)
    
    stats_summary <- diff_data %>% group_by(tissue) %>%  summarise(
        n = n(),
        mean = mean(diff, na.rm = TRUE),
        median = median(diff, na.rm = TRUE),
        sd = sd(diff, na.rm = TRUE),
        se = sd / sqrt(n),
        .groups = "drop"
      )
    
    tissue_colors <- c("LnF" = "#87CEEB", "LpF" = "#CD5C5C")
    
    p_box <- ggplot(diff_data, aes(x = tissue, y = diff, fill = tissue)) +
   
      geom_boxplot(
        width = 0.6,
        alpha = 0.7,
        outlier.shape = NA,
        linewidth = 0.8
      ) +
      geom_jitter(
        aes(color = tissue),
        width = 0.2,
        height = 0,
        size = 2.5,
        alpha = 0.7,
        shape = 19
      ) +
      geom_point(
        data = stats_summary,
        aes(x = tissue, y = mean),
        shape = 23, 
        size = 4,
        fill = "black",
        color = "black",
        stroke = 1.2
      ) +
      geom_text(
        data = stats_summary,
        aes(x = tissue, y = mean, label = sprintf("Mean: %.3f", mean)),
        vjust = 1.8,
        size = 3.2,
        fontface = "bold",
        color = "darkblue"
      ) +
      scale_fill_manual(values = tissue_colors, name = "Tissue Type") +
      scale_color_manual(values = tissue_colors, guide = "none") +
      labs(
        title = "Comparison of Difference Values: LnF vs LpF",
        subtitle = sprintf("Wilcoxon rank-sum test: p = %.4f", wilcox_test$p.value),
        x = "Tissue Type",
        y = "Difference (MIMER mean - Random mean)",
        caption = sprintf("Sample size: LnF = %d, LpF = %d", 
                         stats_summary$n[stats_summary$tissue == "LnF"],
                         stats_summary$n[stats_summary$tissue == "LpF"])
      ) +
      ylim(
        min(diff_data$diff, na.rm = TRUE) * 1.1,
        max(diff_data$diff, na.rm = TRUE) * 1.1
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(size = 13, face = "bold", margin = margin(r = 10)),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
        plot.caption = element_text(hjust = 0.5, size = 10, color = "gray50"),
        plot.margin = margin(15, 15, 15, 15)
      )
    
    p_box <- p_box + 
      stat_compare_means(
        method = "wilcox.test",
        comparisons = list(c("LnF", "LpF")),
        label = "p.signif",
        symnum.args = list(
          cutpoints = c(0, 0.001, 0.01, 0.05, 1),
          symbols = c("***", "**", "*", "NS")
        ),
        vjust = 0.2,
        size = 5
      )
    
    print(p_box)
    
    cat("\n=== Statistical Summary ===\n")
    print(stats_summary)
    cat("\n=== Wilcoxon Rank-Sum Test Results ===\n")
    print(wilcox_test)
    
  } else {
    cat("Error: At least one group has zero samples\n")
  }
} else {
  cat("Error: file_stats data frame not found. Please run the previous code first.\n")
}

#Figure 5c-d -----
library("Startrac")
library("tictoc")
library("ggpubr")
library("ggplot2")
library("ComplexHeatmap")
library("RColorBrewer")
library("circlize")
tcr_meata <- read.csv("tcr_meta_filter.csv",header = T)
cd8 <- tcr_use[tcr_use$stype == 'CD8',]
out1 <- Startrac.run(cd8, proj="ESCC",verbose=F)
level_cd8 <- c("CD8_C1_CCR7","CD8_C2_IL7R","CD8_C3_TCF7","CD8_C4_ANXA1" ,"CD8_C5_GZMK" , "CD8_C6_PDCD1"  ,
              "CD8_C7_CX3CR1","CD8_C8_KIR", "CD8_C9_KLRC2", "CD8_C10_CD16" ,"CD8_C11_ISG15","CD8_C12_CD161"
)
out2 = as.data.table(out1@cluster.sig.data)[aid!=out1@proj,][order(majorCluster),]
out3 =out2[out2$index == "expa",]
out3$majorCluster <- factor(out3$majorCluster,
                            levels = level_cd8 )
library(ggpubr)
my_comparisons = list(
  c("CD8_C4_ANXA1","CD8_C6_PDCD1")
)
ggboxplot(out3,
          x="majorCluster",y="value",palette = col_cd8$color,
          color = "majorCluster", add = "point", outlier.colour=NULL) +
  stat_compare_means(comparisons = my_comparisons,paired = T,method ="t.test")+
  facet_wrap(~index,ncol=1,scales = "free_y") +
  theme(axis.text.x=element_text(angle = 60,hjust = 1),
        legend.position = 'none',
  )
ggboxplot(as.data.table(obj@pIndex.sig.migr)[aid!=obj@proj,][order(majorCluster),],
          x="majorCluster",y="value",palette = col_cd8$color,
          color = "majorCluster", add = "point", outlier.colour=NULL) +
  stat_compare_means(comparisons = my_comparisons,
                     paired = T,method ="t.test")+
  facet_wrap(~index,ncol=1,scales = "free_y") +
  theme(axis.text.x=element_text(angle = 60,hjust = 1),
        legend.position = 'none')  

#CD4T
cd4 <- tcr_use[tcr_use$stype == 'CD4',]
out1 <- Startrac.run(cd4, proj="ESCC",verbose=F)
level_cd4 <- c("CD4_C1_CCR7" ,"CD4_C2_TCF7" ,"CD4_C3_ANXA1","CD4_C4_GPR183","CD4_C5_RTKN2" ,
               "CD4_C6_CD25","CD4_C7_OX40", "CD4_C8_ITGB1","CD4_C9_IFNG" , "CD4_C10_CXCR5" ,
               "CD4_C11_IL17A" ,"CD4_C12_NKG7" , "CD4_C13_ISG15"
)
out2 = as.data.table(out1@cluster.sig.data)[aid!=out1@proj,][order(majorCluster),]
out3 =out2[out2$index == "expa",]
out3$majorCluster <- factor(out3$majorCluster,
                            levels = level_cd8 )
library(ggpubr)
my_comparisons = list(
      c("CD4_C3_ANXA1","CD4_C7_OX40")
)
ggboxplot(out3,
          x="majorCluster",y="value",palette = col_cd8$color,
          color = "majorCluster", add = "point", outlier.colour=NULL) +
  stat_compare_means(comparisons = my_comparisons,paired = T,method ="t.test")+
  facet_wrap(~index,ncol=1,scales = "free_y") +
  theme(axis.text.x=element_text(angle = 60,hjust = 1),
        legend.position = 'none',
  )
ggboxplot(as.data.table(obj@pIndex.sig.migr)[aid!=obj@proj,][order(majorCluster),],
          x="majorCluster",y="value",palette = col_cd8$color,
          color = "majorCluster", add = "point", outlier.colour=NULL) +
  stat_compare_means(comparisons = my_comparisons,
                     paired = T,method ="t.test")+
  facet_wrap(~index,ncol=1,scales = "free_y") +
  theme(axis.text.x=element_text(angle = 60,hjust = 1),
        legend.position = 'none')  

#Figure 5e-f-----
my = subset(sc, L3_C == 'CD4');Idents(my)<-my$L4_C
gene<-c('S1PR5','S1PR1','CCR8','CXCR3','CXCR6')
my@active.ident <- factor(my@active.ident, levels=sort(names(table(my@active.ident))))
DotPlot(my,features = gene,dot.scale = 6,scale.min = 0)+scale_colour_gradientn(colours = rev(getPalette(10))) +theme_bw() + xlab('') + ylab('') + coord_flip()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1),axis.text = element_text(size = 14,color = 'black'),strip.text.x = element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor = element_blank())

my = subset(sc, L3_C == 'CD8');Idents(my)<-my$L4_C
gene<-c('S1PR5','S1PR1','CXCR6')
my@active.ident <- factor(my@active.ident, levels=sort(names(table(my@active.ident))))
DotPlot(my,features = gene,dot.scale = 6,scale.min = 0)+scale_colour_gradientn(colours = rev(getPalette(10))) +theme_bw() + xlab('') + ylab('') + coord_flip()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1),axis.text = element_text(size = 14,color = 'black'),strip.text.x = element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor = element_blank())

#Figure 5g-----
my = subset(sc, L4_C == 'Mac_C2_SPP1' & metastasis == 'Y');Idents(my)<-my$tissue;my<-subset(my,ident=c('pLN','Tu'));my$tissue<-Idents(my)
gene<-c('MARCO','CCL2','SPP1','CCL7','MRC1','IL1B')
my@active.ident <- factor(my@active.ident, levels=sort(names(table(my@active.ident))))
DotPlot(my,features = gene,dot.scale = 8)+scale_colour_gradientn(colours = rev(getPalette(10))) +theme_bw() + xlab('') + ylab('') + coord_flip()+
  theme(axis.text.x = element_text(angle = 0,hjust = 1),axis.text = element_text(size = 14,color = 'black'),strip.text.x = element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor = element_blank())


#Figure 5h -----
load("ellphonedb.RData")
library(stringr)
library(tidyverse)
library(circlize)
library(reshape2)
com$L = str_split(com$interacting_pair,"_",simplify = T)[,1]
com$R = str_split(com$interacting_pair,"_",simplify = T)[,2]

com$sender = str_split(com$celltype_pairs,"\\.",simplify = T)[,1]
com$recptor = str_split(com$celltype_pairs,"\\.",simplify = T)[,2]

df = com %>% select(c(L,R,sender,mean))

df$L = paste0(df$L,"_",df$sender)

df$mean = as.numeric(df$mean)

rownames(df) <- 1:nrow(df)

fd <- dcast(df,L~R,value.var = "mean")

mat2 <- matrix(rnorm(25), nrow = 5)
rownames(mat2) <- paste0("A", 1:5)
colnames(mat2) <- paste0("C", 1:5)

mat3 <- matrix(rnorm(25), nrow = 5)
rownames(mat3) <- paste0("B", 1:5)
colnames(mat3) <- paste0("C", 1:5)

mat <- matrix(0, nrow = 10, ncol = 10)
rownames(mat) <- c(rownames(mat2), rownames(mat3))
colnames(mat) <- c(colnames(mat1), colnames(mat2))
mat[rownames(mat1), colnames(mat1)] <- mat1
mat[rownames(mat2), colnames(mat2)] <- mat2
mat[rownames(mat3), colnames(mat3)] <- mat3

nm <- unique(unlist(dimnames(mat)))
group <- structure(gsub("\\d", "", nm), names = nm)
grid.col <- structure(
  c(rep("#fb8072", 5), rep("#80b1d3", 5), rep("#fdb462", 5)),
  names = c(paste0("A", 1:5), paste0("B", 1:5), paste0("C", 1:5))
)
group <- structure(gsub("\\d", "", nm), names = nm)
chordDiagram(
  mat, group = group, grid.col = grid.col, 
  annotationTrack = c("grid", "axis"),
  preAllocateTracks = list(
    track.height = mm_h(4),
    track.margin = c(mm_h(4), 0)
  )
)

circos.track(
  track.index = 2, 
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    circos.text(
      mean(xlim), mean(ylim),
      sector.index, cex = 0.6,
      niceFacing = TRUE
    )
  }, 
  bg.border = NA
)

highlight.sector(
  rownames(mat1), track.index = 1, col = "#fb8072", 
  text = "A", cex = 0.8, text.col = "white", 
  niceFacing = TRUE
)
highlight.sector(
  colnames(mat1), track.index = 1, col = "#80b1d3", 
  text = "B", cex = 0.8, text.col = "white", 
  niceFacing = TRUE
)
highlight.sector(
  colnames(mat2), track.index = 1, col = "#fdb462", 
  text = "C", cex = 0.8, text.col = "white", 
  niceFacing = TRUE
)
circos.clear()