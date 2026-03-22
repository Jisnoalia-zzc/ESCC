library(pheatmap)

total = read.csv("nor.csv",header = T, row.names = 1)

annotation = read.csv("annotation.csv",header = T, row.names = 1)

bk <- c(seq(-3,-0.01,by=0.01),seq(0,3,by=0.01))

pheatmap(total, cluster_cols = FALSE, cluster_rows = FALSE,annotation_row = annotation_row,color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),breaks=bk,show_rownames = TRUE,show_colnames = FALSE)
