library(ggplot2)
library(cowplot)
pdf(file="/lustre/home/xhzh/scSpatial/xiaodu/codes/zxh/MIA/HALLMARK_ME_DEV1.pdf", width=18, height=8)
meta <- read.table("HALLMARK_meta.data.txt", header=T, sep="\t")

hallmark <- as.character(colnames(meta[84:133]))


m1 <- margin(t=0.5, r=0.1, b=0.5, l=0.1, unit="cm")
theme3 <- theme(plot.title = element_text(size = 16, face = "bold"), axis.title.y = element_text(size=15), axis.title.x = element_text(size=15), axis.text.x = element_text(size = 12), axis.text.y = element_text(size =14), legend.title = element_blank(), legend.text = element_text(size = 12), legend.position="right", plot.margin=m1, panel.background = element_rect(fill='white', color = 'black'))

plots <- vector("list", 50)
plots1 <- vector("list", 50)

c <- 0
for (h in hallmark) {
    h1 <- gsub("HALLMARK_", "", h)
    c = c + 1
    meta1 <- meta[meta$development != "unknown",]
    meta1$development <- factor(meta1$development, levels = c('Nor', 'Hyp', 'MiD', 'MoD', 'SD&CA', 'ICA', 'MCA'))
    de <- aggregate(meta1[,h], by=list(type=meta1$development),FUN=function(x) c(mean=mean(x), sd=sd(x), n <- length(x), sq=sqrt(n))); de$dev <- de$type; de$group = "DEV"; de$mean <- de$x[,1]; de$sd <- de$x[,2]; de$sem <- de$x[,2]/de$x[,4]

    meta2 <- meta[meta$ME != "unknown",]
    meta2$ME <-factor(meta2$ME,levels = c('Nor-ME','Hyp-ME','MiD-ME','MoD-ME','SD&CA-ME','ICA-ME','MCA-ME'))
    me <- aggregate(meta2[,h], by=list(type=meta2$ME),FUN=function(x) c(mean=mean(x), sd=sd(x), n <- length(x), sq=sqrt(n))); me$dev <- gsub("-ME", "", me$type); me$group = "ME"; me$mean <- me$x[,1]; me$sd <- me$x[,2]; me$sem <- me$x[,2]/me$x[,4]

    merge <- rbind(de, me)
    merge$dev <- factor(merge$dev, levels = c('Nor', 'Hyp', 'MiD', 'MoD', 'SD&CA', 'ICA', 'MCA'))

    plots[[c]] <- ggplot(data=merge, aes(x=dev, y=mean, group=group)) + geom_point(aes(color=group)) + geom_line(aes(color=group)) + geom_errorbar(aes(ymin = mean-sem, ymax = mean+sem), width = 0.1, alpha=0.7) + labs(title=h1, x="Dev stages", y="Score (mean+-sem)") + theme3

    de1 <-  aggregate(meta1[,h], by=list(type=meta1$development,file=meta1$file),mean); de1$lab <- paste(de1$type, de1$file, sep="-")

    me1 <- aggregate(meta2[,h], by=list(type=meta2$ME, file=meta2$file),mean); me1$dev <- gsub("-ME", "", me1$type); me1$lab <- paste(me1$dev, me1$file, sep="-")

    de2 <- de1[de1$lab %in% me1$lab,]
    merge1 <- cbind(de2, me1)
    merge1$ratio <- merge1[,3] / merge1[,7]
    merge1$dev <- factor(merge1$dev, levels = c('Nor', 'Hyp', 'MiD', 'MoD', 'SD&CA', 'ICA', 'MCA'))
    de2me <- aggregate(merge1$ratio, by=list(dev1=merge1$dev), FUN=function(x) c(mean=mean(x), sd=sd(x), n <- length(x), sq=sqrt(n)));
    de2me$mean <- de2me$x[,1]; de2me$sd <- de2me$x[,2]; de2me$sem <- de2me$x[,2]/de2me$x[,4]
    plots1[[c]] <- ggplot(data=de2me, aes(x=dev1, y=mean, group=1)) + geom_point() + geom_line() + geom_errorbar(aes(ymin = mean-sem, ymax = mean+sem), width = 0.1, alpha = 0.7) + labs(title=h1, x="Dev stages", y="Ratio+-sem") + theme3
}

for ( i in c(1,7,13,19,25,31,37,43) ) {
    j = i + 5
    p <- plot_grid(plotlist = plots[c(i:j)], ncol = 3, nrow = 2)
    print(p)
}
plot_grid(plotlist = plots[c(49,50)], ncol = 3, nrow = 2)

for ( i in c(1,7,13,19,25,31,37,43) ) {
    j = i + 5
    p <- plot_grid(plotlist = plots1[c(i:j)], ncol = 3, nrow = 2)
    print(p)
}
plot_grid(plotlist = plots1[c(49,50)], ncol = 3, nrow = 2)
dev.off()
