###S9A
library(data.table)
library(tibble)
library(dplyr)
library(ggplot2)

# 1. 
setwd('C:/Users/Ron/Desktop/Figs/S9')
total_df <- fread('C:/Users/Ron/Desktop/Figs/S9/S9_a_TLS_Group_and_mimer_mark_on_slice_data.csv')

# 2. 
all_slices <- unique(total_df$slice_name)

# 3. 
library(patchwork)
plots_list <- list()

for (slice in all_slices) {
  slice_data <- total_df %>% filter(slice_name == slice)
  
  if (all(slice_data$has_tls)) {
    p <- ggplot(slice_data, aes(x = imagecol, y = imagerow)) +
      geom_point(aes(color = plot_group, alpha = plot_alpha), size = 0.8) +
      geom_point(
        data = subset(slice_data, mimer == "MIMER"),
        aes(shape = "MIMER"),
        color = "black",
        size = 0.5
      ) +
      scale_color_manual(
        name = "TLS Group",
        values = c( "G1" = "#E41A1C", "G2" = "#377EB8","G3" = "#4DAF4A", "Not TLS" = "gray80")
        
      ) +
      scale_alpha_identity() +
      scale_shape_manual(values = c("MIMER" = 3), name = "") +
      theme_minimal() +
      labs(title = paste(slice, "(with TLS)"), x = "X", y = "Y") +
      coord_fixed()
  } else {
    p <- ggplot(slice_data, aes(x = imagecol, y = imagerow)) +
      geom_point(aes(color = plot_group), size = 0.8, alpha = 0.4) +
      geom_point(
        data = subset(slice_data, mimer == "MIMER"),
        aes(shape = "MIMER"),
        color = "black",
        size = 0.5
      ) +
      scale_color_manual(values = c("Not TLS" = "gray80")) +
      scale_shape_manual(values = c("MIMER" = 3), name = "") +
      theme_minimal() +
      labs(title = paste(slice, "(without TLS)"), x = "X", y = "Y") +
      coord_fixed()
  }
  
  plots_list[[slice]] <- p
}

# 4.
grid_plot <- wrap_plots(plots_list, ncol = 4, nrow = ceiling(length(plots_list)/4))
print(grid_plot)


######ME #######
load("/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhangzhcihao/spatial_pathology/sub.dist_ME.Rdata")
sub.dist.all = do.call(rbind,sub.dist.ME)

sub.dist.all$cellID = sub.dist.all$id

load("/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhangzhcihao/spatial_pathology/sp_meta_info.Rdata")
sub.dist.all %>% left_join(meta,by = "cellID") -> plot_df
# table(plot_df$distinct_area)
colnames(plot_df)
# HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
plot_df$ME = factor(plot_df$ME,
                    levels =  c('Nor-ME','Hyp-ME','MiD-ME','MoD-ME','SD&CA-ME','ICA-ME','MCA-ME'))
hallmark_paths = c("pro_fibrotic_signature",'ECM',
                   'Anti_inflammatory','Immunosuppression',
                   'pro_metastasis',"pro_fibrotic_signature",
                   'ECM','Anti_inflammatory')
hallmark_list = lapply(hallmark_paths, function(path){
  ggplot(plot_df,aes(x=min,y=plot_df[[path]] ,color=plot_df[[path]] ) )+
    geom_point_rast(size=0.1) +
    scale_color_distiller(palette = "Spectral")+
    stat_cor(size = 2)+
    geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2",span=1, size = 0.5)+
    facet_grid(~ME)+
    labs(y=path,x="Distance to Mimer (Spot)")+theme_classic() +
    theme(
      legend.key.size = unit(2, "mm"),
      panel.spacing = unit(1, "mm") ,
      axis.line = element_line(size = 0.2),
      axis.ticks = element_line(size = 0.2),
      text = element_text(size = 5), 
      axis.text = element_text(size = 6),
      legend.text = element_text(size = 5),
      #plot.title = element_text(size = 10),
      strip.text.x = element_text(size = 6),
      strip.background = element_rect(size = 0.2)
      
    )->p1
  p1[["labels"]][["colour"]] = path
  return(p1)
})
names(hallmark_list) = hallmark_paths

library(patchwork)
wrap_plots(hallmark_list,ncol = 1)->p
ggsave(
  filename = "/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhangzhcihao/spatial_pathology/distinct/plot_A4_3.pdf",  # 支持 PDF/PNG/TIFF 等格式
  plot = p,device = "pdf",width = 210,height = 297,units = "mm", dpi = 300 )

##### S9B
# S9b
library(dplyr)
library(ggplot2)
library(tidyr)

# 1. 
setwd('C:/Users/Ron/Desktop/Figs/S9')
TLS_mimer_sorted <- read.csv('C:/Users/Ron/Desktop/Figs/S9/S9_b_TLS_mimer_sorted_wide.csv')

# 2. 
TLS_mimer_long <- TLS_mimer_sorted %>%
  select(patient, TLS_ID, TLS_group, MIMER_prop, Others_prop) %>%
  pivot_longer(
    cols = c(MIMER_prop, Others_prop),
    names_to = "cell_type",
    values_to = "proportion"
  ) %>%
  mutate(
    cell_type = factor(
      cell_type, 
      levels = c("Others_prop", "MIMER_prop"),
      labels = c("Others", "MIMER")
    ),
    TLS_ID = factor(TLS_ID, levels = unique(TLS_ID))
  )

# 3. 
group_counts <- table(TLS_mimer_sorted$TLS_group)

cum_counts <- cumsum(group_counts)
xintercepts <- cum_counts[-length(cum_counts)] + 0.5

label_x <- as.numeric(cum_counts - group_counts/2)

group_labels <- names(group_counts)

# 4. 
options(repr.plot.width=15, repr.plot.height=10)
#pdf(file = 'S9_b_Mimer_propotions_by_TLS_group_from_saved.pdf', 8, 5)

p <- ggplot(TLS_mimer_long, aes(x = TLS_ID, y = proportion, fill = cell_type)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c("MIMER" = "#bd2628", "Others" = "#6495ed"),
                    name = "Cell Type") +
  labs(title = "MIMER Proportions by TLS Group", x = "TLS ID", y = "Proportion") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) + geom_vline(
    xintercept = xintercepts,
    linetype = "dashed", color = "gray50", alpha = 0.7
  ) + annotate(
    "text",
    x = label_x, 
    y = -0.05,
    label = group_labels,
    size = 3,
    color = "black"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.05)))

print(p)
#dev.off()


#### S9D
# S9d
library(data.table)
library(tidyverse)
library(ggplot2)
library(patchwork)

setwd('C:/Users/Ron/Desktop/Figs/S9/')

annotation_df_mark_mimer_pos_G3_add_TLS_flare <- fread('C:/Users/Ron/Desktop/Figs/S9/S9_d_annotation_df_mark_mimer_pos_G3_add_TLS_flare.csv')
all_cols <- names(annotation_df_mark_mimer_pos_G3_add_TLS_flare)
start_idx <- which(all_cols == "cluster_size")
end_idx <- which(all_cols == "peripheral_Tumorsuppression_avg")
plot_vars <- c(all_cols[start_idx:end_idx],'core_TLS_G3_avg','peri_TLS_G3_avg')

plot_data <- annotation_df_mark_mimer_pos_G3_add_TLS_flare %>%
  select(TLS_ID, is_mimer_pos_tls, all_of(plot_vars)) %>%
  filter(is_mimer_pos_tls %in% c("other_G3", "mimer_G3"))

create_single_plot <- function(var_name, data) {
  var_data <- data %>%
    select(TLS_ID, is_mimer_pos_tls, value = all_of(var_name)) %>%
    filter(!is.na(value))
  p_value <- NA
  sig_label <- "NA"
  
  if (length(unique(var_data$is_mimer_pos_tls)) == 2) {
    test_result <- tryCatch(
      wilcox.test(value ~ is_mimer_pos_tls, data = var_data),
      error = function(e) NULL
    )
    
    if (!is.null(test_result)) {
      p_value <- test_result$p.value
      sig_label <- case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    }
  }
  
  y_range <- range(var_data$value, na.rm = TRUE)
  y_max <- y_range[2]
  y_min <- y_range[1]
  y_span <- y_max - y_min
  
  y_limit_upper <- y_max + 0.2 * y_span
  
  p <- ggplot(var_data, aes(x = is_mimer_pos_tls, y = value, fill = is_mimer_pos_tls)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.2, show.legend = FALSE) +
    geom_text(
      x = 1.5,
      y = y_max + 0.1 * y_span,
      label = sig_label,
      size = 5,
      fontface = "bold",
      inherit.aes = FALSE
    ) +
    labs(
      title = var_name,
      y = "Value",
      x = ""
    ) +
    scale_fill_brewer(palette = "Set1") +
    coord_cartesian(ylim = c(y_min, y_limit_upper)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
    )
  
  return(p)
}

plot_list <- lapply(plot_vars, function(var) {
  create_single_plot(var, plot_data)
})

n_cols <- 6
combined_plot <- wrap_plots(plot_list, ncol = n_cols) +
  plot_annotation(
    title = "Comparison of TLS Features: other_G3 vs mimer_G3",
    caption = "*** p<0.001, ** p<0.01, * p<0.05, ns not significant"
  )

options(repr.plot.width = 18, repr.plot.height = 3 * ceiling(length(plot_vars)/n_cols))
print(combined_plot)

##### S9e
#CN9
md = read.table('score.txt',sep='\t',header=T,check.names = F)
res = md[md$Signature %in% 'Activated B',]
ggboxplot(res, x="type", y="val", color="type", palette="npg", add="jitter", bxp.errorbar=F, ylab='TLS scores', add.params=list(size=5)) +
  theme_classic() + stat_compare_means(comparisons=list(c('With','Without')), method='t.test', paired=F, size=5) +
  theme(axis.text=element_text(size=15,colour="black"),axis.title=element_text(size=18),axis.title.x=element_blank(),legend.position='none')
ggsave("Activated.B.CN9.pdf", width=4, height=6, device=cairo_pdf, dpi=300)



