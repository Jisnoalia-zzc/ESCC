#### E7A
load("/realspace/project/proj_ESCC_STW_ZWM_2022_01/zzc/pcc/samples_corr-20240708.RData")
res
#
n_CM = 7
pdf("cellratio_person_coorelation_ward_D2-20240708.pdf",width = 12,height = 6)
res1 = pheatmap(res$r,clustering_distance_rows = 'correlation',clustering_distance_cols = 'correlation',
                clustering_method = 'ward.D2',color = colorRampPalette(rev(RColorBrewer::brewer.pal(11,'RdBu')))(100),
                breaks=seq(-1,1,length.out=101),silent = F,border_color = NA,
                # display_numbers =res0,display_numbers = T,
                number_color = 'black',cutree_cols = n_CM,cutree_rows =n_CM )

res1
dev.off()
#
mp = as.data.frame(readxl::read_xlsx('SampleInfo.xlsx', col_names=T)); mp = tibble::column_to_rownames(mp, var='patient')
mc = c('CD8_C6_PDCD1','FB_C3_COL1A1','Endo_C3_RGCC','Mac_C3_C1QC','FB_C5_PDGFRB','Plasma','Neutrophil','CD4_C7_OX40','Mac_C2_SPP1','FB_C4_APOE',
       'CD4_C6_CD25','DC_C3_CD1A','Mac_C6_MKI67','B_C5_ISG15','CD8_C8_KIR','DC_C5_IL3RA','CD4_C13_ISG15','CD8_C11_ISG15','CD4_C10_CXCR5',
       'B_C6_BCL6','B_C7_MKI67','CD4_C9_IFNG','DC_C4_LAMP3','DC_C1_CLEC9A','Mac_C5_CXCL10','CD4_C11_IL17A','CD8_C10_CD16','CD8_C4_ANXA1',
       'Endo_C4_CCL21','B_C2_DUSP4','DC_C2_CD1C','CD4_C12_NKG7','CD4_C8_ITGB1','Mac_C4_LYVE1','Mono_C1_CD14','FB_C1_CFD','FB_C2_IGF1',
       'CD4_C2_TCF7','CD8_C1_CCR7','CD8_C3_TCF7','CD8_C9_KLRC2','CD4_C4_GPR183','CD4_C1_CCR7','CD4_C5_RTKN2','B_C4_TCL1A','CD8_C5_GZMK',
       'Endo_C2_FBLN5','FB_C6_ACTA2','Mac_C1_NLRP3','Mono_C2_CD16','B_C3_GPR183','CD8_C7_CX3CR1','CD4_C3_ANXA1','Endo_C1_ACKR1','CD8_C12_CD161',
       'B_C1_CCR7','CD8_C2_IL7R')
colnames(abun)[56] = 'CD8_C6_PDCD1'; myd = t(abun); myd = myd[match(mc,rownames(myd)),]
colnames(myd)=c("1019","1022","1104","1119","1123","0118","1201","1204","1209","1210","1229","0128","0308","0826")
pheatmap(myd, cluster_rows=F, clustering_distance_cols='correlation', clustering_method='ward.D2', border_color=NA, number_color='black',
         color=colorRampPalette(RColorBrewer::brewer.pal(9,'Reds'))(100), breaks=seq(0,1,length.out=101), silent=F, annotation_col=mp)


## E7B
library(dplyr); library(tidyr); library(ggplot2); library(ggpubr)

pair_data <- read.csv("E7_b_M1_M7_TF_only_pair_data.csv", stringsAsFactors = FALSE)

plot_module_violin <- function(module_name, data) {
  
 mod_data <- data %>% filter(Module == module_name)
  
  if (nrow(mod_data) < 2) {
    warning(paste(module_name, "样本数不足，跳过绘图"))
    return(NULL)
  }
  
  long_data <- mod_data %>%
    pivot_longer(cols = c(MIMER, Permutation), names_to = "Group",values_to = "Moran_I" ) %>%
    mutate(
      GroupLabel = ifelse(Group == "MIMER",module_name, "Permutation\n(mean of 1000 times)"),
      GroupLabel = factor(GroupLabel,levels = c(module_name, "Permutation\n(mean of 1000 times)"))
    )
  
  wilcox_test <- wilcox.test(mod_data$MIMER, mod_data$Permutation, paired = TRUE)
  p_value <- wilcox_test$p.value
  
  p <- ggplot(long_data, aes(x = GroupLabel, y = Moran_I, fill = GroupLabel)) + geom_violin(alpha = 0.7, width = 1.2) +
    geom_boxplot(width = 0.15, fill = "white", alpha = 0.8, outlier.shape = NA) +
    geom_line(aes(group = File), color = "gray80", alpha = 0.5, size = 0.8) + geom_point(aes(color = GroupLabel), size = 1, alpha = 0.6) +
    scale_fill_manual(values = setNames(c("#E41A1C", "#377EB8"), c(module_name, "Permutation\n(mean of 1000 times)") )) +
    scale_color_manual(values = setNames(c("#E41A1C", "#377EB8"), c(module_name, "Permutation\n(mean of 1000 times)"))) +
    theme_bw(base_size = 11) +
    labs(title = paste("Comparison of Moran's I:", module_name, "vs Permutation"),
      subtitle = paste0("Wilcoxon signed-rank test, p = ", format(p_value, digits = 3, scientific = TRUE)), x = "", y = "Moran's I") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold", margin = margin(b = 10)),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray30", margin = margin(b = 15)),
      axis.title = element_text(size = 18, face = "bold"),
      axis.text.x = element_text(size = 16, face = "bold", color = "black", margin = margin(t = 5)),
      axis.text.y = element_text(size = 16, color = "black"),axis.title.y = element_text(margin = margin(r = 10)),
      panel.grid.major.x = element_blank(),panel.grid.minor = element_blank(), legend.position = "none", plot.margin = margin(20, 20, 20, 20)) +
    stat_compare_means(paired = TRUE,method = "wilcox.test",label = "p.format",label.x = 1.5,label.y = max(long_data$Moran_I, na.rm = TRUE) - 
        (max(long_data$Moran_I, na.rm = TRUE) - min(long_data$Moran_I, na.rm = TRUE)) * 0.1,size = 6, vjust = -0.5 )
  
  return(p)
}

modules <- unique(pair_data$Module)
for (mod in modules) {cat("绘制模块:", mod, "\n")
  p <- plot_module_violin(mod, pair_data)
    print(p)
}


## E7C
library(ggplot2);library(scales);library(data.table);library(dplyr);library(tidyr);library(gridExtra)

Module_score_add_me_dev_tissue <- fread('D:/Postdoc/Projects/Single_cell/ZWM/code/2026_mimer代码整理/Figs_final/E7_c_k_Module_score_add_me_dev_tissue.csv')
setDT(Module_score_add_me_dev_tissue)

modules <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7")
module_names <- c("Module 1", "Module 2", "Module 3", "Module 4", "Module 5", "Module 6", "Module 7")

analysis_types <- c("development", "ME")

cat("原始development列的取值分布:\n")
print(table(Module_score_add_me_dev_tissue$development))

cat("\n原始ME列的取值分布:\n")
print(table(Module_score_add_me_dev_tissue$ME))

type_defs <- list(
  development = list(
    type_col = "development", all_types = c("Nor", "Hyp", "MiD", "MoD", "SD&CA", "ICA", "MCA"),nor_type = "Nor" ),
  ME = list(
    type_col = "ME", all_types = c("Nor-ME", "Hyp-ME", "MiD-ME", "MoD-ME", "SD&CA-ME", "ICA-ME", "MCA-ME"),nor_type = "Nor-ME" )
)

for (type in names(type_defs)) {
  if (type == "ME") {
    type_defs[[type]]$filtered_types <- type_defs[[type]]$all_types[!type_defs[[type]]$all_types %in% c("unknown", "MCA-ME")]
  } else { type_defs[[type]]$filtered_types <- type_defs[[type]]$all_types[!type_defs[[type]]$all_types %in% c("unknown", "MCA")]  }
}

common_types <- c("Nor", "Hyp", "MiD", "MoD", "SD&CA", "ICA")

threshold_percentile <- 0.95  
threshold_name <- "Top 5%"

all_modules_results <- list()

for (module_idx in seq_along(modules)) {
  module <- modules[module_idx]
  module_name <- module_names[module_idx]
  
  cat(sprintf("\n=== 分析模块: %s (%s) ===\n", module, module_name))
  data_copy <- copy(Module_score_add_me_dev_tissue)
  module_scores <- data_copy[[module]]
  threshold_value <- quantile(module_scores, threshold_percentile, na.rm = TRUE)
  
  cat(sprintf("模块 %s 的top5%%阈值: %.6f\n", module, threshold_value))
  cat(sprintf("模块 %s 分数范围: %.6f 到 %.6f\n", module,min(module_scores, na.rm = TRUE), max(module_scores, na.rm = TRUE)))
  
  data_copy[, top5_status := ifelse(get(module) >= threshold_value, "Top5", "Non-Top5")]
  
  module_results <- list()
  
  for (analysis_type in analysis_types) {
    cat(sprintf("\n分析类型: %s\n", analysis_type))
    
    type_info <- type_defs[[analysis_type]]
    type_col <- type_info$type_col
    filtered_types <- type_info$filtered_types
    nor_type <- type_info$nor_type
    
    base_data <- data_copy[get(type_col) %in% filtered_types]
    
    if (analysis_type == "ME") {
      base_data[, common_type := gsub("-ME", "", get(type_col))] } else { base_data[, common_type := get(type_col)]  }
    
    type_files <- base_data[,.(files = unique(file)), by = common_type ]
    setnames(type_files, "common_type", "type")
    
    type_n_files <- type_files[, .(n_files = .N), by = type]
    type_n_files$type <- factor(type_n_files$type, levels = common_types)
    type_n_files <- type_n_files[order(type)]
    
    sample_stats_list <- list()
    
    for (current_common_type in common_types) {files_with_type <- unique(base_data[common_type == current_common_type, file])
      
      if (length(files_with_type) > 0) {
       for (current_file in files_with_type) {
          type_cells <- base_data[file == current_file & common_type == current_common_type]
          top5_cells <- sum(type_cells$top5_status == "Top5", na.rm = TRUE)
          total_cells <- nrow(type_cells)
          
          prop <- ifelse(total_cells > 0, top5_cells / total_cells, 0)
          
          sample_stats_list[[paste(current_common_type, current_file, sep="_")]] <- data.frame(
            file = current_file,common_type = current_common_type,top5_cells = top5_cells,total_cells = total_cells,prop = prop,
            threshold = threshold_name, module = module, analysis_type = analysis_type )  }
      } else {     cat(sprintf("  警告: 通用类型 %s 在%s数据中未找到\n", current_common_type, analysis_type))    }   }
    
    if (length(sample_stats_list) > 0) {     sample_stats <- rbindlist(sample_stats_list)  } else {
      sample_stats <- data.table( file = character(), common_type = factor(levels = common_types),
        top5_cells = integer(), total_cells = integer(),prop = numeric(),threshold = character(),module = character(), analysis_type = character())}
    
    sample_stats$common_type <- factor(sample_stats$common_type, levels = common_types)
    
    mean_stats <- sample_stats[, 
                               .(mean_prop = mean(prop, na.rm = TRUE),
                                 se_prop = sd(prop, na.rm = TRUE) / sqrt(.N),
                                 median_prop = median(prop, na.rm = TRUE),
                                 sd_prop = sd(prop, na.rm = TRUE),
                                 n_samples = .N,
                                 n_nonzero = sum(prop > 0)), 
                               by = common_type
    ]
    
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
    
    mean_stats[, zero_percent := ifelse(n_samples > 0, 
                                        round((n_samples - n_nonzero) / n_samples * 100, 1), 
                                        0)]
    mean_stats[, threshold := threshold_name]
    mean_stats[, module := module]  
    mean_stats[, analysis_type := analysis_type]
    
    mean_stats$common_type <- factor(mean_stats$common_type, levels = common_types)
    mean_stats <- mean_stats[order(common_type)]
    
    nor_data <- sample_stats[common_type == "Nor", prop]
    n_nor <- length(nor_data)
    
    comparison_results <- list()
    
    for (current_common_type in common_types) {
      if (current_common_type != "Nor") {
       
        current_data <- sample_stats[common_type == current_common_type, prop]
        n_current <- length(current_data)
        
        if (n_nor >= 2 && n_current >= 2) {
          test_result <- wilcox.test(current_data, nor_data, exact = FALSE)
          p_value <- test_result$p.value
        } else {
          p_value <- NA
        }
        
        significance <- ifelse(is.na(p_value), "NA",
                               ifelse(p_value < 0.001, "***",
                                      ifelse(p_value < 0.01, "**",
                                             ifelse(p_value < 0.05, "*", "ns"))))
        
        comparison_results[[current_common_type]] <- list(
          p_value = p_value,
          significance = significance,
          n_nor = n_nor,
          n_current = n_current
        )
      }
    }
    
    module_results[[analysis_type]] <- list(
      sample_stats = sample_stats,
      mean_stats = mean_stats,
      comparison_results = comparison_results,
      n_nor = n_nor
    )
  }
  
  all_modules_results[[module]] <- list(
    development = module_results[["development"]],
    ME = module_results[["ME"]],
    module_name = module_name
  )
}

plot_list <- list()

base_title <- "Top5% Spot Proportion Across Types"
x_label <- "Type"

for (module_idx in seq_along(modules)) {
  module <- modules[module_idx]
  module_name <- module_names[module_idx]
  
  results <- all_modules_results[[module]]
  
  if (!is.null(results)) {
    dev_mean_stats <- results$development$mean_stats
    me_mean_stats <- results$ME$mean_stats
    
    dev_mean_stats$type_with_n <- paste0(dev_mean_stats$common_type, " (n=", dev_mean_stats$n_samples, ")")
    me_mean_stats$type_with_n <- paste0(me_mean_stats$common_type, " (n=", me_mean_stats$n_samples, ")")
    
    combined_data <- rbindlist(list(dev_mean_stats, me_mean_stats), fill = TRUE)
    combined_data$analysis_type <- factor(combined_data$analysis_type, levels = c("development", "ME"))
    
    combined_data$common_type <- factor(combined_data$common_type, levels = common_types)
    
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
    
    annotation_data <- data.frame()
    
    for (current_type in common_types) {
      if (current_type != "Nor") {
        comp_result <- results$development$comparison_results[[current_type]]
        
        if (!is.null(comp_result) && !is.na(comp_result$p_value)) {
          type_pos <- which(common_types == current_type)
          
          type_mean <- dev_mean_stats[common_type == current_type, mean_prop]
          type_se <- dev_mean_stats[common_type == current_type, se_prop]
          
          if (length(type_mean) > 0 && length(type_se) > 0) {
            y_pos <- type_mean + type_se
            
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
    
    for (current_type in common_types) {
      if (current_type != "Nor") {
        comp_result <- results$ME$comparison_results[[current_type]]
        
        if (!is.null(comp_result) && !is.na(comp_result$p_value)) {
         
          type_pos <- which(common_types == current_type)
          
          
          type_mean <- me_mean_stats[common_type == current_type, mean_prop]
          type_se <- me_mean_stats[common_type == current_type, se_prop]
          
          if (length(type_mean) > 0 && length(type_se) > 0) {
            y_pos <- type_mean + type_se
            
           
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
    
    y_max_dev <- max(dev_mean_stats$mean_prop + dev_mean_stats$se_prop, na.rm = TRUE, default = 0)
    y_max_me <- max(me_mean_stats$mean_prop + me_mean_stats$se_prop, na.rm = TRUE, default = 0)
    y_max <- max(y_max_dev, y_max_me, annotation_data$y_pos, na.rm = TRUE, default = 0)
    y_limit <- ifelse(y_max > 0, y_max * 1.15, 0.1)
    
    combined_data$x_pos <- as.numeric(combined_data$common_type)
    
   
    combined_data$x_pos_adjusted <- combined_data$x_pos
    combined_data$x_pos_adjusted[combined_data$analysis_type == "development"] <- 
      combined_data$x_pos[combined_data$analysis_type == "development"] - 0.1
    combined_data$x_pos_adjusted[combined_data$analysis_type == "ME"] <- 
      combined_data$x_pos[combined_data$analysis_type == "ME"] + 0.1
    
    
    if (nrow(annotation_data) > 0) {
      annotation_data$x_pos_adjusted <- annotation_data$x_pos
      annotation_data$x_pos_adjusted[annotation_data$analysis_type == "development"] <- 
        annotation_data$x_pos[annotation_data$analysis_type == "development"] - 0.1
      annotation_data$x_pos_adjusted[annotation_data$analysis_type == "ME"] <- 
        annotation_data$x_pos[annotation_data$analysis_type == "ME"] + 0.1
    }
    
    
    p <- ggplot(combined_data, aes(x = x_pos_adjusted, y = mean_prop, group = analysis_type, color = analysis_type)) +
      
      geom_point(size = 3) +
      
      geom_line(linewidth = 1) +
      
      geom_errorbar(
        aes(ymin = mean_prop - se_prop, ymax = mean_prop + se_prop),
        width = 0.2,
        linewidth = 0.8
      )
    
    
    if (nrow(annotation_data) > 0) {
      
      sig_annotations <- annotation_data[annotation_data$significance %in% c("*", "**", "***"), ]
      
      if (nrow(sig_annotations) > 0) {
        
        y_offset <- max(combined_data$mean_prop + combined_data$se_prop, na.rm = TRUE) * 0.1
        sig_annotations$y_pos <- sig_annotations$y_pos + y_offset
        
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
    
    p <- p +
      scale_x_continuous(
        breaks = 1:length(common_types),
        labels = type_labels,
        limits = c(0.5, length(common_types) + 0.5)
      ) +
      scale_y_continuous(
        labels = percent_format(accuracy = 1),
        limits = c(0, y_limit),
        expand = expansion(mult = c(0.05, 0.2))  # 为显著性标记留出更多空间
      ) +
      scale_color_discrete(
        name = "Analysis Type",
        labels = c("development", "ME")
      ) +
      labs(
        title = paste(module_name, "\n", base_title),
        x = x_label, 
        y = "Proportion of Top5% Spots"
      ) +
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

if (length(plot_list) > 0) {
  n_cols <- 3
  n_rows <- ceiling(length(plot_list) / n_cols)
  
  combined_plot <- grid.arrange(
    grobs = plot_list,
    ncol = n_cols,
    nrow = n_rows,
    top = "Combined Analysis - All Modules (Top 5% Threshold)"
  )
  
  print(combined_plot)
}

cat("\n==========  ==========\n")

for (module_idx in seq_along(modules)) {
  module <- modules[module_idx]
  module_name <- module_names[module_idx]
  
  cat(sprintf("\n模块: %s (%s)\n", module, module_name))
  
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
  
    cat("===  ===\n")
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


# E7d
library(ggplot2);library(scales);library(data.table);library(dplyr);library(tidyr);library(patchwork);library(cowplot)

mimer_median_enrich_development_ME <- fread('D:/Postdoc/Projects/Single_cell/ZWM/code/2026_mimer代码整理/Figs_final/E7_d_mimer_median_enrich_development_ME.csv')

cell_types <- c("Endo_C3_RGCC", "Mac_C2_SPP1", "FB_C3_COL1A1",  "CD8_C6_CD39","CD4_C7_OX40", "Mac_C3_C1QC", "Neutrophil", "FB_C5_PDGFRB", "Plasma" )
cell_types_cols <- gsub("_", ".", cell_types)
cell_prop_columns <- cell_types_cols

setDT(mimer_median_enrich_development_ME)

mimer_median_enrich_development_ME[, ME := as.character(ME)]
mimer_median_enrich_development_ME[, development := as.character(development)]

stage_order <- c("Nor", "Hyp", "MiD", "MoD", "SD&CA", "ICA")
me_levels <- paste0(stage_order, "-ME")
dev_levels <- stage_order

class_colors <- c("ME" = "#E41A1C", "Development" = "#377EB8")

plot_list <- list()

# --- Helper Function for Statistics and Comparisons ---
calculate_stats_and_comparisons_cell_prop <- function(data_subset, stage_col, levels_to_use, column_name, class_name) {
  # Filter data based on the specified column and levels
  filtered_data <- data_subset[get(stage_col) %in% levels_to_use]
  filtered_data[, type := get(stage_col)]
  filtered_data$type <- factor(filtered_data$type, levels = levels_to_use)
  
  if(!column_name %in% names(filtered_data)) {
    stop(paste("Column", column_name, "not found in data"))
  }
  
  file_means <- filtered_data[
    , .(
      mean_value = mean(get(column_name), na.rm = TRUE)
    ), 
    by = .(file, type)
  ]
  
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
  
  if(class_name == "ME") {
    nor_stage <- "Nor-ME"
  } else {
    nor_stage <- "Nor"
  }
  
  nor_data <- file_means[type == nor_stage, mean_value]
  
  comp_results <- list()
  
  for(stage in levels_to_use) {
    if(stage != nor_stage) {
      stage_data <- file_means[type == stage, mean_value]
      
      n_nor <- length(nor_data)
      n_stage <- length(stage_data)
      
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
    file_means = file_means  ))
}

create_x_labels <- function(me_stats, dev_stats, stage_order, me_levels, dev_levels) {
  x_labels <- character(length(stage_order))
  
  for (i in seq_along(stage_order)) {
    stage <- stage_order[i]
    
    me_stage_name <- paste0(stage, "-ME")
    me_n <- if (me_stage_name %in% me_stats$type) {
      me_stats[type == me_stage_name, n_samples]
    } else {
      0
    }
    
    dev_n <- if (stage %in% dev_stats$type) {
      dev_stats[type == stage, n_samples]
    } else {
      0
    }
    
    if (length(me_n) == 1 && length(dev_n) == 1) {
      x_labels[i] <- sprintf("%s\nME:n=%d\nDev:n=%d", stage, me_n, dev_n)
    } else {
      x_labels[i] <- stage  # 如果获取不到样本数，只显示阶段名称
    }
  }
  
  return(x_labels)
}

for (column_name in cell_prop_columns) {
  friendly_name <- column_name
  cat(sprintf("\nProcessing column: %s\n", friendly_name))
  
  me_result <- calculate_stats_and_comparisons_cell_prop(
    mimer_median_enrich_development_ME, "ME", me_levels, column_name, "ME"
  )
  me_stats <- me_result$stats
  me_comps <- me_result$comparisons
  me_file_means <- me_result$file_means
  
  dev_result <- calculate_stats_and_comparisons_cell_prop(
    mimer_median_enrich_development_ME, "development", dev_levels, column_name, "Development"
  )
  dev_stats <- dev_result$stats
  dev_comps <- dev_result$comparisons
  dev_file_means <- dev_result$file_means
  
  combined_stats <- rbindlist(list(me_stats, dev_stats))
  combined_stats[, base_stage := gsub("-ME$", "", type)] # Remove -ME suffix
  combined_stats$base_stage <- factor(combined_stats$base_stage, levels = stage_order)
  
  x_labels <- create_x_labels(me_stats, dev_stats, stage_order, me_levels, dev_levels)
  names(x_labels) <- stage_order
  
  annotation_list <- list()
  
  for(stage_name in names(me_comps)) {
    comp <- me_comps[[stage_name]]
    if(!is.na(comp$p.value)) {
      stage_pos_in_order <- which(stage_order == comp$base_stage)
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
  
  for(stage_name in names(dev_comps)) {
    comp <- dev_comps[[stage_name]]
    if(!is.na(comp$p.value)) {
      stage_pos_in_order <- which(stage_order == comp$base_stage)
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
  
  if (nrow(annotation_data) > 0) {
    annotation_data[, x_offset := ifelse(class == "ME", -0.2, 0.2)]
    annotation_data[, x_pos := x_pos_index + x_offset]
    
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

if (length(plot_list) > 0) {
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
  
  options(repr.plot.width = 18, repr.plot.height = 12)
  
print(final_plot)
} else {
  cat("No plots generated. Please check your data and target columns.\n")
}

