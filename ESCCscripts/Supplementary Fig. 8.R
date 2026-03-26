#### S8A
# S8a
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(data.table)

setwd('C:/Users/Ron/Desktop/Figs/S8/')
# 加载数据
meta_data_1 <- fread('C:/Users/Ron/Desktop/Figs/S8/S8_a_mimer_5types_above_median_split_tissue_filter_ME_dev.metadata_part1.csv')
meta_data_2 <- fread('C:/Users/Ron/Desktop/Figs/S8/S8_a_mimer_5types_above_median_split_tissue_filter_ME_dev.metadata_part2.csv')
meta_data <- rbind(meta_data_1,meta_data_2)

# 定义5类细胞类型及其对应的above_median列
cell_types <- list(
  "CD8_C6_CD39" = "CD8_C6_CD39_above_median",
  "CD4_C7_OX40" = "CD4_C7_OX40_above_median", 
  "Mac_C2_SPP1" = "Mac_C2_SPP1_above_median",
  "FB_C3_COL1A1" = "FB_C3_COL1A1_above_median",
  "Endo_C3_RGCC" = "Endo_C3_RGCC_above_median"
)

# 定义分组顺序
cell_types_order <- c("CD8_C6_CD39", "CD4_C7_OX40", "Mac_C2_SPP1", "FB_C3_COL1A1", "Endo_C3_RGCC")
cell_order_mimer <- paste0(cell_types_order, "_MIMER")
cell_order_other <- paste0(cell_types_order, "_other")

# 定义颜色方案
cell_type_colors <- c(
  "CD8_C6_CD39" = "#E31A1C",  # 红色
  "CD4_C7_OX40" = "#1F78B4",  # 蓝色
  "Mac_C2_SPP1" = "#33A02C",  # 绿色
  "FB_C3_COL1A1" = "#6A3D9A",  # 紫色
  "Endo_C3_RGCC" = "#FF7F00"   # 橙色
)

group_colors <- c(
  "MIMER" = "#bd2628",  # 红色
  "other" = "#cccccc"   # 灰色
)

# 1. 计算热图矩阵的函数
calculate_heatmap_matrix <- function(filter_col) {
  # 创建5类细胞类型的热图矩阵
  cell_type_groups <- list()
  
  for(cell_name in names(cell_types)) {
    col_name <- cell_types[[cell_name]]
    
    cell_data <- meta_data %>%
      filter(!get(filter_col) %in% c('unknown')) %>%
      mutate(group = case_when(
        mimer == 'MIMER' & get(col_name) == TRUE ~ paste0(cell_name, "_MIMER"),
        mimer == 'Others' & get(col_name) == TRUE ~ paste0(cell_name, "_other"),
        TRUE ~ NA_character_
      )) %>%
      filter(!is.na(group)) %>%
      group_by(group) %>%
      summarise(across(
        starts_with("signature."),
        mean,
        na.rm = TRUE
      ))
    
    cell_type_groups[[cell_name]] <- cell_data
  }
  
  combined_cell_data <- do.call(rbind, cell_type_groups)
  heatmap_matrix <- as.matrix(combined_cell_data[, -1])
  rownames(heatmap_matrix) <- combined_cell_data$group
  
  # 按分组顺序重新排列行
  cell_order <- c(cell_order_mimer, cell_order_other)
  cell_order <- cell_order[cell_order %in% rownames(heatmap_matrix)]
  heatmap_matrix <- heatmap_matrix[cell_order, ]
  
  return(heatmap_matrix)
}

# 2. 计算两个矩阵
me_matrix <- calculate_heatmap_matrix("ME")
dev_matrix <- calculate_heatmap_matrix("development")

# 3. 基于ME数据计算signature的顺序
# 对ME矩阵的列（signature）进行聚类
signature_hclust <- hclust(dist(t(me_matrix)))
signature_order <- colnames(me_matrix)[signature_hclust$order]

# 4. 计算全局最大最小值
global_min <- min(min(me_matrix, na.rm = TRUE), min(dev_matrix, na.rm = TRUE))
global_max <- max(max(me_matrix, na.rm = TRUE), max(dev_matrix, na.rm = TRUE))

# 5. 设置统一的颜色方案（以0为中点）
getPalette = colorRampPalette(rev(brewer.pal(9, "RdBu")))
colors = getPalette(100)  # 100个颜色

# 创建以0为中心的breakpoints序列
if(global_min < 0 && global_max > 0) {
  # 如果最小值小于0且最大值大于0，将范围分成两段
  # 从最小值到0分成50份
  bk1 <- seq(global_min, 0, length.out = 50)
  # 从0到最大值分成50份
  bk2 <- seq(0, global_max, length.out = 51)  # 从0开始，避免重复0
  # 合并，去掉一个0
  bk <- c(bk1, bk2[-1])
} else if(global_min >= 0) {
  # 如果最小值大于等于0，所有值都非负
  bk <- seq(global_min, global_max, length.out = 101)
} else if(global_max <= 0) {
  # 如果最大值小于等于0，所有值都非正
  bk <- seq(global_min, global_max, length.out = 101)
}

# 6. 创建行注释的函数
create_row_annotation <- function(heatmap_matrix) {
  cell_type_labels <- sub("_(MIMER|other)$", "", rownames(heatmap_matrix))
  group_labels <- ifelse(grepl("_MIMER$", rownames(heatmap_matrix)), "MIMER", "other")
  
  row_annotation <- data.frame(
    Cell_Type = factor(cell_type_labels, levels = cell_types_order),
    Group = factor(group_labels, levels = c("MIMER", "other"))
  )
  rownames(row_annotation) <- rownames(heatmap_matrix)
  
  return(row_annotation)
}

# 7. 绘制热图的函数
plot_heatmap <- function(heatmap_matrix, output_filename, plot_title, row_annotation) {
  # 按signature顺序重新排列矩阵的列
  available_signatures <- colnames(heatmap_matrix)
  signature_order_filtered <- signature_order[signature_order %in% available_signatures]
  
  # 如果signature_order_filtered不包含所有可用的signature，添加缺失的
  missing_signatures <- setdiff(available_signatures, signature_order_filtered)
  if(length(missing_signatures) > 0) {
    signature_order_filtered <- c(signature_order_filtered, missing_signatures)
  }
  
  # 按顺序重新排列矩阵的列
  heatmap_matrix_ordered <- heatmap_matrix[, signature_order_filtered, drop = FALSE]
  
  # 设置绘图尺寸
  options(repr.plot.width = 10, repr.plot.height = 14)
  
  # 绘制热图
  pheatmap(t(heatmap_matrix_ordered),
           filename = output_filename,
           annotation_col = row_annotation,
           main = plot_title,
           
           # 颜色与图例
           color = colors,
           border_color = NA,
           breaks = bk,  # 使用以0为中心的breakpoints
           
           # 聚类与树状图
           cluster_rows = FALSE,  # 不聚类行
           cluster_cols = FALSE,  # 不聚类列
           treeheight_row = 0,
           treeheight_col = 0,
           
           # 单元格与字体
           cellwidth = 25,
           cellheight = 20,
           fontsize_row = 10,
           fontsize_col = 12,
           angle_col = 45,
           
           # 图例与标签
           show_rownames = TRUE,
           show_colnames = TRUE,
           
           # 添加颜色区分
           annotation_colors = list(
             Cell_Type = cell_type_colors,
             Group = group_colors
           ),
           
           # 确保注释颜色正确显示
           annotation_names_col = TRUE
  )
  
  cat("已生成热图，保存为:", output_filename, "\n")
}

# 8. 分别绘制两个热图
# 绘制ME热图
me_row_annotation <- create_row_annotation(me_matrix)
plot_heatmap(me_matrix, 
             'cell_type_MIMER_vs_other_comparison_ME.pdf',
             "Signature Scores: MIMER vs Other (by cell type in ME)",
             me_row_annotation)

# 绘制development热图
dev_row_annotation <- create_row_annotation(dev_matrix)
plot_heatmap(dev_matrix,
             'cell_type_MIMER_vs_other_comparison_epitu.pdf',
             "Signature Scores: MIMER vs Other (by cell type in Epi/Tumor)",
             dev_row_annotation)


## S8b
# ============================================================
# Refactored CODEX functional scoring pipeline for MIMER analysis
# ------------------------------------------------------------
# Main fixes compared with the original draft:
# 1) Explicit CODEX subtype -> MIMER supertype mapping
# 2) Collapse mapped CODEX subtypes into 5 supertype classes
# 3) Standardize markers within each patient
# 4) Calculate pathway scores only when marker coverage is sufficient
# 5) Aggregate to patient level before paired testing
# 6) Compute paired Wilcoxon tests independently and apply BH-FDR
# 7) Export Excel tables and plot only FDR-corrected significant results
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(openxlsx)
  library(stringr)
  library(tibble)
  library(scales)
  library(qs) ##zxh
})

# -----------------------------
# 0. User-configurable section
# -----------------------------
input_codex_rdata <- "/lustre/home/xhzh/scSpatial/xiaodu/codes/zxh/MIMER/codex_datasets.Rdata"
input_classification_rdata <- "/lustre/home/xhzh/scSpatial/xiaodu/codes/zxh/MIMER/classification_results.Rdata"

output_dir <- "/lustre/home/xhzh/scSpatial/xiaodu/codes/zxh/MIMER"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Score calculation parameters
min_marker_fraction <- 0.67   # at least 67% of markers must exist in the panel
min_marker_abs <- 1           # but never require fewer than 1 marker
min_pairs_for_test <- 3       # minimum number of paired patients for formal test
fdr_alpha <- 0.05             # plot only FDR-significant pathways

# ----------------------------------------------------
# 1. Explicit CODEX subtype -> MIMER supertype mapping
# ----------------------------------------------------
# Final analysis keeps exactly 5 modified CODEX classes.
subtype_to_supertype <- tibble::tribble(
  ~codex_subtype,       ~mimer_supertype,
  "CD4T-aTreg",         "CD4_C7_0X40",
  "Tcells.cc",          "CD4_C7_0X40",
  "PDCD1+Tex",          "CD8_C6_CD39",
  "LAG3+Tex",           "CD8_C6_CD39",
  "PDCD1+Tex2",         "CD8_C6_CD39",
  "PLVAP+Endo",         "Endo_C3_RGCC",
  "CollagenIV+Endo",    "Endo_C3_RGCC",
  "COL1A1+CAF",         "FB_C3_COL1A1",
  "POSTN+CAF",          "FB_C3_COL1A1",
  "ISG15+Mac",          "Mac_C2_SPP1",
  "SPP1+Mac",           "Mac_C2_SPP1"
)

# ---------------------------------------
# 2. Pathway definition at supertype level
# ---------------------------------------
pathway_defs <- list(
  "CD8_C6_CD39" = list(
    Tex_ExhaustionScore   = c("PD-1", "TOX", "LAG3"),
    Tex_TCF1              = c("TCF-1"),
    Tex_EffectorScore     = c("GZMB", "IFNG", "CD57"),
    Tex_ProliferationScore= c("Ki67", "PCNA", "HistoneH3-pSer28")
  ),
  "CD4_C7_0X40" = list(
    Treg_CD39              = c("CD39"),
    Treg_ActivationScore   = c("ICOS", "OX40"),
    Treg_ProliferationScore= c("Ki67", "PCNA", "HistoneH3-pSer28")
  ),
  "Mac_C2_SPP1" = list(
    TAM_ImmunosuppressionScore = c("PD-L1", "IDO1", "VISTA", "CD39"),
    TAM_M2LikeScore            = c("CD163", "CD206"),
    TAM_AntigenPresentationScore = c("HLA-DR")
  ),
  "Endo_C3_RGCC" = list(
    Endo_PermeabilityScore     = c("PLVAP-PV-1", "Caveolin"),
    Endo_ImmuneRegulationScore = c("CD39", "PD-L1", "IDO1", "HLA-DR"),
    Endo_ProliferationScore    = c("Ki67", "PCNA", "HistoneH3-pSer28")
  ),
  "FB_C3_COL1A1" = list(
    CAF_SMA                   = c("SMA"),
    CAF_ECMRemodelingScore    = c("Periostin", "COL-1", "CollagenIV"),
    CAF_ImmuneRegulationScore = c("PD-L1", "IDO1"),
    Stromal_ISG15_ResponseScore = c("ISG15")
  )
)

# ---------------------------------------
# 3. Helper functions
# ---------------------------------------

# Compatible extraction for Seurat v4 / v5
get_expr_matrix <- function(obj, assay = "Akoya") {
  out <- tryCatch(
    Seurat::GetAssayData(obj, assay = assay, layer = "data"),
    error = function(e) {
      Seurat::GetAssayData(obj, assay = assay, slot = "data")
    }
  )
  return(out)
}

# Safe z-score inside one patient
safe_zscore <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  } else {
    return((x - m) / s)
  }
}

# Required marker count for one pathway
required_marker_n <- function(total_markers,
                              min_fraction = min_marker_fraction,
                              min_abs = min_marker_abs) {
  max(min_abs, ceiling(total_markers * min_fraction))
}

# Paired Wilcoxon at patient level
run_paired_tests <- function(patient_scores,
                             min_pairs = min_pairs_for_test) {

  split_list <- split(
    patient_scores,
    list(patient_scores$mimer_supertype, patient_scores$pathway),
    drop = TRUE
  )

  stats_list <- lapply(split_list, function(df_one) {

    wide_df <- df_one %>%
      dplyr::select(patient_id, condition, mean_score) %>%
      dplyr::distinct() %>%
      tidyr::pivot_wider(names_from = condition, values_from = mean_score)

    wide_df <- wide_df %>%
      dplyr::filter(!is.na(other) & !is.na(MIMER))

    n_pairs <- nrow(wide_df)

    if (n_pairs < min_pairs) {
      return(tibble(
        mimer_supertype = unique(df_one$mimer_supertype),
        pathway = unique(df_one$pathway),
        n_pairs = n_pairs,
        median_diff = if (n_pairs > 0) median(wide_df$MIMER - wide_df$other) else NA_real_,
        mean_diff = if (n_pairs > 0) mean(wide_df$MIMER - wide_df$other) else NA_real_,
        p_value = NA_real_,
        test_note = paste0("Skipped: n_pairs < ", min_pairs)
      ))
    }

    wt <- suppressWarnings(stats::wilcox.test(
      x = wide_df$MIMER,
      y = wide_df$other,
      paired = TRUE,
      exact = FALSE,
      conf.int = FALSE
    ))

    tibble(
      mimer_supertype = unique(df_one$mimer_supertype),
      pathway = unique(df_one$pathway),
      n_pairs = n_pairs,
      median_diff = median(wide_df$MIMER - wide_df$other),
      mean_diff = mean(wide_df$MIMER - wide_df$other),
      p_value = wt$p.value,
      test_note = "Paired Wilcoxon test"
    )
  })

  out <- dplyr::bind_rows(stats_list) %>%
    dplyr::mutate(
      p_adj = stats::p.adjust(p_value, method = "BH"),
      sig_fdr = !is.na(p_adj) & p_adj < fdr_alpha
    ) %>%
    dplyr::arrange(p_adj, p_value)

  return(out)
}

# Plot only FDR-significant pathways
plot_significant_results <- function(patient_scores, paired_stats, output_dir) {

  sig_stats <- paired_stats %>%
    dplyr::filter(sig_fdr)

  if (nrow(sig_stats) == 0) {
    message("No FDR-significant pathways detected. No plot exported.")
    return(invisible(NULL))
  }

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

    p <- ggplot(df_sub, aes(x = condition, y = mean_score, color = condition)) +
      geom_line(aes(group = patient_id), color = "grey80", linewidth = 0.5) +
      geom_boxplot(aes(fill = condition), width = 0.35, alpha = 0.18, outlier.shape = NA) +
      geom_point(size = 2.2, alpha = 0.85) +
      facet_wrap(~ pathway_label, scales = "free_y") +
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
        subtitle = "Only pathways passing BH-FDR are displayed",
        x = NULL,
        y = "Patient-level mean pathway score (within-patient z-scored markers)"
      )

    ggplot2::ggsave(
      filename = file.path(output_dir, paste0("plot_", this_super, "_FDR_sig.pdf")),
      plot = p,
      width = 10,
      height = 6.5
    )
  }
}

# ---------------------------------------
# 4. Load data
# ---------------------------------------
load(input_codex_rdata)
load(input_classification_rdata)
codex = qs::qread('/lustre/home/xhzh/scSpatial/xiaodu/codes/zxh/MIMER/final_codex.qs')

# codex object selection: adapt here if your object name is different
if (!exists("codex")) {
  stop("Object `codex` was not found after loading input_codex_rdata.")
}
if (!exists("classification_results")) {
  stop("Object `classification_results` was not found after loading input_classification_rdata.")
}

# ---------------------------------------
# 5. Prepare metadata and collapse subtype
# ---------------------------------------
codex@meta.data <- codex@meta.data %>% rename_with(~ "cell_id_bk", any_of("cell_id")) ##zxh
meta_tbl <- codex@meta.data %>%
  tibble::rownames_to_column("cell_id")

#required_meta_cols <- c("Patient", "subCelltype.x")
required_meta_cols <- c("Patient", "subCelltype") ##zxh
missing_meta_cols <- setdiff(required_meta_cols, colnames(meta_tbl))
if (length(missing_meta_cols) > 0) {
  stop("Missing metadata columns in codex@meta.data: ",
       paste(missing_meta_cols, collapse = ", "))
}

required_cls_cols <- c("cell_id", "final_classification")
missing_cls_cols <- setdiff(required_cls_cols, colnames(classification_results))
if (length(missing_cls_cols) > 0) {
  stop("Missing columns in classification_results: ",
       paste(missing_cls_cols, collapse = ", "))
}

mimer_cell_ids <- classification_results %>%
  dplyr::filter(final_classification == "MIMER_group") %>%
  dplyr::pull(cell_id) %>%
  unique()

meta_tbl <- meta_tbl %>%
  dplyr::transmute(
    cell_id = cell_id,
    patient_id = Patient,
    #codex_subtype = subCelltype.x
    codex_subtype = subCelltype ##zxh
  ) %>%
  dplyr::left_join(subtype_to_supertype, by = "codex_subtype") %>%
  dplyr::filter(!is.na(mimer_supertype)) %>%
  dplyr::mutate(
    analysis_subtype = mimer_supertype,   # collapse CODEX subtype to 5 final classes
    condition = ifelse(cell_id %in% mimer_cell_ids, "MIMER", "other")
  )

if (nrow(meta_tbl) == 0) {
  stop("No mapped cells remained after subtype_to_supertype filtering.")
}

# ---------------------------------------
# 6. Extract expression and standardize markers within patient
# ---------------------------------------
expr_mat <- get_expr_matrix(codex, assay = "Akoya")
expr_mat <- as.matrix(expr_mat)

all_markers <- unique(unlist(unname(unlist(pathway_defs, recursive = FALSE))))
available_markers <- intersect(all_markers, rownames(expr_mat))

if (length(available_markers) == 0) {
  stop("None of the pathway markers are present in the Akoya assay.")
}

cell_keep <- intersect(meta_tbl$cell_id, colnames(expr_mat))
meta_tbl <- meta_tbl %>% dplyr::filter(cell_id %in% cell_keep)
expr_mat <- expr_mat[available_markers, meta_tbl$cell_id, drop = FALSE]

expr_df <- as.data.frame(t(expr_mat))
expr_df$cell_id <- rownames(expr_df)

score_input <- meta_tbl %>%
  dplyr::left_join(expr_df, by = "cell_id")

# Standardize every marker inside each patient
score_input <- score_input %>%
  dplyr::group_by(patient_id) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(available_markers), safe_zscore)) %>%
  dplyr::ungroup()

# ---------------------------------------
# 7. Cell-level pathway score calculation
# ---------------------------------------
coverage_log <- list()
cell_score_list <- list()

for (this_super in names(pathway_defs)) {

  super_cells <- score_input %>%
    dplyr::filter(mimer_supertype == this_super)

  if (nrow(super_cells) == 0) next

  pathways_this_super <- pathway_defs[[this_super]]

  for (this_pathway in names(pathways_this_super)) {

    marker_set <- pathways_this_super[[this_pathway]]
    markers_present <- intersect(marker_set, colnames(super_cells))
    markers_required <- required_marker_n(length(marker_set))

    coverage_log[[paste(this_super, this_pathway, sep = "__")]] <- tibble(
      mimer_supertype = this_super,
      pathway = this_pathway,
      defined_markers_n = length(marker_set),
      required_markers_n = markers_required,
      available_markers_n = length(markers_present),
      available_markers = paste(markers_present, collapse = "; "),
      coverage_pass = length(markers_present) >= markers_required
    )

    if (length(markers_present) < markers_required) {
      next
    }

    tmp_df <- super_cells %>%
      dplyr::transmute(
        cell_id,
        patient_id,
        condition,
        mimer_supertype,
        pathway = this_pathway,
        pathway_score = rowMeans(dplyr::across(dplyr::all_of(markers_present)), na.rm = TRUE),
        marker_count_used = length(markers_present)
      )

    cell_score_list[[paste(this_super, this_pathway, sep = "__")]] <- tmp_df
  }
}

marker_coverage <- dplyr::bind_rows(coverage_log)

if (length(cell_score_list) == 0) {
  stop("No pathway could be scored after marker coverage filtering.")
}

cell_scores <- dplyr::bind_rows(cell_score_list)

# ---------------------------------------
# 8. Aggregate to patient level
# ---------------------------------------
patient_scores <- cell_scores %>%
  dplyr::group_by(patient_id, condition, mimer_supertype, pathway) %>%
  dplyr::summarise(
    mean_score = mean(pathway_score, na.rm = TRUE),
    median_score = median(pathway_score, na.rm = TRUE),
    n_cells = dplyr::n(),
    marker_count_used = dplyr::first(marker_count_used),
    .groups = "drop"
  ) %>%
  dplyr::arrange(mimer_supertype, pathway, patient_id, condition)

# Keep only truly paired patients for statistical testing
paired_patient_scores <- patient_scores %>%
  dplyr::group_by(patient_id, mimer_supertype, pathway) %>%
  dplyr::filter(dplyr::n_distinct(condition) == 2) %>%
  dplyr::ungroup()

# ---------------------------------------
# 9. Independent paired test + BH-FDR
# ---------------------------------------
paired_stats <- run_paired_tests(
  patient_scores = paired_patient_scores,
  min_pairs = min_pairs_for_test
)

# ---------------------------------------
# 10. Export tables to Excel
# ---------------------------------------
xlsx_path <- file.path(output_dir, "codex_function_results.xlsx")

wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "mapping_table")
openxlsx::writeData(wb, "mapping_table", subtype_to_supertype)

openxlsx::addWorksheet(wb, "marker_coverage")
openxlsx::writeData(wb, "marker_coverage", marker_coverage)

openxlsx::addWorksheet(wb, "patient_level_scores")
openxlsx::writeData(wb, "patient_level_scores", patient_scores)

openxlsx::addWorksheet(wb, "paired_stats")
openxlsx::writeData(wb, "paired_stats", paired_stats)

openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)

# ---------------------------------------
# 11. Export plots (only FDR-significant)
# ---------------------------------------
plot_significant_results(
  patient_scores = paired_patient_scores,
  paired_stats = paired_stats,
  output_dir = output_dir
)

# ---------------------------------------
# 12. Save RDS copies for reproducibility
# ---------------------------------------
saveRDS(subtype_to_supertype, file.path(output_dir, "subtype_to_supertype.rds"))
saveRDS(marker_coverage, file.path(output_dir, "marker_coverage.rds"))
saveRDS(patient_scores, file.path(output_dir, "patient_level_scores.rds"))
saveRDS(paired_stats, file.path(output_dir, "paired_stats.rds"))

message("Refactored pipeline completed.")
message("Excel: ", xlsx_path)

