library(ggrastr)
library(ggpubr)
library(tidyverse)
load("/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhangzhichao/spatial_pathology/sub.dist.Rdata")
sub.dist.all = do.call(rbind,sub.dist)

sub.dist.all$cellID = sub.dist.all$id
load("/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhangzhichao/spatial_pathology/sp_meta_info.Rdata")
# sub.dist.all %>% left_join(meta_info,by = "cellID") -> plot_df
sub.dist.all %>% left_join(meta,by = "cellID") -> plot_df
# table(plot_df$distinct_area)
colnames(plot_df)
# HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
plot_df$development = factor(plot_df$development,
                             levels =  c('Nor','Hyp','MiD','MoD','SD&CA','ICA','MCA'))
hallmark_cols = which(grepl("^HALLMARK",colnames(plot_df)))
head(plot_df)
colnames(plot_df)
use_cols = c("min","Cancer","development",'cellID',colnames(plot_df)[hallmark_cols[1:4]])
hallmark_paths = c("HALLMARK_ANGIOGENESIS",'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
                   'HALLMARK_KRAS_SIGNALING_UP','HALLMARK_TGF_BETA_SIGNALING',
                   'HALLMARK_MYOGENESIS','HALLMARK_TNFA_SIGNALING_VIA_NFKB',
                   'HALLMARK_PANCREAS_BETA_CELLS','HALLMARK_COAGULATION')
# 
# hallmark_paths = c("pro_fibrotic_signature",'ECM',
#                    'Anti_inflammatory','Immunosuppression',
#                    'pro_metastasis')
hallmark_list = lapply(hallmark_paths, function(path){
  ggplot(plot_df,aes(x=min,y=plot_df[[path]] ,color=Cancer) )+
    geom_point_rast(size=0.1) +
    scale_color_distiller(palette = "Spectral")+
    stat_cor(size = 2)+
    geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2",span=1, size = 0.5)+
    facet_grid(~development)+
    labs(y=path,x="Distance to Mimer (Spot)")+theme_classic() +
    theme(
      legend.key.size = unit(2, "mm"),
      panel.spacing = unit(1, "mm") ,
      axis.line = element_line(size = 0.2),
      axis.ticks = element_line(size = 0.2),
      # 全局文本大小（包括标题、坐标轴标签等）
      text = element_text(size = 5),  # 6pt
      # 坐标轴刻度文字
      axis.text = element_text(size = 6),
      # 图例文字
      legend.text = element_text(size = 5),
      # 标题文字
      #plot.title = element_text(size = 10),
      strip.text.x = element_text(size = 6),
      strip.background = element_rect(size = 0.2)
      
    )->p1
  return(p1)
})
names(hallmark_list) = hallmark_paths

library(patchwork)
wrap_plots(hallmark_list,ncol = 1)->p



ggsave(
  filename = "/realspace/project/proj_ESCC_STW_ZWM_2022_01/zhangzhichao/spatial_pathology/distinct/plot_A4_2.pdf",  # 支持 PDF/PNG/TIFF 等格式
  plot = p,
  device = "pdf",            # 保存为 PDF
  width = 210,               # A4 宽度 (mm)
  height = 297,              # A4 高度 (mm)
  units = "mm",              # 单位设为毫米
  dpi = 300                  # 分辨率（默认 300 DPI）
)

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

