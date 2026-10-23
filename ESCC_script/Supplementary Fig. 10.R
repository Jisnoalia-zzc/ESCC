#### S10a
library(ggsci)
library(dplyr)
library(ggplot2)
library(tibble)
library(patchwork)

metadata <- fread("D:/Postdoc/Projects/Single_cell/ZWM/code/2026_mimer´úÂëÕûÀí/Figs_final/Sup10_ab_scSurvival_metadata.csv")
data <- matrix(metadata$attention, ncol = 1)
set.seed(42)
kmeans_result <- kmeans(data, centers = 2, nstart = 10)
labels <- kmeans_result$cluster
cluster_centers <- kmeans_result$centers
atten_thr <- mean(cluster_centers)
cat("attention cutoff:", atten_thr, "\n")
metadata <- metadata %>%
  mutate(
    surv.risk.group = case_when(
      hazard_adj > 0 & attention > atten_thr ~ 'higher',
      hazard_adj < 0 & attention > atten_thr ~ 'lower',
      TRUE ~ 'inattentive'
    ) %>% factor(levels = c('higher', 'inattentive', 'lower'))
  )

df <- metadata %>%
  as_tibble() %>%
  group_by(Patient, subCelltype) %>%
  summarise(hazard_adj = mean(hazard_adj, na.rm = TRUE), .groups = "drop")

order_sub <- df %>%
  group_by(subCelltype) %>%
  summarise(mean_hazard = mean(hazard_adj)) %>%
  arrange(desc(mean_hazard)) %>%
  pull(subCelltype)

df$subCelltype <- factor(df$subCelltype, levels = order_sub)

df_bar <- metadata %>%
  as_tibble() %>%
  mutate(subCelltype = factor(subCelltype, levels = order_sub)) %>%
  count(subCelltype, surv.risk.group) %>%
  group_by(subCelltype) %>%
  mutate(
    proportion = n / sum(n),
    percentage = round(proportion * 100, 1)
  ) %>%
  ungroup()

risk_colors <- c(
  'higher' = '#E41A1C', 
  'lower' = '#377EB8', 
  'inattentive' = '#999999')

sci_palette <- c(
  pal_npg()(10),
  pal_lancet()(10),
  pal_jco()(10),
  pal_nejm()(8),
  pal_jama()(10),
  pal_locuszoom()(10)
)

if (length(order_sub) > length(sci_palette)) {
  sci_palette <- rep(sci_palette, 
                     ceiling(length(order_sub) / length(sci_palette)))
}
sci_palette <- sci_palette[1:length(order_sub)]

p_violin <- ggplot(df, aes(x = subCelltype, y = hazard_adj, fill = subCelltype)) +
  geom_violin(
    trim = FALSE,
    scale = "width",
    alpha = 0.8,
    linewidth = 0.3
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.6
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.6
  ) +
  scale_fill_manual(values = sci_palette) +
  labs(
    x = NULL,
    y = "Adjusted Hazard Score",
    title = "Distribution of hazard_adj by subCelltype"
  ) +
  theme_classic(base_size = 14) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p_bar <- ggplot(df_bar, aes(x = subCelltype, y = proportion, fill = surv.risk.group)) +
  geom_bar(stat = "identity", width = 1, color = "white", size = 0.3) +
  geom_text(
    aes(label = ifelse(percentage >= 5, paste0(percentage, "%"), "")),
    position = position_stack(vjust = 0.5),
    size = 2.5,
    color = "white",
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = risk_colors,
    name = "Risk Group",
    labels = c("Higher Risk", "Lower Risk", "Inattentive")
  ) +
  labs(
    x = NULL,
    y = "Percentage (%)",
    title = "Survival Risk Group Distribution by subCelltype"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.y = element_line(color = "grey80", linewidth = 0.3)
  )

combined_plot <- p_bar / p_violin +
  plot_layout(heights = c(1, 1))

options(repr.plot.width = 15, repr.plot.height = 8)
combined_plot



#### S10b
# ====================  ====================
library(ggsci)
library(dplyr)
library(ggplot2)
library(tibble)
library(patchwork)
library(rstatix)
library(ggpubr)
library(ggsignif)
library(tidyr)


# ====================  ====================
subcelltype_mapping <- c(
  "PDCD1+Tex2"   = "PDCD1+Tex",
  "LAG3+Tex"     = "PDCD1+Tex",
  "CD4T-aTreg"   = "OX4T+Tregs",
  "Tcells.cc"    = "OX4T+Tregs",
  "SPP1+Mac"     = "SPP1+Macs",
  "ISG15+Mac"    = "SPP1+Macs",
  "COL1A1+CAF"   = "COL1A1+CAFs",
  "POSTN+CAF"    = "COL1A1+CAFs",
  "PLVAP+Endo"   = "RGCC+Endos",
  "CollagenIV+Endo" = "RGCC+Endos"
)

metadata_filtered <- fread("Sup10_b_scSurvival_combined.csv")
# ===== =====
df_violin <- metadata_filtered %>%
  as_tibble() %>%
  group_by(Sample, mimer_mapped, three_group) %>%
  summarise(
    hazard_adj = mean(hazard_adj, na.rm = TRUE),
    surv.risk.group = first(surv.risk.group),
    cell_count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    three_group = factor(
      three_group,
      levels = c("Mimer+_Cell+", "Mimer+_Cell-", "Mimer-_All")
    )
  )

unique_regions <- sort(unique(df_violin$mimer_mapped))

order_sub <- c()
for (region in unique_regions) {
  order_sub <- c(
    order_sub,
    paste0(region, "_Mimer+_Cell+"),
    paste0(region, "_Mimer+_Cell-"),
    paste0(region, "_Mimer-_All")
  )
}

df_violin <- df_violin %>%
  mutate(
    region_three_group = factor(
      paste(mimer_mapped, three_group, sep = "_"),
      levels = order_sub
    )
  )

# ===== =====
df_bar <- metadata_filtered %>%
  mutate(
    region_three_group = factor(
      paste(mimer_mapped, three_group, sep = "_"),
      levels = order_sub
    )
  ) %>%
  count(region_three_group, surv.risk.group) %>%
  group_by(region_three_group) %>%
  mutate(
    proportion = n / sum(n),
    percentage = round(proportion * 100, 1)
  ) %>%
  ungroup() %>%
  mutate(
    mimer_mapped = gsub("(_Mimer\\+_Cell\\+|__Mimer\\+_Cell-|__Mimer-_All)$", "", region_three_group)
  )

# =====  =====
risk_colors <- c(
  'higher' = '#E41A1C',
  'lower' = '#377EB8',
  'inattentive' = '#999999'
)

three_group_colors <- c(
  "Mimer+_Cell+" = "#E41A1C",
  "Mimer+_Cell-" = "#377EB8",
  "Mimer-_All"   = "#4DAF4A"
)

color_vector <- c()
for (region in unique_regions) {
  color_vector[paste0(region, "_Mimer+_Cell+")] <- three_group_colors["Mimer+_Cell+"]
  color_vector[paste0(region, "_Mimer+_Cell-")] <- three_group_colors["Mimer+_Cell-"]
  color_vector[paste0(region, "_Mimer-_All")]   <- three_group_colors["Mimer-_All"]
}

# ===== =====
p_violin <- ggplot(df_violin, aes(x = region_three_group, y = hazard_adj, fill = region_three_group)) +
  geom_violin(
    trim = FALSE,
    scale = "width",
    alpha = 0.7,
    linewidth = 0.3
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.8,
    linewidth = 0.5
  ) +
  geom_jitter(
    size = 2,
    alpha = 0.7,
    shape = 21,
    stroke = 0.3,
    position = position_jitter(width = 0.15, height = 0.1),
    aes(fill = region_three_group)
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.6
  ) +
  scale_fill_manual(values = color_vector) +
  labs(
    x = NULL,
    y = "Adjusted Hazard Score",
    title = "Hazard Score by Mapped Cell Type and Three-Group Status"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold",
      size = 12
    ),
    axis.text.y = element_text(size = 12),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
  )

# =====  =====
signif_comparisons <- list()
signif_y_positions <- c()

y_range <- range(df_violin$hazard_adj, na.rm = TRUE)
y_max <- y_range[2]
y_min <- y_range[1]
y_spread <- diff(y_range)

current_y_offset <- y_max + 0.05 * y_spread

for (region in unique_regions) {
  sub_data <- df_violin %>% filter(mimer_mapped == region)
  
  group_a <- paste0(region, "_Mimer+_Cell+")
  group_b <- paste0(region, "_Mimer+_Cell-")
  data_ab <- sub_data %>% filter(three_group %in% c("Mimer+_Cell+", "Mimer+_Cell-"))
  if (n_distinct(data_ab$three_group) == 2) {
    p_val1 <- wilcox.test(hazard_adj ~ three_group, data = data_ab)$p.value
    signif_comparisons[[length(signif_comparisons) + 1]] <- c(group_a, group_b)
    signif_y_positions <- c(signif_y_positions, current_y_offset)
    current_y_offset <- current_y_offset + 0.08 * y_spread
  }
  
  group_c <- paste0(region, "_Mimer-_All")
  data_ac <- sub_data %>% filter(three_group %in% c("Mimer+_Cell+", "Mimer-_All"))
  if (n_distinct(data_ac$three_group) == 2) {
    p_val2 <- wilcox.test(hazard_adj ~ three_group, data = data_ac)$p.value
    signif_comparisons[[length(signif_comparisons) + 1]] <- c(group_a, group_c)
    signif_y_positions <- c(signif_y_positions, current_y_offset)
    current_y_offset <- current_y_offset + 0.08 * y_spread
  }
}

if (length(signif_comparisons) > 0) {
  p_violin <- p_violin +
    geom_signif(
      comparisons = signif_comparisons,
      map_signif_level = function(p) {
        if (p < 0.001) "***"
        else if (p < 0.01) "**"
        else if (p < 0.05) "*"
        else "ns"
      },
      y_position = signif_y_positions,
      tip_length = 0.02,
      textsize = 4,
      vjust = 0.5,
      color = "black"
    )
}

# =====  =====
p_bar <- ggplot(df_bar, aes(x = region_three_group, y = proportion, fill = surv.risk.group)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", size = 0.3) +
  geom_text(
    aes(label = ifelse(percentage >= 5, paste0(percentage, "%"), "")),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "white",
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = risk_colors,
    name = "Risk Group",
    labels = c("Higher Risk", "Lower Risk", "Inattentive")
  ) +
  labs(
    x = NULL,
    y = "Cell Percentage (%)",
    title = "Cell Composition by Risk Group (Mapped Types)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.y = element_line(color = "grey80", linewidth = 0.3),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold")
  )

# =====  =====
combined_plot <- p_bar / p_violin +
  plot_layout(heights = c(1, 1.5))

options(repr.plot.width = 14, repr.plot.height = 10)
print(combined_plot)

#------------------------------
#------part 2. high risk prop---
#-------------------------------

# =====  =====
df_prop <- metadata_filtered %>%
  filter(!is.na(surv.risk.group)) %>%
  group_by(Sample, mimer_mapped, three_group) %>%
  summarise(
    total_cells = n(),
    high_risk_cells = sum(surv.risk.group == "higher", na.rm = TRUE),
    prop_high_risk = high_risk_cells / total_cells,
    .groups = "drop"
  ) %>%
  mutate(
    three_group = factor(three_group,
                         levels = c("Mimer+_Cell+", "Mimer+_Cell-", "Mimer-_All"))
  )

# =====  =====
unique_regions <- sort(unique(df_prop$mimer_mapped))
order_sub <- c()
for (region in unique_regions) {
  order_sub <- c(order_sub,
                 paste0(region, "_Mimer+_Cell+"),
                 paste0(region, "_Mimer+_Cell-"),
                 paste0(region, "_Mimer-_All"))
}
df_prop <- df_prop %>%
  mutate(
    region_three_group = factor(paste(mimer_mapped, three_group, sep = "_"),
                                levels = order_sub)
  )

# =====  =====
three_group_colors <- c("Mimer+_Cell+" = "#E41A1C",
                        "Mimer+_Cell-" = "#377EB8",
                        "Mimer-_All"   = "#4DAF4A")
color_vector <- c()
for (region in unique_regions) {
  color_vector[paste0(region, "_Mimer+_Cell+")] <- three_group_colors["Mimer+_Cell+"]
  color_vector[paste0(region, "_Mimer+_Cell-")] <- three_group_colors["Mimer+_Cell-"]
  color_vector[paste0(region, "_Mimer-_All")]   <- three_group_colors["Mimer-_All"]
}

p_box <- ggplot(df_prop, aes(x = region_three_group, y = prop_high_risk, fill = region_three_group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7, linewidth = 0.4) +
  geom_jitter(size = 2, alpha = 0.6, shape = 21, stroke = 0.3,
              position = position_jitter(width = 0.15, height = 0.01)) +
  scale_fill_manual(values = color_vector) +
  labs(
    x = NULL,
    y = "Proportion of High Risk Cells",
    title = "High Risk Proportion by Mapped Cell Type and Three-Group Status"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
  )

signif_comparisons <- list()
signif_labels <- character()
signif_y_positions <- numeric()

y_range <- range(df_prop$prop_high_risk, na.rm = TRUE)
y_max <- y_range[2]
y_spread <- diff(y_range)
current_y_offset <- y_max + 0.03 * y_spread

for (region in unique_regions) {
  sub_data <- df_prop %>% filter(mimer_mapped == region)
  
  grp_a <- paste0(region, "_Mimer+_Cell+")
  grp_b <- paste0(region, "_Mimer+_Cell-")
  data_ab <- sub_data %>% filter(three_group %in% c("Mimer+_Cell+", "Mimer+_Cell-"))
  if (n_distinct(data_ab$three_group) == 2 && nrow(data_ab) >= 4) {
    wt <- wilcox.test(prop_high_risk ~ three_group, data = data_ab, exact = FALSE)
    pval <- wt$p.value
    label <- ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", "ns")))
    signif_comparisons[[length(signif_comparisons) + 1]] <- c(grp_a, grp_b)
    signif_labels <- c(signif_labels, label)
    signif_y_positions <- c(signif_y_positions, current_y_offset)
    current_y_offset <- current_y_offset + 0.04 * y_spread
  }
  
  grp_c <- paste0(region, "_Mimer-_All")
  data_ac <- sub_data %>% filter(three_group %in% c("Mimer+_Cell+", "Mimer-_All"))
  if (n_distinct(data_ac$three_group) == 2 && nrow(data_ac) >= 4) {
    wt <- wilcox.test(prop_high_risk ~ three_group, data = data_ac, exact = FALSE)
    pval <- wt$p.value
    label <- ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", "ns")))
    signif_comparisons[[length(signif_comparisons) + 1]] <- c(grp_a, grp_c)
    signif_labels <- c(signif_labels, label)
    signif_y_positions <- c(signif_y_positions, current_y_offset)
    current_y_offset <- current_y_offset + 0.04 * y_spread
  }
}

if (length(signif_comparisons) > 0) {
  p_box <- p_box +
    geom_signif(
      comparisons = signif_comparisons,
      annotations = signif_labels,
      y_position = signif_y_positions,
      tip_length = 0.015,
      textsize = 3.5,
      vjust = 0.5,
      color = "black"
    )
}

# ===== =====
options(repr.plot.width = 12, repr.plot.height = 6)
print(p_box)



#### S10cd
# ============================================================
# Full R pipeline for MIMER-ecotype integration, analysis, Excel export,
# and figure generation
#
# Inputs (expected in the working directory):
#   1) kmeans_full_pipeline_results.xlsx
#   2) Figure 4k l_MIMER clinical data .xlsx
#
# Main outputs:
#   1) mimer_ecotype_integrated_table.xlsx
#   2) mimer_ecotype_analysis_results_from_R.xlsx
#   3) mimer_ecotype_summary_figure.png
#   4) mimer_ecotype_summary_figure.pdf
#   5) figure_panel_a_cn_composition.png
#   6) figure_panel_b_mimer_by_ecotype.png
#   7) figure_panel_c_km_ecotype.png
#   8) figure_panel_d_sequential_cox.png
# ============================================================

# ---------------------------
# 0) Package setup
# ---------------------------
required_pkgs <- c(
  "readxl", "dplyr", "tidyr", "survival", "writexl",
  "ggplot2", "patchwork", "scales"
)

install_if_missing <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
  }
}

install_if_missing(required_pkgs)

library(readxl)
library(dplyr)
library(tidyr)
library(survival)
library(writexl)
library(ggplot2)
library(patchwork)
library(scales)

# ---------------------------
# 1) File paths
# ---------------------------
kmeans_file <- "~/Downloads/mimer_ecotype/kmeans_full_pipeline_results.xlsx"
mimer_file  <- "~/Downloads/mimer_ecotype/Figure 4k l_MIMER clinical data .xlsx"

integrated_output <- "~/Downloads/mimer_ecotype/mimer_ecotype_output2/mimer_ecotype_integrated_table.xlsx"
results_output    <- "~/Downloads/mimer_ecotype/mimer_ecotype_output2/mimer_ecotype_analysis_results_from_R.xlsx"

panel_a_file <- "~/Downloads/mimer_ecotype/mimer_ecotype_output2/figure_panel_a_cn_composition.png"
panel_b_file <- "~/Downloads/mimer_ecotype/mimer_ecotype_output2/figure_panel_b_mimer_by_ecotype.png"
panel_c_file <- "~/Downloads/mimer_ecotype/mimer_ecotype_output2/figure_panel_c_km_ecotype.png"
panel_d_file <- "~/Downloads/mimer_ecotype/mimer_ecotype_output2/figure_panel_d_sequential_cox.png"

summary_png  <- "~/Downloads/mimer_ecotype/mimer_ecotype_output2/mimer_ecotype_summary_figure.png"
summary_pdf  <- "~/Downloads/mimer_ecotype/mimer_ecotype_output2/mimer_ecotype_summary_figure.pdf"

# ---------------------------
# 2) Helper functions
# ---------------------------
safe_numeric <- function(x) suppressWarnings(as.numeric(x))

glm_or_table <- function(model) {
  sm <- summary(model)$coefficients
  ci <- suppressMessages(confint.default(model))
  data.frame(
    Variable   = rownames(sm),
    Coef_logit = sm[, "Estimate"],
    SE         = sm[, "Std. Error"],
    OR         = exp(sm[, "Estimate"]),
    CI95_low   = exp(ci[, 1]),
    CI95_high  = exp(ci[, 2]),
    P_value    = sm[, "Pr(>|z|)"],
    row.names = NULL,
    check.names = FALSE
  )
}

cox_hr_table <- function(model) {
  sm <- summary(model)
  ci <- sm$conf.int
  data.frame(
    Variable   = rownames(sm$coefficients),
    Coef_logHR = sm$coefficients[, "coef"],
    SE         = sm$coefficients[, "se(coef)"],
    HR         = ci[, "exp(coef)"],
    CI95_low   = ci[, "lower .95"],
    CI95_high  = ci[, "upper .95"],
    P_value    = sm$coefficients[, "Pr(>|z|)"],
    row.names = NULL,
    check.names = FALSE
  )
}

km_median <- function(time, event) {
  fit <- survfit(Surv(time, event) ~ 1)
  tb <- summary(fit)$table
  if ("median" %in% names(tb)) {
    return(unname(tb["median"]))
  }
  return(NA_real_)
}

extract_survfit_df <- function(fit) {
  s <- summary(fit)
  out <- data.frame(
    time = s$time,
    surv = s$surv,
    strata = s$strata,
    n_risk = s$n.risk,
    n_event = s$n.event,
    n_censor = s$n.censor,
    row.names = NULL
  )
  out$strata <- gsub("^Subtype=", "", out$strata)
  out
}

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, formatC(p, format = "e", digits = 2), sprintf("%.3f", p)))
}

# ---------------------------
# 3) Read input workbooks
# ---------------------------
sheet_names_kmeans <- excel_sheets(kmeans_file)
sheet_names_mimer  <- excel_sheets(mimer_file)

subtype_assignments <- read_xlsx(kmeans_file, sheet = "Subtype_assignments")
feature_summary     <- read_xlsx(kmeans_file, sheet = "Subtype_feature_summary")
benchmark_df        <- read_xlsx(kmeans_file, sheet = "Benchmark_k_selection")
mimer_df            <- read_xlsx(mimer_file, sheet = sheet_names_mimer[1])

eco_order   <- c("Ecotype_1", "Ecotype_2", "Ecotype_3")
final_order <- c("S1_-dominant", "S2_-dominant", "S3_-dominant")
final_to_eco <- c("S1_-dominant" = "Ecotype_1",
                  "S2_-dominant" = "Ecotype_2",
                  "S3_-dominant" = "Ecotype_3")

# ---------------------------
# 4) Create integrated table
# ---------------------------
cn_cols <- grep("^CN[0-9]+$", names(subtype_assignments), value = TRUE)

integrated_df <- mimer_df %>%
  left_join(
    subtype_assignments %>%
      select(
        Patient_ID, Sample, Cluster_raw, Final_Tumor_subtype, Subtype, all_of(cn_cols)
      ),
    by = c("Patient_ID", "Sample")
  ) %>%
  mutate(
    MIMER_high       = ifelse(tolower(MIMER) == "high", 1L, 0L),
    Age_2_num        = safe_numeric(`Age_2 (1, >=median; 0, <median)`),
    Grade_2_num      = safe_numeric(Grade_2),
    AJCC_Stage_2_num = safe_numeric(AJCC_Stage_2),
    time_month       = safe_numeric(`Survival_time (month)`),
    event            = safe_numeric(`Survival_status (1,death; 0, alive)`),
    Subtype          = factor(Subtype, levels = eco_order)
  )

subtype_assignments <- subtype_assignments[subtype_assignments$Sample %in% unique(integrated_df$Sample),]
# Export integrated table
write_xlsx(list(Integrated_Data = integrated_df), path = integrated_output)

# ---------------------------
# 5) Analysis 1: MIMER distribution across ecotypes
# ---------------------------
mimer_table <- table(integrated_df$Subtype, integrated_df$MIMER)
chi_res <- chisq.test(mimer_table)

mimer_by_ecotype <- as.data.frame.matrix(mimer_table)
mimer_by_ecotype$Subtype <- rownames(mimer_by_ecotype)
mimer_by_ecotype <- mimer_by_ecotype %>%
  relocate(Subtype) %>%
  mutate(
    Total = high + low,
    Pct_high = high / Total,
    Pct_high_display = sprintf("%.1f%%", 100 * Pct_high)
  )

chi_expected <- as.data.frame.matrix(chi_res$expected)
chi_expected$Subtype <- rownames(chi_expected)
chi_expected <- chi_expected %>% relocate(Subtype)

# ---------------------------
# 6) Analysis 2: Logistic regression
#    MIMER_high ~ Subtype + covariates
# ---------------------------
logistic_df <- integrated_df %>%
  filter(
    complete.cases(MIMER_high, Subtype, Age_2_num, Grade_2_num, AJCC_Stage_2_num)
  )

logit_fit <- glm(
  MIMER_high ~ Subtype + Age_2_num + Grade_2_num + AJCC_Stage_2_num,
  data = logistic_df,
  family = binomial(),
  control = glm.control(maxit = 100)
)

logit_results <- glm_or_table(logit_fit)

# ---------------------------
# 7) Analysis 3: Survival models
# ---------------------------
# 7a) Original ecotype survival from the subtype assignment table
surv_original_df <- subtype_assignments %>%
  mutate(
    time_month = safe_numeric(`Survival_time (months)`),
    event = safe_numeric(`Survival_status (1, dead; 0, alive)`),
    Subtype = factor(Subtype, levels = eco_order)
  ) %>%
  filter(complete.cases(time_month, event, Subtype))

km_fit_original <- survfit(Surv(time_month, event) ~ Subtype, data = surv_original_df)
km_df_original <- extract_survfit_df(km_fit_original)

logrank_original <- survdiff(Surv(time_month, event) ~ Subtype, data = surv_original_df)
logrank_p_original <- pchisq(
  q = logrank_original$chisq,
  df = length(logrank_original$n) - 1,
  lower.tail = FALSE
)

# 7b) Merged-cohort Cox models
cox_df <- integrated_df %>%
  filter(
    complete.cases(time_month, event, Subtype, MIMER_high, Age_2_num, Grade_2_num, AJCC_Stage_2_num)
  )

cox_fit_1 <- coxph(
  Surv(time_month, event) ~ Subtype + Age_2_num + Grade_2_num + AJCC_Stage_2_num,
  data = cox_df
)

cox_fit_2 <- coxph(
  Surv(time_month, event) ~ MIMER_high + Age_2_num + Grade_2_num + AJCC_Stage_2_num,
  data = cox_df
)

cox_fit_3 <- coxph(
  Surv(time_month, event) ~ Subtype + MIMER_high + Age_2_num + Grade_2_num + AJCC_Stage_2_num,
  data = cox_df
)

cox_fit_4 <- coxph(
  Surv(time_month, event) ~ Subtype * MIMER_high + Age_2_num + Grade_2_num + AJCC_Stage_2_num,
  data = cox_df
)

cox_results_1 <- cox_hr_table(cox_fit_1)
cox_results_2 <- cox_hr_table(cox_fit_2)
cox_results_3 <- cox_hr_table(cox_fit_3)
cox_results_4 <- cox_hr_table(cox_fit_4)

model_comparison <- data.frame(
  Model = c(
    "Model_1_Ecotype_plus_covariates",
    "Model_2_MIMER_plus_covariates",
    "Model_3_Ecotype_plus_MIMER_plus_covariates",
    "Model_4_Interaction_model"
  ),
  N = nrow(cox_df),
  Events = sum(cox_df$event),
  Parameters = c(length(coef(cox_fit_1)), length(coef(cox_fit_2)),
                 length(coef(cox_fit_3)), length(coef(cox_fit_4))),
  LogLik = c(as.numeric(logLik(cox_fit_1)), as.numeric(logLik(cox_fit_2)),
             as.numeric(logLik(cox_fit_3)), as.numeric(logLik(cox_fit_4))),
  AIC = c(AIC(cox_fit_1), AIC(cox_fit_2), AIC(cox_fit_3), AIC(cox_fit_4)),
  check.names = FALSE
)

# Likelihood-ratio tests computed manually
lrt_13_chisq <- 2 * (as.numeric(logLik(cox_fit_3)) - as.numeric(logLik(cox_fit_1)))
lrt_13_df    <- length(coef(cox_fit_3)) - length(coef(cox_fit_1))
lrt_13_p     <- pchisq(lrt_13_chisq, df = lrt_13_df, lower.tail = FALSE)

lrt_34_chisq <- 2 * (as.numeric(logLik(cox_fit_4)) - as.numeric(logLik(cox_fit_3)))
lrt_34_df    <- length(coef(cox_fit_4)) - length(coef(cox_fit_3))
lrt_34_p     <- pchisq(lrt_34_chisq, df = lrt_34_df, lower.tail = FALSE)

lrt_tests <- data.frame(
  Comparison = c(
    "Model_3 vs Model_1 (add MIMER_high)",
    "Model_4 vs Model_3 (add interaction)"
  ),
  Chi_square = c(lrt_13_chisq, lrt_34_chisq),
  df = c(lrt_13_df, lrt_34_df),
  P_value = c(lrt_13_p, lrt_34_p),
  check.names = FALSE
)

# Effect attenuation
logHR_model1 <- cox_results_1$Coef_logHR[cox_results_1$Variable == "SubtypeEcotype_2"]
logHR_model3 <- cox_results_3$Coef_logHR[cox_results_3$Variable == "SubtypeEcotype_2"]

effect_attenuation <- data.frame(
  Contrast = "Ecotype_2 vs Ecotype_1",
  Model1_logHR = logHR_model1,
  Model3_logHR = logHR_model3,
  Attenuation_fraction = (abs(logHR_model1) - abs(logHR_model3)) / abs(logHR_model1),
  check.names = FALSE
) %>%
  mutate(
    Attenuation_percent = sprintf("%.1f%%", 100 * Attenuation_fraction)
  )

# Stratified descriptive summary
stratified_summary <- do.call(
  rbind,
  lapply(eco_order, function(eco) {
    do.call(
      rbind,
      lapply(c("high", "low"), function(mim_flag) {
        g <- integrated_df %>%
          filter(Subtype == eco, MIMER == mim_flag)
        if (nrow(g) == 0) return(NULL)
        data.frame(
          Subtype = eco,
          MIMER = mim_flag,
          N = nrow(g),
          Events = sum(g$event, na.rm = TRUE),
          Median_survival_months_KM = km_median(g$time_month, g$event),
          row.names = NULL,
          check.names = FALSE
        )
      })
    )
  })
)

# ---------------------------
# 8) Overview sheet
# ---------------------------
benchmark_k3 <- benchmark_df %>% filter(k == 3) %>% slice(1)

overview <- data.frame(
  Metric = c(
    "Input_kmeans_workbook",
    "Input_mimer_workbook",
    "Merged_samples",
    "Complete_cases_logistic",
    "Complete_cases_cox",
    "Total_events_cox",
    "Original_ecotype_logrank_p",
    "Merged_chi_square_MIMER_by_ecotype",
    "Merged_chi_square_p_value",
    "Strongest_logistic_signal",
    "MIMER_survival_signal",
    "Ecotype_effect_after_MIMER",
    "Add_MIMER_LRT_p",
    "Interaction_LRT_p"
  ),
  Value = c(
    kmeans_file,
    mimer_file,
    nrow(integrated_df),
    nrow(logistic_df),
    nrow(cox_df),
    sum(cox_df$event),
    signif(logrank_p_original, 4),
    signif(unname(chi_res$statistic), 4),
    signif(chi_res$p.value, 4),
    "Ecotype_2 vs Ecotype_1 OR expected ~0.03 if reproduced",
    "MIMER_high HR expected >1 in adjusted model",
    "Ecotype_2 HR attenuates after adding MIMER_high",
    signif(lrt_13_p, 4),
    signif(lrt_34_p, 4)
  ),
  check.names = FALSE
)

readme_sheet <- data.frame(
  Item = c(
    "Analysis scope",
    "Integrated table",
    "Cohort used for panel c",
    "Cohort used for logistic/Cox",
    "Covariates used",
    "Main interpretation"
  ),
  Value = c(
    "Integrates the uploaded ecotype workbook with the uploaded MIMER clinical workbook by Patient_ID + Sample.",
    integrated_output,
    "Original ecotype cohort from Subtype_assignments in kmeans_full_pipeline_results.xlsx.",
    "Merged cohort after joining the two uploaded workbooks.",
    "Age_2_num + Grade_2_num + AJCC_Stage_2_num.",
    "Differential MIMER enrichment may partially explain ecotype-associated survival heterogeneity."
  ),
  check.names = FALSE
)

# ---------------------------
# 9) Export analysis workbook
# ---------------------------
results_list <- list(
  README = readme_sheet,
  Integrated_Data = integrated_df,
  MIMER_by_Ecotype = mimer_by_ecotype,
  ChiSquare_Expected = chi_expected,
  Logistic_MIMER_high = logit_results,
  Cox_Model_1 = cox_results_1,
  Cox_Model_2 = cox_results_2,
  Cox_Model_3 = cox_results_3,
  Cox_Model_4_Interaction = cox_results_4,
  Model_Comparison = model_comparison,
  LRT_Tests = lrt_tests,
  Effect_Attenuation = effect_attenuation,
  Stratified_Summary = stratified_summary,
  KM_Original_Ecotypes = km_df_original,
  Overview = overview
)

write_xlsx(results_list, path = results_output)

# ---------------------------
# 10) Plotting
# ---------------------------

# Panel A: CN composition heatmap by ecotype
feature_plot_df <- feature_summary %>%
  mutate(Subtype = final_to_eco[Final_Tumor_subtype]) %>%
  filter(!is.na(Subtype)) %>%
  select(Subtype, all_of(cn_cols))

top_cns <- feature_plot_df %>%
  summarise(across(all_of(cn_cols), ~ var(.x, na.rm = TRUE))) %>%
  pivot_longer(cols = everything(), names_to = "CN", values_to = "Var") %>%
  arrange(desc(Var)) %>%
  slice(1:10) %>%
  pull(CN)

heatmap_df <- feature_plot_df %>%
  select(Subtype, all_of(top_cns)) %>%
  pivot_longer(cols = -Subtype, names_to = "CN", values_to = "MeanProportion") %>%
  mutate(
    Subtype = factor(Subtype, levels = eco_order),
    CN = factor(CN, levels = top_cns)
  )

panel_a <- ggplot(heatmap_df, aes(x = CN, y = Subtype, fill = MeanProportion)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", MeanProportion)), size = 3.0) +
  labs(
    title = "a  CN composition distinguishes tumor ecotypes",
    x = NULL, y = NULL, fill = "Mean proportion"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

# Panel B: MIMER-high proportion by ecotype
bar_df <- integrated_df %>%
  group_by(Subtype) %>%
  summarise(
    High = sum(MIMER_high, na.rm = TRUE),
    Total = n(),
    Proportion = High / Total,
    .groups = "drop"
  ) %>%
  mutate(Subtype = factor(Subtype, levels = eco_order))

panel_b <- ggplot(bar_df, aes(x = Subtype, y = Proportion)) +
  geom_col(width = 0.65) +
  geom_text(
    aes(
      label = paste0(High, "/", Total, "\n", sprintf("%.1f%%", 100 * Proportion)),
      y = pmin(Proportion + 0.06, 0.96)
    ),
    size = 3.5
  ) +
  annotate(
    "text", x = 3, y = 0.05,
    label = paste0("Chi-square p = ", formatC(chi_res$p.value, format = "e", digits = 2)),
    hjust = 1, size = 3.5
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.0)) +
  labs(
    title = "b  MIMER-high is unevenly distributed across ecotypes",
    x = NULL, y = "MIMER-high proportion"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# Panel C: Kaplan-Meier by ecotype (original ecotype cohort)
panel_c <- ggplot(km_df_original, aes(x = time, y = surv, color = strata)) +
  geom_step(linewidth = 0.9) +
  annotate(
    "text", x = max(km_df_original$time, na.rm = TRUE) * 0.98, y = 0.05,
    label = paste0("Log-rank p = ", formatC(logrank_p_original, format = "f", digits = 4)),
    hjust = 1, size = 3.5
  ) +
  coord_cartesian(ylim = c(0, 1.02)) +
  labs(
    title = "c  Tumor ecotypes show distinct survival outcomes",
    x = "Time (months)",
    y = "Overall survival probability",
    color = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

panel_c_u <- ggsurvplot(
  km_fit_original,
  data = surv_original_df,
  pval = TRUE,
  conf.int = TRUE, 
  fun = "pct",
  risk.table = TRUE,
  risk.table.height = 0.25,
  linetype = "strata",
  legend = "bottom",
  legend.labs = c("Ecotype_1", "Ecotype_2", "Ecotype_3"),
  legend.title = "Subtype",
  title = "c  Tumor ecotypes show distinct survival outcomes",
  xlab = "Time (months)",
  ylab = "Overall survival probability",
)
pdf("~/Downloads/mimer_ecotype/mimer_ecotype_output2/panel_c_u.pdf", width = 8, height = 6)
print(panel_c_u)
dev.off()
# Panel D: Sequential Cox forest plot
forest_df <- bind_rows(
  cox_results_1 %>%
    filter(Variable %in% c("SubtypeEcotype_2", "SubtypeEcotype_3")) %>%
    mutate(Label = c("Model 1: Ecotype_2", "Model 1: Ecotype_3")),
  cox_results_2 %>%
    filter(Variable == "MIMER_high") %>%
    mutate(Label = "Model 2: MIMER-high"),
  cox_results_3 %>%
    filter(Variable %in% c("SubtypeEcotype_2", "SubtypeEcotype_3", "MIMER_high")) %>%
    mutate(Label = c("Model 3: Ecotype_2", "Model 3: Ecotype_3", "Model 3: MIMER-high"))
) %>%
  mutate(
    Label = factor(Label, levels = rev(c(
      "Model 1: Ecotype_2", "Model 1: Ecotype_3",
      "Model 2: MIMER-high",
      "Model 3: Ecotype_2", "Model 3: Ecotype_3", "Model 3: MIMER-high"
    )))
  )

panel_d <- ggplot(forest_df, aes(x = HR, y = Label)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_errorbarh(aes(xmin = CI95_low, xmax = CI95_high), width = 0.18) +
  geom_point(size = 2.2) +
  geom_text(
    aes(
      label = paste0("HR ", sprintf("%.2f", HR), ", p=", fmt_p(P_value))
    ),
    position = "dodge",
    vjust = -1.2, size = 3
  ) +
  scale_x_log10() +
  labs(
    title = "d  MIMER attenuates the ecotype-associated survival effect",
    x = "Hazard ratio (log scale)",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# Save single panels
ggsave(panel_a_file, panel_a, width = 8.6, height = 3.8, dpi = 300)
ggsave(panel_b_file, panel_b, width = 6.6, height = 4.2, dpi = 300)
ggsave(panel_c_file, panel_c, width = 6.8, height = 4.4, dpi = 300)
ggsave(panel_d_file, panel_d, width = 8, height = 4.8, dpi = 300)

# Combined figure
combined_figure <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_annotation(
    title = "Differential MIMER enrichment partially explains ecotype-associated survival heterogeneity"
  )

ggsave(summary_png, combined_figure, width = 16, height = 10, dpi = 300)
ggsave(summary_pdf, combined_figure, width = 16, height = 10)

# ---------------------------
# 11) Console summary
# ---------------------------
cat("Done.\n")
cat("Integrated table: ", integrated_output, "\n", sep = "")
cat("Analysis workbook: ", results_output, "\n", sep = "")
cat("Summary figure PNG: ", summary_png, "\n", sep = "")
cat("Summary figure PDF: ", summary_pdf, "\n", sep = "")
cat("Key merged-cohort chi-square p: ", formatC(chi_res$p.value, format = "e", digits = 3), "\n", sep = "")
cat("Original ecotype log-rank p: ", formatC(logrank_p_original, format = "f", digits = 4), "\n", sep = "")

