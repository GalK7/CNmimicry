# Purpose:
#   Compare ALD/DLD burden and gene-level ALD/DLD enrichment between
#   IDH-mutant and IDH-wildtype LGG samples.
#
# Inputs:
#   - Core workflow ALD/DLD call files from result directories defined in config.R.
#   - Cached or downloaded TCGA LGG somatic mutation data for IDH status calling.
#
# Outputs:
#   - IDH-stratified burden and gene-enrichment CSV files written to analysis/outputs/idh.

source("config.R")

# ---- settings ----
cohort  <- "LGG"
prev_thr <- prevalence_thresholds[1]
prev_tag <- format(prev_thr, nsmall = 2)

ald_tpm_filter_enabled <- !is.null(ald_min_amplified_median_tpm)
if (ald_tpm_filter_enabled &&
    (!is.numeric(ald_min_amplified_median_tpm) ||
     length(ald_min_amplified_median_tpm) != 1L ||
     !is.finite(ald_min_amplified_median_tpm) ||
     ald_min_amplified_median_tpm < 0)) {
  stop("ald_min_amplified_median_tpm must be NULL or one finite nonnegative number.")
}
ald_tpm_tag <- if (ald_tpm_filter_enabled) {
  paste0("aldTPMgt", format(ald_min_amplified_median_tpm, trim = TRUE, scientific = FALSE))
} else {
  "aldTPMoff"
}
ald_tpm_suffix <- if (ald_tpm_filter_enabled) paste0("_", ald_tpm_tag) else ""

sig_tag <- if (isTRUE(use_gistic_significant)) {
  paste0("with_SIG_", gistic_sig_gene_source)
} else {
  "no_SIG"
}
gthr_tag <- if (isTRUE(use_gistic_thresholded_intersection)) "with_GISTICthr" else "no_GISTICthr"
cnv_tag <- if (isTRUE(use_cnv_median_normalization)) "cnvnorm" else "rawcnv"

burden_type <- "ALD"
burden_label <- if (burden_type == "ALD") "ALD" else if (burden_type == "DLD") "DLD" else NA_character_
force_rebuild_maf_cache <- FALSE

get_maf_cached_ensemble <- function(gdc_project, force_rebuild = FALSE) {
  maf_dir <- file.path(workflow_root, "analysis", "cache", "maf_ensemble_masked")
  maf_gdc_dir <- file.path(workflow_root, "analysis", "GDCdata")
  dir.create(maf_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(maf_gdc_dir, recursive = TRUE, showWarnings = FALSE)
  
  maf_path <- file.path(maf_dir, paste0(gdc_project, "_ensemble_masked_maf.rds"))
  if (file.exists(maf_path) && !isTRUE(force_rebuild)) {
    message("[IDH_LGG] Using cached MAF: ", maf_path)
    return(readRDS(maf_path))
  }
  
  message("[IDH_LGG] Downloading MAF for ", gdc_project)
  query_maf <- TCGAbiolinks::GDCquery(
    project       = gdc_project,
    data.category = "Simple Nucleotide Variation",
    data.type     = "Masked Somatic Mutation",
    access        = "open",
    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
  )
  
  TCGAbiolinks::GDCdownload(query_maf, directory = maf_gdc_dir)
  maf_data <- TCGAbiolinks::GDCprepare(query_maf, directory = maf_gdc_dir)
  
  saveRDS(maf_data, maf_path)
  message("[IDH_LGG] Saved MAF cache: ", maf_path)
  maf_data
}

if (isTRUE(use_gistic_thresholded_intersection)) {
  results_dir_base <- results_dir
  results_dir_sig_base <- results_dir_sig
  results_dir_cnvnorm_base <- results_dir_cnvnorm
} else {
  results_dir_base <- results_dir_no_gistic_thr
  results_dir_sig_base <- results_dir_sig_no_gistic_thr
  results_dir_cnvnorm_base <- results_dir_cnvnorm_no_gistic_thr
}

results_base <- if (isTRUE(use_cnv_median_normalization)) {
  file.path(results_dir_cnvnorm_base, sig_tag)
} else {
  if (isTRUE(use_gistic_significant)) results_dir_sig_base else results_dir_base
}

if (ald_tpm_filter_enabled) results_base <- file.path(results_base, ald_tpm_tag)
res_dir <- file.path(results_base, cohort, paste0("prev_", format(prev_thr, nsmall = 2)))


# ---- load workflow outputs ----
if (burden_type == "ALD") {
  burden_file <- "ald_per_sample.rds"
  frac_col <- "ald_frac"
  n_col <- "n_ald_genes"
} else if (burden_type == "DLD") {
  burden_file <- "dld_per_sample.rds"
  frac_col <- "dld_frac"
  n_col <- "n_dld_genes"
} else {
  stop("burden_type must be 'ALD' or 'DLD'")
}

burden_path <- file.path(res_dir, burden_file)
if (!file.exists(burden_path)) stop("Missing: ", burden_path)

df <- readRDS(burden_path)
df$patient <- substr(df$Sample, 1, 12)
# ---- query MAF + call IDH status (hotspots) ----
if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) stop("Install TCGAbiolinks")
if (!requireNamespace("dplyr", quietly = TRUE)) stop("Install dplyr")

library(TCGAbiolinks)
library(dplyr)

maf_data <- get_maf_cached_ensemble(
  gdc_project = "TCGA-LGG",
  force_rebuild = force_rebuild_maf_cache
)

maf_data$patient <- substr(maf_data$Tumor_Sample_Barcode, 1, 12)

idh_hotspots <- maf_data |>
  filter(
    (Hugo_Symbol == "IDH1" & startsWith(HGVSp_Short, "p.R132")) |
      (Hugo_Symbol == "IDH2" & startsWith(HGVSp_Short, "p.R172"))
  ) |>
  distinct(patient)

idh_mut_patients <- idh_hotspots$patient


# ------------------------------------------------------------------------------
# ---- PDGFRA ALD enrichment test ----

ald_calls_file <- file.path(res_dir, "ald_calls_long.rds")
if (!file.exists(ald_calls_file)) stop("Missing: ", ald_calls_file)

ald_calls_long <- readRDS(ald_calls_file)

pdgfra_ald <- ald_calls_long %>%
  filter(Gene == "PDGFRA") %>%
  mutate(patient = substr(Sample, 1, 12)) %>%
  distinct(patient)

# all LGG patients in this analysis
all_patients <- unique(df$patient)

# binary PDGFRA-ALD indicator
pdgfra_flag <- all_patients %in% pdgfra_ald$patient

# IDH status aligned to same patient universe
idh_status_all <- ifelse(all_patients %in% idh_mut_patients, "IDH-mut", "IDH-wt")

# 2x2 contingency table
cont_table <- table(
  IDH_status = idh_status_all,
  PDGFRA_ALD = pdgfra_flag
)

print("Contingency table (IDH x PDGFRA-ALD):")
print(cont_table)

# Fisher exact test
fisher_res <- fisher.test(cont_table)

print("Fisher exact test:")
print(fisher_res)

# ------------------------------------------------------------------------------

df$IDH_status <- ifelse(df$patient %in% idh_mut_patients, "IDH-mut", "IDH-wt")
df$IDH_status <- factor(df$IDH_status, levels = c("IDH-wt", "IDH-mut"))


# ---- save per-sample data for plotting ----
plot_df <- df %>%
  transmute(
    cohort = cohort,
    prev_thr = prev_thr,
    model = sig_tag,
    Sample,
    patient,
    IDH_status,
    burden_frac = .data[[frac_col]],
    burden_n = .data[[n_col]],
    burden_type = burden_label
  )


# ---- stats ----
n_per_group <- table(df$IDH_status)
med <- tapply(df[[frac_col]], df$IDH_status, median, na.rm = TRUE)
iqr <- tapply(df[[frac_col]], df$IDH_status, IQR, na.rm = TRUE)

ratio_mut_wt <- unname(med["IDH-mut"] / med["IDH-wt"])
x <- df[[frac_col]][df$IDH_status == "IDH-wt"]
y <- df[[frac_col]][df$IDH_status == "IDH-mut"]
wilcox_p <- wilcox.test(x, y, exact = FALSE)$p.value

print(n_per_group)
print(med)
print(iqr)
print(ratio_mut_wt)
print(wilcox_p)

# ---- save output ----
out <- data.frame(
  cohort = cohort,
  prev_thr = prev_thr,
  n_wt  = unname(n_per_group["IDH-wt"]),
  n_mut = unname(n_per_group["IDH-mut"]),
  med_wt  = unname(med["IDH-wt"]),
  med_mut = unname(med["IDH-mut"]),
  ratio_mut_wt = ratio_mut_wt,
  wilcox_p = wilcox_p
)

output_dir <- file.path(workflow_root, "analysis", "outputs", "idh")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

plot_file <- file.path(
  output_dir,
  paste0("LGG_IDH_", burden_label, "_per_sample_", gthr_tag, "_", sig_tag, "_", cnv_tag, "_prev_", prev_tag, ald_tpm_suffix, ".csv")
)
write.csv(plot_df, plot_file, row.names = FALSE)
message("Saved per-sample plot data: ", plot_file)


out_file <- file.path(
  output_dir,
  paste0("IDH_LGG_", burden_label, "_burden_test_", gthr_tag, "_", sig_tag, "_", cnv_tag, "_prev_", prev_tag, ald_tpm_suffix, ".csv")
)
write.csv(out, out_file, row.names = FALSE)
message("Saved: ", out_file)



# ------------------------------------------------------------------
# ---- Gene-level ALD / DLD enrichment in IDH-mut LGG ----
all_patients <- unique(df$patient)
idh_status_df <- data.frame(
  patient = all_patients,
  IDH_status = ifelse(all_patients %in% idh_mut_patients, "IDH-mut", "IDH-wt")
)

test_gene_enrichment <- function(calls_long, feature_label) {
  
  genes <- sort(unique(calls_long$Gene))
  
  res <- lapply(genes, function(g) {
    gene_patients <- unique(substr(calls_long$Sample[calls_long$Gene == g], 1, 12))
    
    tmp <- idh_status_df
    tmp$feature_positive <- tmp$patient %in% gene_patients
    
    cont_table <- table(tmp$IDH_status, tmp$feature_positive)
    fisher_res <- fisher.test(cont_table)
    
    n_wt  <- sum(tmp$IDH_status == "IDH-wt")
    n_mut <- sum(tmp$IDH_status == "IDH-mut")
    
    n_feature_wt  <- sum(tmp$feature_positive[tmp$IDH_status == "IDH-wt"])
    n_feature_mut <- sum(tmp$feature_positive[tmp$IDH_status == "IDH-mut"])
    
    data.frame(
      gene = g,
      feature = feature_label,
      n_IDH_wt = n_wt,
      n_IDH_mut = n_mut,
      n_positive_IDH_wt = n_feature_wt,
      n_positive_IDH_mut = n_feature_mut,
      pct_positive_IDH_wt = 100 * n_feature_wt / n_wt,
      pct_positive_IDH_mut = 100 * n_feature_mut / n_mut,
      p_value = fisher_res$p.value
    )
  })
  
  res <- bind_rows(res) %>%
    mutate(
      FDR = p.adjust(p_value, method = "BH"),
      direction = ifelse(pct_positive_IDH_mut > pct_positive_IDH_wt,
                         "enriched_in_IDH_mut",
                         "depleted_in_IDH_mut")
    ) %>%
    arrange(FDR, p_value)
  
  res
}

# ---- ALD gene-level enrichment ----
ald_calls_file <- file.path(res_dir, "ald_calls_long.rds")
if (!file.exists(ald_calls_file)) stop("Missing: ", ald_calls_file)

ald_calls_long <- readRDS(ald_calls_file)

ald_gene_results <- test_gene_enrichment(ald_calls_long, "ALD") %>%
  filter(direction == "enriched_in_IDH_mut", FDR < 0.05) %>%
  select(-direction)

ald_gene_file <- file.path(
  output_dir,
  paste0("LGG_IDH_ALD_gene_enrichment_", gthr_tag, "_", sig_tag, "_", cnv_tag, "_prev_", prev_tag, ald_tpm_suffix, ".csv")
)
write.csv(ald_gene_results, ald_gene_file, row.names = FALSE)
message("Saved: ", ald_gene_file)

# ---- DLD gene-level enrichment ----
dld_calls_file <- file.path(res_dir, "dld_calls_long.rds")
if (!file.exists(dld_calls_file)) stop("Missing: ", dld_calls_file)

dld_calls_long <- readRDS(dld_calls_file)

dld_gene_results <- test_gene_enrichment(dld_calls_long, "DLD") %>%
  filter(direction == "enriched_in_IDH_mut", FDR < 0.05) %>%
  select(-direction)

dld_gene_file <- file.path(
  output_dir,
  paste0("LGG_IDH_DLD_gene_enrichment_", gthr_tag, "_", sig_tag, "_", cnv_tag, "_prev_", prev_tag, ald_tpm_suffix, ".csv")
)
write.csv(dld_gene_results, dld_gene_file, row.names = FALSE)
message("Saved: ", dld_gene_file)
