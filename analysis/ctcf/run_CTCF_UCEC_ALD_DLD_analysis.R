# Purpose:
#   Compare ALD/DLD burden and gene-level ALD/DLD enrichment between
#   CTCF-mutant and CTCF-wildtype UCEC samples.
#
# Inputs:
#   - Core workflow ALD/DLD call files from result directories defined in config.R.
#   - Cached or downloaded TCGA UCEC somatic mutation data for CTCF status calling.
#   - TCGA UCEC subtype annotations used for MSI and POLE-related sensitivity checks.
#
# Outputs:
#   - CTCF-stratified burden and gene-enrichment CSV files written to analysis/outputs/ctcf.

source("config.R")

if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) stop("Install TCGAbiolinks")
if (!requireNamespace("dplyr", quietly = TRUE)) stop("Install dplyr")
if (!requireNamespace("tibble", quietly = TRUE)) stop("Install tibble")

library(TCGAbiolinks)
library(dplyr)
library(tibble)

# ---- settings ----
cohort  <- "UCEC"
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
force_rebuild_maf_cache <- FALSE

sig_tag <- if (isTRUE(use_gistic_significant)) {
  paste0("with_SIG_", gistic_sig_gene_source)
} else {
  "no_SIG"
}
gthr_tag <- if (isTRUE(use_gistic_thresholded_intersection)) "with_GISTICthr" else "no_GISTICthr"
cnv_tag  <- if (isTRUE(use_cnv_median_normalization)) "cnvnorm" else "rawcnv"

get_maf_cached_ensemble <- function(gdc_project, force_rebuild = FALSE) {
  maf_dir <- file.path(workflow_root, "analysis", "cache", "maf_ensemble_masked")
  maf_gdc_dir <- file.path(workflow_root, "analysis", "GDCdata")
  dir.create(maf_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(maf_gdc_dir, recursive = TRUE, showWarnings = FALSE)
  
  maf_path <- file.path(maf_dir, paste0(gdc_project, "_ensemble_masked_maf.rds"))
  if (file.exists(maf_path) && !isTRUE(force_rebuild)) {
    message("[CTCF_UCEC] Using cached MAF: ", maf_path)
    return(readRDS(maf_path))
  }
  
  message("[CTCF_UCEC] Downloading MAF for ", gdc_project)
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
  message("[CTCF_UCEC] Saved MAF cache: ", maf_path)
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

res_base <- if (isTRUE(use_cnv_median_normalization)) {
  file.path(results_dir_cnvnorm_base, sig_tag)
} else {
  if (isTRUE(use_gistic_significant)) results_dir_sig_base else results_dir_base
}

if (ald_tpm_filter_enabled) res_base <- file.path(res_base, ald_tpm_tag)
res_dir <- file.path(res_base, cohort, paste0("prev_", format(prev_thr, nsmall = 2)))

ald_file <- file.path(res_dir, "ald_per_sample.rds")
dld_file <- file.path(res_dir, "dld_per_sample.rds")

if (!file.exists(ald_file)) stop("Missing: ", ald_file)
if (!file.exists(dld_file)) stop("Missing: ", dld_file)

ald_df <- readRDS(ald_file)
dld_df <- readRDS(dld_file)

ald_df$patient <- substr(ald_df$Sample, 1, 12)
dld_df$patient <- substr(dld_df$Sample, 1, 12)

# Align ALD sample universe to DLD so both summaries use the same tumors.
all_samples <- tibble::tibble(Sample = unique(dld_df$Sample))

ald_df <- all_samples %>%
  left_join(ald_df, by = "Sample") %>%
  mutate(
    ald_frac = dplyr::coalesce(ald_frac, 0),
    n_ald_genes = dplyr::coalesce(n_ald_genes, 0L),
    patient = substr(Sample, 1, 12)
  )

# ---- MAF: CTCF mutation status ----
maf <- get_maf_cached_ensemble(
  gdc_project = "TCGA-UCEC",
  force_rebuild = force_rebuild_maf_cache
)

maf$patient <- substr(maf$Tumor_Sample_Barcode, 1, 12)

ctcf_mut_patients <- maf |>
  filter(Hugo_Symbol == "CTCF") |>
  distinct(patient) |>
  pull(patient)

ald_df$CTCF_status <- ifelse(ald_df$patient %in% ctcf_mut_patients, "CTCF-mut", "CTCF-wt")
dld_df$CTCF_status <- ifelse(dld_df$patient %in% ctcf_mut_patients, "CTCF-mut", "CTCF-wt")

ald_df$CTCF_status <- factor(ald_df$CTCF_status, levels = c("CTCF-wt","CTCF-mut"))
dld_df$CTCF_status <- factor(dld_df$CTCF_status, levels = c("CTCF-wt","CTCF-mut"))

# ---- UCEC subtype (IntegrativeCluster) ----
subtype_cache_file <- file.path(
  workflow_root, "analysis", "cache", "subtypes",
  "UCEC_integrative_cluster.csv"
)

if (file.exists(subtype_cache_file)) {
  sub <- read.csv(subtype_cache_file, stringsAsFactors = FALSE)
} else {
  sub <- TCGAquery_subtype("UCEC") |>
    transmute(
      patient = patient,
      integrative_cluster = IntegrativeCluster
    )

  dir.create(
    dirname(subtype_cache_file),
    recursive = TRUE,
    showWarnings = FALSE
  )
  write.csv(sub, subtype_cache_file, row.names = FALSE)
}

ald_df <- ald_df |> left_join(sub, by = "patient")
dld_df <- dld_df |> left_join(sub, by = "patient")

# ---- helper: summarize + wilcox ----
summarize_compare <- function(df, value_col) {
  n <- table(df$CTCF_status)
  
  med <- tapply(df[[value_col]], df$CTCF_status, median, na.rm = TRUE)
  iqr <- tapply(df[[value_col]], df$CTCF_status, IQR,    na.rm = TRUE)
  
  ratio <- unname(med["CTCF-mut"] / med["CTCF-wt"])
  p <- wilcox.test(df[[value_col]] ~ df$CTCF_status, exact = FALSE)$p.value
  
  list(n=n, med=med, iqr=iqr, ratio=ratio, p=p)
}

# ---- ALD comparisons ----
res_ald_all <- summarize_compare(ald_df, "ald_frac")

ald_msi <- ald_df |> filter(integrative_cluster == "MSI")
res_ald_msi <- summarize_compare(ald_msi, "ald_frac")

ald_no_pole <- ald_df |> filter(integrative_cluster != "POLE")
res_ald_no_pole <- summarize_compare(ald_no_pole, "ald_frac")

# ---- DLD comparisons ----
res_dld_all <- summarize_compare(dld_df, "dld_frac")

dld_msi <- dld_df |> filter(integrative_cluster == "MSI")
res_dld_msi <- summarize_compare(dld_msi, "dld_frac")

dld_no_pole <- dld_df |> filter(integrative_cluster != "POLE")
res_dld_no_pole <- summarize_compare(dld_no_pole, "dld_frac")

# ---- save tidy summary ----
output_dir <- file.path(workflow_root, "analysis", "outputs", "ctcf")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

out <- data.frame(
  cohort = cohort,
  prev_thr = prev_thr,
  
  ALD_full_cohort_ratio = res_ald_all$ratio,
  ALD_full_cohort_p     = res_ald_all$p,
  ALD_MSI_ratio = res_ald_msi$ratio,
  ALD_MSI_p     = res_ald_msi$p,
  ALD_POLE_excluded_ratio = res_ald_no_pole$ratio,
  ALD_POLE_excluded_p     = res_ald_no_pole$p,
  
  DLD_full_cohort_ratio = res_dld_all$ratio,
  DLD_full_cohort_p     = res_dld_all$p,
  DLD_MSI_ratio = res_dld_msi$ratio,
  DLD_MSI_p     = res_dld_msi$p,
  DLD_POLE_excluded_ratio = res_dld_no_pole$ratio,
  DLD_POLE_excluded_p     = res_dld_no_pole$p,
  
  stringsAsFactors = FALSE
)

out_file <- file.path(
  output_dir,
  paste0("CTCF_UCEC_ALD_DLD_burden_test_", gthr_tag, "_", sig_tag, "_", cnv_tag, "_prev_", prev_tag, ald_tpm_suffix, ".csv")
)
write.csv(out, out_file, row.names = FALSE)
message("Saved: ", out_file)

# ---- print key results ----
print(res_ald_all); print(res_ald_msi); print(res_ald_no_pole)
print(res_dld_all); print(res_dld_msi); print(res_dld_no_pole)

# ------------------------------------------------------------------
# ---- Gene-level ALD / DLD enrichment in CTCF-mut UCEC ----

all_patients <- unique(dld_df$patient)

ctcf_status_df <- data.frame(
  patient = all_patients,
  CTCF_status = ifelse(all_patients %in% ctcf_mut_patients, "CTCF-mut", "CTCF-wt")
)

test_gene_enrichment <- function(calls_long, feature_label) {
  
  genes <- sort(unique(calls_long$Gene))
  
  res <- lapply(genes, function(g) {
    gene_patients <- unique(substr(calls_long$Sample[calls_long$Gene == g], 1, 12))
    
    tmp <- ctcf_status_df
    tmp$feature_positive <- tmp$patient %in% gene_patients
    
    cont_table <- table(tmp$CTCF_status, tmp$feature_positive)
    fisher_res <- fisher.test(cont_table)
    
    n_wt  <- sum(tmp$CTCF_status == "CTCF-wt")
    n_mut <- sum(tmp$CTCF_status == "CTCF-mut")
    
    n_feature_wt  <- sum(tmp$feature_positive[tmp$CTCF_status == "CTCF-wt"])
    n_feature_mut <- sum(tmp$feature_positive[tmp$CTCF_status == "CTCF-mut"])
    
    data.frame(
      gene = g,
      feature = feature_label,
      n_CTCF_wt = n_wt,
      n_CTCF_mut = n_mut,
      n_positive_CTCF_wt = n_feature_wt,
      n_positive_CTCF_mut = n_feature_mut,
      pct_positive_CTCF_wt = 100 * n_feature_wt / n_wt,
      pct_positive_CTCF_mut = 100 * n_feature_mut / n_mut,
      p_value = fisher_res$p.value
    )
  })
  
  res <- bind_rows(res) %>%
    mutate(
      FDR = p.adjust(p_value, method = "BH"),
      direction = ifelse(pct_positive_CTCF_mut > pct_positive_CTCF_wt,
                         "enriched_in_CTCF_mut",
                         "depleted_in_CTCF_mut")
    ) %>%
    arrange(FDR, p_value)
  
  res
}

# ------------------------------------------------------------------------------
# ---- ALD gene-level enrichment ----
ald_calls_file <- file.path(res_dir, "ald_calls_long.rds")
if (!file.exists(ald_calls_file)) stop("Missing: ", ald_calls_file)

ald_calls_long <- readRDS(ald_calls_file)

ald_gene_results <- test_gene_enrichment(ald_calls_long, "ALD") %>%
  filter(direction == "enriched_in_CTCF_mut", FDR < 0.05) %>%
  select(-direction)

ald_gene_file <- file.path(
  output_dir,
  paste0("UCEC_CTCF_ALD_gene_enrichment_", gthr_tag, "_", sig_tag, "_", cnv_tag, "_prev_", prev_tag, ald_tpm_suffix, ".csv")
)
write.csv(ald_gene_results, ald_gene_file, row.names = FALSE)
message("Saved: ", ald_gene_file)

# ---- DLD gene-level enrichment ----
dld_calls_file <- file.path(res_dir, "dld_calls_long.rds")
if (!file.exists(dld_calls_file)) stop("Missing: ", dld_calls_file)

dld_calls_long <- readRDS(dld_calls_file)

dld_gene_results <- test_gene_enrichment(dld_calls_long, "DLD") %>%
  filter(direction == "enriched_in_CTCF_mut", FDR < 0.05) %>%
  select(-direction)

dld_gene_file <- file.path(
  output_dir,
  paste0("UCEC_CTCF_DLD_gene_enrichment_", gthr_tag, "_", sig_tag, "_", cnv_tag, "_prev_", prev_tag, ald_tpm_suffix, ".csv")
)
write.csv(dld_gene_results, dld_gene_file, row.names = FALSE)
message("Saved: ", dld_gene_file)

message("ALD enriched genes saved: ", nrow(ald_gene_results))
message("DLD enriched genes saved: ", nrow(dld_gene_results))
