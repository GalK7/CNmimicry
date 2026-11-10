# Purpose:
#   Summarize mutation-type composition in DLD versus other-diploid
#   gene-sample pairs across cohorts.
#
# Inputs:
#   - Core workflow DLD call files and deletion gene lists from result directories
#     defined in config.R.
#   - Prepared copy-number data and cached or downloaded TCGA somatic mutation data.
#
# Outputs:
#   - Mutation-type composition CSV and workbook files written to analysis/outputs/mutation.

source("config.R")

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
source(file.path("R", "prepare_data.R"))

stopifnot(requireNamespace("dplyr", quietly = TRUE))
stopifnot(requireNamespace("readr", quietly = TRUE))
stopifnot(requireNamespace("tidyr", quietly = TRUE))
stopifnot(requireNamespace("tibble", quietly = TRUE))
stopifnot(requireNamespace("TCGAbiolinks", quietly = TRUE))
stopifnot(requireNamespace("openxlsx", quietly = TRUE))

library(dplyr)
library(tidyr)
library(tibble)

# ----------------------------
# Mutation classes
# ----------------------------
silent_classes <- c("Silent")

missense_classes <- c(
  "Missense_Mutation",
  "In_Frame_Del",
  "In_Frame_Ins"
)

truncating_classes <- c(
  "Nonsense_Mutation",
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "Splice_Site",
  "Nonstop_Mutation",
  "Translation_Start_Site"
)

collapse_mutation_class <- function(vc) {
  if (any(vc %in% truncating_classes)) return("truncating")
  if (any(vc %in% missense_classes))   return("missense")
  if (any(vc %in% silent_classes))     return("silent")
  return(NA_character_)
}

# ----------------------------
# MAF caching
# ----------------------------
get_maf_cached_ensemble <- function(gdc_project, cache_dir, gdc_dir = NULL) {
  maf_dir <- file.path(workflow_root, "analysis", "cache", "maf_ensemble_masked")
  maf_gdc_dir <- file.path(workflow_root, "analysis", "GDCdata")
  dir.create(maf_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(maf_gdc_dir, recursive = TRUE, showWarnings = FALSE)
  
  maf_path <- file.path(maf_dir, paste0(gdc_project, "_ensemble_masked_maf.rds"))
  if (file.exists(maf_path)) return(readRDS(maf_path))
  
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
  maf_data
}

# ----------------------------
# Match run_all.R results_dir logic
# ----------------------------
if (isTRUE(use_gistic_thresholded_intersection)) {
  results_dir_base <- results_dir
  results_dir_sig_base <- results_dir_sig
  results_dir_cnvnorm_base <- results_dir_cnvnorm
} else {
  results_dir_base <- results_dir_no_gistic_thr
  results_dir_sig_base <- results_dir_sig_no_gistic_thr
  results_dir_cnvnorm_base <- results_dir_cnvnorm_no_gistic_thr
}

gistic_sig_folder_tag <- if (isTRUE(use_gistic_significant)) {
  paste0("with_SIG_", gistic_sig_gene_source)
} else {
  "no_SIG"
}

results_dir_use <- if (isTRUE(use_cnv_median_normalization)) {
  file.path(results_dir_cnvnorm_base, gistic_sig_folder_tag)
} else {
  if (isTRUE(use_gistic_significant)) results_dir_sig_base else results_dir_base
}
if (ald_tpm_filter_enabled) results_dir_use <- file.path(results_dir_use, ald_tpm_tag)

# ----------------------------
# One cohort-level mutation-type composition summary
# ----------------------------
run_mutation_type_composition_one <- function(cdef, prev_thr, results_dir_use, cache_dir) {
  cohort <- cdef$cohort
  gdc_project <- cdef$gdc_project
  gdc_dir <- cdef$gdc_dir
  
  out_dir <- file.path(results_dir_use, cohort, paste0("prev_", prev_thr))
  
  dld_calls_path <- file.path(out_dir, "dld_calls_long.rds")
  dld_ps_path    <- file.path(out_dir, "dld_per_sample.rds")
  del_genes_path <- file.path(out_dir, "deleted_genes.rds")
  
  if (!file.exists(dld_calls_path)) stop("Missing: ", dld_calls_path)
  if (!file.exists(dld_ps_path))    stop("Missing: ", dld_ps_path)
  if (!file.exists(del_genes_path)) stop("Missing: ", del_genes_path)
  
  dld_per_sample <- readRDS(dld_ps_path)
  dld_calls <- readRDS(dld_calls_path) %>%
    distinct(Sample, Gene)
  
  deleted_genes <- readRDS(del_genes_path)
  
  prepared <- prepare_data(
    cohort = cohort,
    gdc_project = gdc_project,
    gdc_dir = gdc_dir,
    cache_dir = cache_dir,
    er_pos_only = isTRUE(cdef$er_pos_only)
  )
  
  cnv_mat <- prepared$cnv_mat
  
  if (isTRUE(use_cnv_median_normalization)) {
    denom <- apply(cnv_mat, 2, median, na.rm = TRUE)
    denom[is.na(denom) | denom == 0] <- NA_real_
    cnv_mat <- sweep(cnv_mat, 2, denom, "/")
    del_cutoff_use <- cnv_del_cutoff_norm
    dip_center_use <- cnv_diploid_norm
    dip_tol_use <- cnv_diploid_tol_norm
  } else {
    del_cutoff_use <- cnv_del_cutoff
    dip_center_use <- cnv_diploid
    dip_tol_use <- 0
  }
  
  universe_samples <- dld_per_sample$Sample %>% unique()
  universe_samples <- intersect(universe_samples, colnames(cnv_mat))
  deleted_genes <- intersect(deleted_genes, rownames(cnv_mat))
  
  if (length(universe_samples) == 0 || length(deleted_genes) == 0) {
    message(cohort, ": no eligible sample-gene pairs available for mutation-type composition.")
    return(tibble(
      cancer = character(), threshold = numeric(), group = character(),
      mutation_class = character(), n = integer(), total_mut = integer(), pct = numeric()
    ))
  }
  
  is_diploid <- function(x) {
    if (dip_tol_use == 0) {
      x == dip_center_use
    } else {
      abs(x - dip_center_use) <= dip_tol_use
    }
  }
  
  # All eligible sample-gene pairs for deleted genes in this cancer
  cnv_long <- as.data.frame(t(cnv_mat[deleted_genes, universe_samples, drop = FALSE])) %>%
    rownames_to_column("Sample") %>%
    pivot_longer(-Sample, names_to = "Gene", values_to = "cnv_value") %>%
    mutate(
      patient = substr(Sample, 1, 12),
      is_deleted = cnv_value <= del_cutoff_use,
      is_diploid = is_diploid(cnv_value)
    )
  
  dld_pairs <- dld_calls %>%
    distinct(Sample, Gene) %>%
    mutate(is_dld = TRUE)
  
  pair_df <- cnv_long %>%
    left_join(dld_pairs, by = c("Sample", "Gene")) %>%
    mutate(
      is_dld = dplyr::coalesce(is_dld, FALSE),
      group = case_when(
        is_deleted ~ NA_character_,
        is_dld ~ "DLD",
        is_diploid ~ "other_diploid",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(group))
  
  maf_data <- get_maf_cached_ensemble(gdc_project, cache_dir, gdc_dir)
  
  maf_gene_class <- maf_data %>%
    transmute(
      patient = substr(Tumor_Sample_Barcode, 1, 12),
      Gene = Hugo_Symbol,
      Variant_Classification = Variant_Classification
    ) %>%
    group_by(patient, Gene) %>%
    summarise(
      mutation_class = collapse_mutation_class(Variant_Classification),
      .groups = "drop"
    ) %>%
    filter(!is.na(mutation_class))
  
  plot_df <- pair_df %>%
    left_join(maf_gene_class, by = c("patient", "Gene")) %>%
    filter(!is.na(mutation_class)) %>%
    count(group, mutation_class, name = "n") %>%
    complete(
      group = c("DLD", "other_diploid"),
      mutation_class = c("silent", "missense", "truncating"),
      fill = list(n = 0)
    ) %>%
    group_by(group) %>%
    mutate(
      total_mut = sum(n),
      pct = ifelse(total_mut > 0, 100 * n / total_mut, NA_real_)
    ) %>%
    ungroup() %>%
    mutate(
      cancer = cohort,
      threshold = prev_thr
    ) %>%
    select(cancer, threshold, group, mutation_class, n, total_mut, pct)
  
  plot_df
}

# ----------------------------
# Run all cancers / thresholds
# ----------------------------
cancers_to_run <- cancers

output_dir <- file.path(workflow_root, "analysis", "outputs", "mutation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

all_results <- list()

for (thr in prevalence_thresholds) {
  message("\n[mutation_x_DLD_mutation_type_composition] Building outputs for prev_thr=", thr)
  
  wb <- openxlsx::createWorkbook()
  
  for (cancer_name in names(cancers_to_run)) {
    cdef <- cancers_to_run[[cancer_name]]
    message("[mutation_x_DLD_mutation_type_composition]  Running: ", cancer_name, " | prev_thr=", thr)
    
    res <- tryCatch(
      run_mutation_type_composition_one(
        cdef = cdef,
        prev_thr = thr,
        results_dir_use = results_dir_use,
        cache_dir = cache_dir
      ),
      error = function(e) {
        stop("Mutation-type composition failed for ", cancer_name,
             " | prev_thr=", thr, ": ", conditionMessage(e), call. = FALSE)
      }
    )
    
    all_results[[paste0(cancer_name, "_", thr)]] <- res
    openxlsx::addWorksheet(wb, sheetName = cancer_name)
    openxlsx::writeData(wb, sheet = cancer_name, x = res)
  }
  
  gthr_tag <- if (isTRUE(use_gistic_thresholded_intersection)) "with_GISTICthr" else "no_GISTICthr"
  sig_tag  <- gistic_sig_folder_tag
  cnv_tag  <- if (isTRUE(use_cnv_median_normalization)) "CNVnorm" else "rawCNV"
  
  out_xlsx <- file.path(
    output_dir,
    paste0("mutation_x_DLD_mutation_type_composition_", thr, "_", gthr_tag, "_", sig_tag, "_", cnv_tag, ".xlsx")
  )
  
  openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
  message("[mutation_x_DLD_mutation_type_composition] Saved: ", out_xlsx)
}

combined_df <- bind_rows(all_results)

gthr_tag <- if (isTRUE(use_gistic_thresholded_intersection)) "with_GISTICthr" else "no_GISTICthr"
sig_tag  <- gistic_sig_folder_tag
cnv_tag  <- if (isTRUE(use_cnv_median_normalization)) "CNVnorm" else "rawCNV"

out_csv <- file.path(
  output_dir,
  paste0("mutation_x_DLD_mutation_type_composition_all_", gthr_tag, "_", sig_tag, "_", cnv_tag, ".csv")
)

readr::write_csv(combined_df, out_csv)
message("Saved combined CSV: ", out_csv)
