################################################################################
# Purpose
# Run the ALD promoter methylation screen across TCGA cohorts.
#
# Inputs
# - CNmimicry core outputs for ALD calls and ALD-eligible gene sets.
# - TCGA HM450 primary-tumor methylation data, downloaded or loaded through
#   analysis/methylation/R/methylation_helpers.R.
#
# Outputs
# - Per-cohort ALD promoter methylation Wilcoxon result files.
# - Per-cohort input-overlap summaries and promoter methylation matrix caches.
# - Files are written under analysis/methylation/outputs.
################################################################################

source("config.R")
source(file.path("analysis", "methylation", "R", "methylation_helpers.R"))

stopifnot(requireNamespace("TCGAbiolinks", quietly = TRUE))
stopifnot(requireNamespace("SummarizedExperiment", quietly = TRUE))
stopifnot(requireNamespace("dplyr", quietly = TRUE))
stopifnot(requireNamespace("readr", quietly = TRUE))
stopifnot(requireNamespace("minfi", quietly = TRUE))
stopifnot(requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE))

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(readr)

# ------------------------------------------------------------------------------
# Run settings
# ------------------------------------------------------------------------------
cohorts <- names(cancers)

prev_thr <- prevalence_thresholds[1]

promoter_regions <- c("TSS200", "TSS1500", "5'UTR")
methylation_platform <- "Illumina Human Methylation 450"
methylation_sample_type <- "Primary Tumor"
methylation_setting_tag <- paste(
  make_promoter_region_tag(promoter_regions),
  "HM450",
  "PrimaryTumor",
  sep = "_"
)

min_ald <- 10
min_ref <- 10
delta_methylation_cutoff <- 0.20
sample_filter_tag <- paste0("minALD", min_ald, "_minRef", min_ref)
delta_filter_tag <- paste0(
  "delta",
  gsub("\\.", "p", format(delta_methylation_cutoff, nsmall = 2))
)

force_rebuild_methylation_matrix <- FALSE

# ------------------------------------------------------------------------------
# Run ALD promoter methylation screen
# ------------------------------------------------------------------------------
for (cohort in cohorts) {
  
  message("Running ALD promoter methylation screen for ", cohort, " ...")
  
  out_dir <- resolve_methylation_output_dir(cohort)
  
  ctx <- resolve_methylation_context(cohort, prev_thr, apply_ald_tpm_filter = TRUE)
  
  out_file <- file.path(
    out_dir,
    paste0(
      cohort,
      "_ALD_promoter_methylation_wilcox_",
      ctx$gthr_tag, "_",
      ctx$sig_tag, "_",
      ctx$cnv_tag,
      "_methylation_", methylation_setting_tag,
      "_", sample_filter_tag,
      "_", delta_filter_tag,
      "_prev_", format(prev_thr, nsmall = 2),
      ctx$ald_tpm_suffix,
      ".csv"
    )
  )
  
  input_summary_file <- file.path(
    out_dir,
    paste0(
      cohort,
      "_ALD_promoter_methylation_input_summary_",
      ctx$gthr_tag, "_",
      ctx$sig_tag, "_",
      ctx$cnv_tag,
      "_methylation_", methylation_setting_tag,
      "_", sample_filter_tag,
      "_", delta_filter_tag,
      "_prev_", format(prev_thr, nsmall = 2),
      ctx$ald_tpm_suffix,
      ".csv"
    )
  )
  
  promoter_obj <- tryCatch(
    build_promoter_methylation_matrix(
      cohort = cohort,
      promoter_regions = promoter_regions,
      platform = methylation_platform,
      sample_type = methylation_sample_type,
      force_rebuild = force_rebuild_methylation_matrix
    ),
    error = function(e) {
      message("[methylation] Skipping ", cohort, ": ", conditionMessage(e))
      
      res_df <- data.frame(
        gene = character(),
        n_ald = integer(),
        n_other_diploid = integer(),
        median_methylation_ald = numeric(),
        median_methylation_other = numeric(),
        delta_methylation = numeric(),
        p_value = numeric(),
        FDR = numeric(),
        passes_fdr_delta = logical(),
        stringsAsFactors = FALSE
      )
      
      input_summary <- data.frame(
        cohort = cohort,
        prev_thr = prev_thr,
        gistic_thresholded_intersection = ctx$gthr_tag,
        gistic_significant = ctx$sig_tag,
        cnv_mode = ctx$cnv_tag,
        ald_tpm_filter_enabled = ctx$ald_tpm_filter_enabled,
        ald_min_amplified_median_tpm = if (ctx$ald_tpm_filter_enabled) ald_min_amplified_median_tpm else NA_real_,
        ald_tpm_tag = ctx$ald_tpm_tag,
        promoter_regions = paste(promoter_regions, collapse = ";"),
        methylation_platform = methylation_platform,
        methylation_sample_type = methylation_sample_type,
        min_ald = min_ald,
        min_ref = min_ref,
        delta_methylation_cutoff = delta_methylation_cutoff,
        methylation_status = "no_matching_methylation_data",
        methylation_message = conditionMessage(e),
        n_methylation_patients = 0L,
        n_cnv_patients = NA_integer_,
        n_ald_call_patients = NA_integer_,
        n_methylation_cnv_overlap_patients = 0L,
        n_methylation_ald_call_overlap_patients = 0L,
        n_amplified_genes_in_cnmimicry_results = NA_integer_,
        n_ald_eligible_genes_in_cnmimicry_results = NA_integer_,
        n_genes_with_promoter_methylation = 0L,
        n_amplified_genes_with_promoter_methylation = 0L,
        n_ald_eligible_genes_with_promoter_methylation = 0L,
        n_genes_tested_after_min_sample_filters = 0L,
        stringsAsFactors = FALSE
      )
      
      write.csv(res_df, out_file, row.names = FALSE)
      message("Saved empty result: ", out_file)
      write.csv(input_summary, input_summary_file, row.names = FALSE)
      message("Saved input summary: ", input_summary_file)
      
      NULL
    }
  )
  
  if (is.null(promoter_obj)) {
    next
  }
  
  promoter_mat <- promoter_obj$promoter_methylation_matrix
  
  message("Promoter methylation matrix dimensions:")
  print(dim(promoter_mat))
  
  ald_inputs <- load_ald_methylation_inputs(
    cohort = cohort,
    prev_thr = prev_thr
  )
  
  res_df <- run_ald_promoter_methylation_screen(
    promoter_mat = promoter_mat,
    ald_inputs = ald_inputs,
    min_ald = min_ald,
    min_ref = min_ref
  )
  
  input_summary <- summarize_ald_methylation_inputs(
    promoter_mat = promoter_mat,
    ald_inputs = ald_inputs,
    res_df = res_df,
    promoter_regions = promoter_regions,
    platform = methylation_platform,
    sample_type = methylation_sample_type,
    min_ald = min_ald,
    min_ref = min_ref,
    delta_methylation_cutoff = delta_methylation_cutoff
  )
  
  if (nrow(res_df) > 0) {
    res_df <- res_df %>%
      mutate(
        passes_fdr_delta = FDR < 0.05 &
          abs(delta_methylation) >= delta_methylation_cutoff
      )
    
    message("Genes tested: ", nrow(res_df))
    message("FDR < 0.05: ", sum(res_df$FDR < 0.05, na.rm = TRUE))
    message("FDR < 0.05 and |delta methylation| >= ", delta_methylation_cutoff, ": ",
            sum(res_df$passes_fdr_delta, na.rm = TRUE))
    message("FDR < 0.05 and delta methylation >= ", delta_methylation_cutoff, ": ",
            sum(res_df$FDR < 0.05 & res_df$delta_methylation >= delta_methylation_cutoff, na.rm = TRUE))
    message("FDR < 0.05 and delta methylation <= -", delta_methylation_cutoff, ": ",
            sum(res_df$FDR < 0.05 & res_df$delta_methylation <= -delta_methylation_cutoff, na.rm = TRUE))
  } else {
    message("No valid genes tested.")
    res_df <- data.frame(
      gene = character(), n_ald = integer(), n_other_diploid = integer(),
      median_methylation_ald = numeric(), median_methylation_other = numeric(),
      delta_methylation = numeric(), p_value = numeric(), FDR = numeric(),
      passes_fdr_delta = logical(), stringsAsFactors = FALSE
    )
  }
  
  write.csv(res_df, out_file, row.names = FALSE)
  message("Saved: ", out_file)
  
  write.csv(input_summary, input_summary_file, row.names = FALSE)
  message("Saved input summary: ", input_summary_file)
}
