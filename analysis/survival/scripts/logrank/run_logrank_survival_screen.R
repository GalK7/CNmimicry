# Purpose: Screen ALD/DLD and true AMP/DEL associations with overall survival across cohorts.
# Inputs: Core calls/gene lists, prepared expression/CNV data, and clinical OS.
# Outputs: Per-cohort results and summaries under analysis/survival/outputs/logrank,
#   with a genes_in_peak_version subfolder for peak-restricted runs and an aldTPMgt<value>
#   subfolder when TPM filtering is enabled.

source("config.R")
source(file.path("analysis", "survival", "R", "survival_helpers.R"))
source(file.path("analysis", "survival", "R", "logrank_helpers.R"))

stopifnot(requireNamespace("TCGAbiolinks", quietly = TRUE))
stopifnot(requireNamespace("dplyr", quietly = TRUE))
stopifnot(requireNamespace("survival", quietly = TRUE))
stopifnot(requireNamespace("readr", quietly = TRUE))

library(TCGAbiolinks)
library(dplyr)
library(survival)
library(readr)

# ------------------------------------------------------------------------------
# Run settings
# ------------------------------------------------------------------------------
cohorts <- names(cancers)
modes <- c("DLD", "DEL", "ALD", "AMP")
prev_thr <- prevalence_thresholds[1]
ald_tpm <- get_survival_ald_tpm_settings()
min_pos <- 10
min_ref <- 10
sample_filter_tag <- paste0("minPos", min_pos, "_minRef", min_ref)
fdr_cutoffs <- c(0.05, 0.10)
force_rebuild_survival_cache <- FALSE

logrank_output_root <- resolve_logrank_output_root()
dir.create(logrank_output_root, recursive = TRUE, showWarnings = FALSE)

sig_tag <- if (isTRUE(use_gistic_significant)) {
  paste0("with_SIG_", gistic_sig_gene_source)
} else {
  "no_SIG"
}
gthr_tag <- if (isTRUE(use_gistic_thresholded_intersection)) "with_GISTICthr" else "no_GISTICthr"
cnv_tag <- if (isTRUE(use_cnv_median_normalization)) "cnvnorm" else "rawcnv"
prev_tag <- format(prev_thr, nsmall = 2)

# ------------------------------------------------------------------------------
# Run log-rank screens
# ------------------------------------------------------------------------------
for (cohort in cohorts) {
  message("[log-rank] Preparing shared inputs for ", cohort, " ...")
  core_inputs <- build_survival_core_inputs(
    cohort = cohort,
    prev_thr = prev_thr,
    force_rebuild_survival_cache = force_rebuild_survival_cache
  )
  inputs_deleted <- attach_deleted_inputs(core_inputs)
  inputs_amplified <- attach_amplified_inputs(core_inputs)

  cohort_out_dir <- file.path(logrank_output_root, cohort)
  dir.create(cohort_out_dir, recursive = TRUE, showWarnings = FALSE)

  for (mode in modes) {
    message("[log-rank] Running ", mode, " log-rank screen for ", cohort, " ...")

    mode_df <- run_logrank_screen(
      mode = mode,
      inputs_deleted = inputs_deleted,
      inputs_amplified = inputs_amplified,
      min_pos = min_pos,
      min_ref = min_ref
    )

    if (nrow(mode_df) == 0) {
      mode_df <- empty_logrank_df(mode)
    }

    comparison_label <- switch(
      mode,
      DLD = "DLD_vs_otherDiploids",
      DEL = "DEL_vs_otherDiploidsExcludingDLD",
      ALD = "ALD_vs_otherDiploids",
      AMP = "AMP_vs_otherDiploidsExcludingALD",
      stop("Unsupported mode: ", mode)
    )

    out_file <- file.path(
      cohort_out_dir,
      paste0(
        cohort, "_", comparison_label, "_OS_logrank_",
        gthr_tag, "_", sig_tag, "_", cnv_tag,
        "_", sample_filter_tag,
        "_prev_", prev_tag,
        ald_tpm$suffix,
        ".csv"
      )
    )

    write.csv(mode_df, out_file, row.names = FALSE)

    message("[log-rank] Saved: ", out_file)
    message("[log-rank] Tested genes: ", nrow(mode_df))
    message("[log-rank] FDR < 0.05: ", count_logrank_sig(mode_df, 0.05))
  }
}

# ------------------------------------------------------------------------------
# Cross-cancer summaries
# ------------------------------------------------------------------------------
for (fdr_cutoff in fdr_cutoffs) {
  fdr_dirname <- paste0("fdr_", format(fdr_cutoff, nsmall = 2))
  fdr_out_dir <- file.path(logrank_output_root, fdr_dirname)
  dir.create(fdr_out_dir, recursive = TRUE, showWarnings = FALSE)

  summary_df <- lapply(cohorts, function(cohort) {
    cohort_out_dir <- file.path(logrank_output_root, cohort)

    read_mode <- function(mode) {
      comparison_label <- switch(
        mode,
        DLD = "DLD_vs_otherDiploids",
        DEL = "DEL_vs_otherDiploidsExcludingDLD",
        ALD = "ALD_vs_otherDiploids",
        AMP = "AMP_vs_otherDiploidsExcludingALD",
        stop("Unsupported mode: ", mode)
      )

      safe_read_logrank_csv(file.path(
        cohort_out_dir,
        paste0(
          cohort, "_", comparison_label, "_OS_logrank_",
          gthr_tag, "_", sig_tag, "_", cnv_tag,
          "_", sample_filter_tag,
          "_prev_", prev_tag,
          ald_tpm$suffix,
          ".csv"
        )
      ))
    }

    dld_df <- read_mode("DLD")
    del_df <- read_mode("DEL")
    ald_df <- read_mode("ALD")
    amp_df <- read_mode("AMP")

    data.frame(
      cohort = cohort,

      dld_tested = nrow(dld_df),
      dld_sig = count_logrank_sig(dld_df, fdr_cutoff),
      dld_logrank_non_ok = count_logrank_non_ok(dld_df),

      del_tested = nrow(del_df),
      del_sig = count_logrank_sig(del_df, fdr_cutoff),
      del_logrank_non_ok = count_logrank_non_ok(del_df),

      ald_tested = nrow(ald_df),
      ald_sig = count_logrank_sig(ald_df, fdr_cutoff),
      ald_logrank_non_ok = count_logrank_non_ok(ald_df),

      amp_tested = nrow(amp_df),
      amp_sig = count_logrank_sig(amp_df, fdr_cutoff),
      amp_logrank_non_ok = count_logrank_non_ok(amp_df),

      total_sig = count_logrank_sig(dld_df, fdr_cutoff) +
        count_logrank_sig(del_df, fdr_cutoff) +
        count_logrank_sig(ald_df, fdr_cutoff) +
        count_logrank_sig(amp_df, fdr_cutoff),
      stringsAsFactors = FALSE
    )
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::arrange(desc(total_sig), cohort)

  out_file <- file.path(
    fdr_out_dir,
    paste0(
      "survival_logrank_crossCancer_summary_",
      gthr_tag, "_", sig_tag, "_", cnv_tag,
      "_", sample_filter_tag,
      "_prev_", prev_tag,
      ald_tpm$suffix,
      ".csv"
    )
  )

  write.csv(summary_df, out_file, row.names = FALSE)
  message("[log-rank] Saved summary: ", out_file)
}
