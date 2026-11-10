# Workflow/make_excel.R

source("config.R")

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Please install openxlsx: install.packages('openxlsx')")
}

# Choose which results folder to scan based on the thresholded-GISTIC-intersection toggle
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

results_dir <- if (isTRUE(use_cnv_median_normalization)) {
  file.path(results_dir_cnvnorm_base, gistic_sig_folder_tag)
} else {
  if (isTRUE(use_gistic_significant)) results_dir_sig_base else results_dir_base
}


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
if (ald_tpm_filter_enabled) results_dir <- file.path(results_dir, ald_tpm_tag)

# Read only the configured cohort/threshold paths, not nested runs of other TPM settings.
summary_files <- unlist(lapply(names(cancers), function(cohort) {
  file.path(results_dir, cohort, paste0("prev_", format(prevalence_thresholds, nsmall = 2)), "summary.csv")
}), use.names = FALSE)
missing_files <- summary_files[!file.exists(summary_files)]
if (length(missing_files) > 0) {
  stop("Missing core summaries. Run run_all.R first:\n", paste(missing_files, collapse = "\n"))
}

df <- do.call(rbind, lapply(summary_files, function(f) read.csv(f, stringsAsFactors = FALSE)))

if (ald_tpm_filter_enabled) {
  required <- c("ald_tpm_filter_enabled", "ald_min_amplified_median_tpm",
                "ald_tpm_tag", "n_ald_eligible_genes", "ald_gene_denominator")
  if (!all(required %in% names(df)) ||
      !isTRUE(all(df$ald_tpm_filter_enabled &
                  df$ald_min_amplified_median_tpm == ald_min_amplified_median_tpm &
                  df$ald_tpm_tag == ald_tpm_tag &
                  df$ald_gene_denominator == "ald_eligible_genes"))) {
    stop("Core summary ALD TPM settings do not match config.R. Run run_all.R for these settings.")
  }
}

# Ensure numeric thresholds (0.1 vs 0.10)
df$threshold <- as.numeric(df$threshold)

# Preserve the cohort order defined in config.R
cancer_order <- names(cancers)

make_threshold_sheet_name <- function(thr) {
  paste0(format(thr * 100, trim = TRUE, scientific = FALSE), "%")
}

build_table <- function(thr, mode = c("ALD", "DLD")) {
  mode <- match.arg(mode)
  
  
  sub <- df[df$threshold == thr, , drop = FALSE]
  sub <- sub[match(cancer_order, sub$cancer), , drop = FALSE]
  
  if (any(is.na(sub$cancer))) {
    stop("Missing cancers for threshold ", thr, ". Did some runs fail or were not produced?")
  }
  
  if (mode == "ALD") {
    out <- rbind(
      "Number of amplified genes"                        = sub$n_amp_genes,
      "Number of samples ALD in >5% of amplified genes"  = sub$n_samples_gt5pct_ald,
      "Number of genes with >20% ALD"                    = sub$n_genes_gt20_ald,
      "Top 10 genes (by ALD %)"                          = sub$top10_genes_ald,
      "Cancer Score (>20%)"                              = sub$cancer_score20_ald,
      "Cancer Score (>50%)"                              = sub$cancer_score50_ald
    )
    if (ald_tpm_filter_enabled) {
      out <- rbind(out[1, , drop = FALSE],
                   "Number of ALD-eligible amplified genes (TPM filter)" = sub$n_ald_eligible_genes,
                   out[-1, , drop = FALSE])
      rownames(out)[3] <- "Number of samples ALD in >5% of ALD-eligible genes"
    }
  } else {
    out <- rbind(
      "Number of deleted genes"                          = sub$n_del_genes,
      "Number of samples DLD in >5% of deleted genes"    = sub$n_samples_gt5pct_dld,
      "Number of genes with >20% DLD"                    = sub$n_genes_gt20_dld,
      "Top 10 genes (by DLD %)"                          = sub$top10_genes_dld,
      "Cancer Score (>20%)"                              = sub$cancer_score20_dld,
      "Cancer Score (>50%)"                              = sub$cancer_score50_dld
    )
  }
  
  colnames(out) <- cancer_order
  out
}

make_settings_table <- function(mode, output_label) {
  data.frame(
    setting = c(
      "output_label",
      "internal_mode",
      "workflow_root",
      "source_results_dir",
      "number_of_summary_files",
      "number_of_cancers",
      "cancers",
      "prevalence_thresholds",
      "use_cnv_median_normalization",
      "cnv_norm_stat",
      "use_gistic_thresholded_intersection",
      "use_gistic_significant",
      "gistic_sig_gene_source",
      "gistic_sig_q_cutoff",
      "ald_tpm_filter_enabled",
      "ald_min_amplified_median_tpm",
      "ald_tpm_tag",
      "gene_denominator_for_burden_and_scores",
      "generated_at"
    ),
    value = c(
      output_label,
      mode,
      workflow_root,
      normalizePath(results_dir, winslash = "/", mustWork = FALSE),
      length(summary_files),
      length(cancer_order),
      paste(cancer_order, collapse = ", "),
      paste(prevalence_thresholds, collapse = ", "),
      as.character(use_cnv_median_normalization),
      cnv_norm_stat,
      as.character(use_gistic_thresholded_intersection),
      as.character(use_gistic_significant),
      if (isTRUE(use_gistic_significant)) gistic_sig_gene_source else NA_character_,
      as.character(gistic_sig_q_cutoff),
      as.character(ald_tpm_filter_enabled),
      if (ald_tpm_filter_enabled) as.character(ald_min_amplified_median_tpm) else "disabled",
      ald_tpm_tag,
      if (mode == "DLD") "deleted genes" else if (ald_tpm_filter_enabled) "ALD-eligible amplified genes (median amplified TPM strictly above the configured minimum)" else "amplified genes",
      as.character(Sys.time())
    ),
    stringsAsFactors = FALSE
  )
}

write_workbook <- function(mode = c("ALD", "DLD"), output_label, out_file) {
  mode <- match.arg(mode)
  wb <- openxlsx::createWorkbook()

  openxlsx::addWorksheet(wb, "settings")
  openxlsx::writeData(wb, "settings", make_settings_table(mode, output_label))
  openxlsx::freezePane(wb, "settings", firstRow = TRUE)
  openxlsx::setColWidths(wb, "settings", cols = 1:2, widths = "auto")
  
  for (thr in prevalence_thresholds) {
    sh <- make_threshold_sheet_name(thr)
    
    tab <- build_table(thr, mode = mode)
    tab_df <- data.frame(
      Summary = rownames(tab),
      as.data.frame(tab, stringsAsFactors = FALSE, check.names = FALSE),
      check.names = FALSE
    )
    
    openxlsx::addWorksheet(wb, sh)
    openxlsx::writeData(wb, sh, tab_df, rowNames = FALSE)
    openxlsx::setColWidths(wb, sh, cols = 1:ncol(tab_df), widths = "auto")
  }
  
  openxlsx::saveWorkbook(wb, out_file, overwrite = TRUE)
  message("Saved: ", out_file)
}

suffix <- if (isTRUE(use_gistic_thresholded_intersection)) {
  if (isTRUE(use_gistic_significant)) paste0("with_GISTIC_with_SIG_", gistic_sig_gene_source) else "with_GISTIC_no_SIG"
} else {
  if (isTRUE(use_gistic_significant)) paste0("no_GISTICthr_with_SIG_", gistic_sig_gene_source) else "no_GISTICthr_no_SIG"
}

if (isTRUE(use_cnv_median_normalization) && cnv_norm_stat == "median") {
  suffix <- paste0(suffix, "_CNVnorm_median")
} else if (isTRUE(use_cnv_median_normalization) && cnv_norm_stat == "mean") {
  suffix <- paste0(suffix, "_CNVnorm_mean")
}

if (ald_tpm_filter_enabled) suffix <- paste0(suffix, "_", ald_tpm_tag)

excel_dir <- file.path(workflow_root, "Excel", suffix)
dir.create(excel_dir, showWarnings = FALSE, recursive = TRUE)

write_workbook("ALD", "ALD", file.path(excel_dir, paste0("ALD_summary_", suffix, ".xlsx")))
write_workbook("DLD", "DLD", file.path(excel_dir, paste0("DLD_summary_", suffix, ".xlsx")))
