# Purpose: Compare AMP versus diploid and AMP+ALD versus remaining diploid survival results across cohorts.
# Inputs: Default AMP results, Supplementary Table 12, core sample groups, and clinical OS.
# Outputs: Workbook under analysis/survival/outputs/logrank, in the configured
#   aldTPMgt<value>/AMP_reference_group_sensitivity subfolder when TPM filtering is enabled,
#   or AMP_reference_group_sensitivity when disabled.

source("config.R")
source(file.path("analysis", "survival", "R", "survival_helpers.R"))
source(file.path("analysis", "survival", "R", "logrank_helpers.R"))

stopifnot(requireNamespace("dplyr", quietly = TRUE))
stopifnot(requireNamespace("openxlsx", quietly = TRUE))
stopifnot(requireNamespace("tidyr", quietly = TRUE))

library(dplyr)

# ------------------------------------------------------------------------------
# AMP reference-group sensitivity analysis
# ------------------------------------------------------------------------------
cohorts <- names(cancers)
prev_thr <- prevalence_thresholds[1]
ald_tpm <- get_survival_ald_tpm_settings()
min_pos <- 10
min_ref <- 10
fdr_cutoff <- 0.05

if (!isTRUE(use_gistic_significant) || !identical(gistic_sig_gene_source, "wide_peak")) {
  stop("This sensitivity analysis must be run with the default wide_peak GISTIC setting.")
}

sample_filter_tag <- paste0("minPos", min_pos, "_minRef", min_ref)
sig_tag <- if (isTRUE(use_gistic_significant)) {
  paste0("with_SIG_", gistic_sig_gene_source)
} else {
  "no_SIG"
}
gthr_tag <- if (isTRUE(use_gistic_thresholded_intersection)) "with_GISTICthr" else "no_GISTICthr"
cnv_tag <- if (isTRUE(use_cnv_median_normalization)) "cnvnorm" else "rawcnv"
prev_tag <- format(prev_thr, nsmall = 2)

main_output_root <- resolve_logrank_output_root()
output_dir <- file.path(main_output_root, "AMP_reference_group_sensitivity")
output_tag <- paste(sig_tag, gthr_tag, cnv_tag, paste0("prev", prev_tag), sep = "_")
output_tag <- paste0(output_tag, ald_tpm$suffix)
output_file <- file.path(output_dir, paste0("AMP_reference_sensitivity_", output_tag, ".xlsx"))

figure7_table <- file.path(
  workflow_root,
  "plots",
  "Figure 7 - ALD DLD survival associations across cancers",
  paste0(
    "Supplementary_Table_12_logrank_",
    gthr_tag, "_", sig_tag, "_", cnv_tag, "_prev", prev_tag,
    "_min", min_pos, "_", min_ref, "_fdr", format(fdr_cutoff, nsmall = 2), ald_tpm$suffix, ".xlsx"
  )
)

current_amp_file <- function(cohort) {
  file.path(
    main_output_root,
    cohort,
    paste0(
      cohort, "_AMP_vs_otherDiploidsExcludingALD_OS_logrank_",
      gthr_tag, "_", sig_tag, "_", cnv_tag, "_",
      sample_filter_tag, "_prev_", prev_tag, ald_tpm$suffix, ".csv"
    )
  )
}

missing_current_files <- vapply(cohorts, current_amp_file, character(1))
missing_current_files <- missing_current_files[!file.exists(missing_current_files)]
if (length(missing_current_files) > 0) {
  stop(
    "Missing default AMP-excluding-ALD log-rank output(s):\n",
    paste(missing_current_files, collapse = "\n")
  )
}
if (!file.exists(figure7_table)) {
  stop("Missing matching Figure 7 survival table used for AMP directions: ", figure7_table,
       "\nFor a TPM-filtered run, update and regenerate Table 12 before running this follow-up.")
}

# ------------------------------------------------------------------------------
# Require Table 12 to match this run before performing AMP analyses.
# ------------------------------------------------------------------------------
validate_table9a_settings <- function(settings, expected) {
  if (!all(c("setting", "value") %in% names(settings))) {
    stop("Table 12 Settings must contain the columns 'setting' and 'value'.",
         call. = FALSE)
  }
  keys <- as.character(settings$setting)
  if (anyNA(keys) || any(!nzchar(keys)) || anyDuplicated(keys)) {
    stop("Table 12 Settings contains missing or duplicate setting names.",
         call. = FALSE)
  }
  missing <- setdiff(names(expected), keys)
  if (length(missing)) {
    stop("Table 12 Settings is missing: ", paste(missing, collapse = ", "),
         ". Regenerate Table 12 with the Figure 7 builder.", call. = FALSE)
  }
  problems <- character()
  for (key in names(expected)) {
    stored <- as.character(settings$value[match(key, keys)])
    wanted <- expected[[key]]
    same <- if (is.numeric(wanted)) {
      number <- suppressWarnings(as.numeric(stored))
      is.finite(number) && isTRUE(all.equal(number, wanted, tolerance = 1e-12))
    } else {
      !is.na(stored) && identical(stored, as.character(wanted))
    }
    if (!same) {
      problems <- c(problems, paste0(key, ": Table 12=", stored,
                                    "; current run=", as.character(wanted)))
    }
  }
  if (length(problems)) {
    stop("Table 12 settings do not match this AMP run:\n",
         paste(problems, collapse = "\n"),
         "\nUse matching run settings and regenerate Table 12 as needed.",
         call. = FALSE)
  }
  invisible(TRUE)
}

if (!"Settings" %in% openxlsx::getSheetNames(figure7_table)) {
  stop("Table 12 has no Settings sheet. Regenerate it with the Figure 7 builder.",
       call. = FALSE)
}
table9a_settings <- openxlsx::read.xlsx(figure7_table, sheet = "Settings")
expected_table9a_settings <- list(
  use_gistic_significant = use_gistic_significant,
  gistic_sig_gene_source = gistic_sig_gene_source,
  use_gistic_thresholded_intersection = use_gistic_thresholded_intersection,
  use_cnv_median_normalization = use_cnv_median_normalization,
  prevalence_threshold = prev_thr,
  min_group = min_pos,
  min_reference = min_ref,
  fdr_cutoff = fdr_cutoff
)
if (ald_tpm$enabled) {
  expected_table9a_settings <- c(expected_table9a_settings, list(
    ald_tpm_filter_enabled = TRUE,
    ald_min_amplified_median_tpm = ald_min_amplified_median_tpm,
    ald_tpm_tag = ald_tpm$tag
  ))
}
validate_table9a_settings(table9a_settings, expected_table9a_settings)

run_amp_reference_screen <- function(inputs, comparison, min_pos, min_ref,
                                     genes_to_test = inputs$amplified_genes_to_test) {
  results <- lapply(genes_to_test, function(gene) {
    if (!gene %in% rownames(inputs$cnv_use)) return(NULL)

    gene_cnv <- inputs$cnv_use[gene, , drop = TRUE]
    amp_patients <- names(gene_cnv)[gene_cnv >= inputs$amp_cutoff]
    amp_patients <- unique(substr(amp_patients, 1, 12))

    diploid_patients <- names(gene_cnv)[inputs$is_diploid(gene_cnv)]
    diploid_patients <- unique(substr(diploid_patients, 1, 12))

    ald_patients <- inputs$ald_calls_patient %>%
      filter(Gene == gene) %>%
      pull(patient) %>%
      unique()

    if (comparison == "AMP_vs_diploids_excluding_ALD") {
      group_patients <- amp_patients
      reference_patients <- setdiff(diploid_patients, ald_patients)
    } else if (comparison == "AMP_vs_all_diploids") {
      group_patients <- amp_patients
      reference_patients <- diploid_patients
    } else if (comparison == "AMP_plus_ALD_vs_other_diploids") {
      group_patients <- union(amp_patients, ald_patients)
      reference_patients <- setdiff(diploid_patients, ald_patients)
    } else {
      stop("Unsupported comparison: ", comparison)
    }

    result <- run_logrank_binary(
      pos_patients = group_patients,
      ref_patients = reference_patients,
      surv_df = inputs$surv_df,
      min_pos = min_pos,
      min_ref = min_ref
    )
    if (is.null(result)) return(NULL)

    data.frame(
      gene = gene,
      n_group = result$n_group,
      n_reference = result$n_ref,
      logrank_direction = result$logrank_direction,
      logrank_p_value = result$logrank_p_value,
      logrank_status = result$logrank_status,
      stringsAsFactors = FALSE
    )
  }) %>%
    bind_rows()

  if (nrow(results) == 0) return(results)

  results %>%
    mutate(FDR = p.adjust(logrank_p_value, method = "BH")) %>%
    arrange(FDR, logrank_p_value)
}

read_current_amp <- function(cohort) {
  current <- read.csv(current_amp_file(cohort), stringsAsFactors = FALSE)
  required_cols <- c("gene", "logrank_p_value", "FDR")
  missing_cols <- setdiff(required_cols, colnames(current))
  if (length(missing_cols) > 0) {
    stop(
      "Default AMP-excluding-ALD file for ", cohort, " is missing: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  current %>%
    transmute(
      cohort = cohort,
      gene = as.character(gene),
      comparison = "AMP_excluding_ALD",
      logrank_direction = NA_character_,
      logrank_p_value = as.numeric(logrank_p_value),
      FDR = as.numeric(FDR)
    )
}

message("[AMP reference sensitivity] Reading default AMP-excluding-ALD results ...")
current_results <- bind_rows(lapply(cohorts, read_current_amp))

figure7_data <- openxlsx::read.xlsx(figure7_table, sheet = "data")

current_direction_table <- figure7_data %>%
  filter(state == "AMP") %>%
  transmute(
    cohort = as.character(cohort),
    gene = as.character(gene),
    logrank_direction = as.character(logrank_direction)
  ) %>%
  distinct(cohort, gene, .keep_all = TRUE)

current_results <- current_results %>%
  select(-logrank_direction) %>%
  left_join(current_direction_table, by = c("cohort", "gene"))

new_results <- lapply(cohorts, function(cohort) {
  message("[AMP reference sensitivity] Preparing ", cohort, " ...")
  core_inputs <- build_survival_core_inputs(
    cohort = cohort,
    prev_thr = prev_thr,
    force_rebuild_survival_cache = FALSE
  )
  inputs <- attach_amplified_inputs(core_inputs)

  bind_rows(lapply(
    c("AMP_vs_all_diploids", "AMP_plus_ALD_vs_other_diploids"),
    function(comparison) {
      message("[AMP reference sensitivity] Running ", cohort, ": ", comparison)
      result <- run_amp_reference_screen(
        inputs = inputs,
        comparison = comparison,
        min_pos = min_pos,
        min_ref = min_ref
      )
      if (nrow(result) == 0) return(result)
      mutate(result, cohort = cohort, comparison = comparison, .before = gene)
    }
  ))
}) %>%
  bind_rows()

all_results <- bind_rows(
  current_results,
  new_results %>%
    select(cohort, gene, comparison, logrank_direction, logrank_p_value, FDR)
) %>%
  mutate(significant = !is.na(FDR) & FDR < fdr_cutoff)

comparison_order <- c(
  "AMP_excluding_ALD",
  "AMP_vs_all_diploids",
  "AMP_plus_ALD_vs_other_diploids"
)

summary_one <- function(df, comparison) {
  x <- df[df$comparison == comparison, , drop = FALSE]
  significant <- x$significant
  data.frame(
    tested = nrow(x),
    significant = sum(significant, na.rm = TRUE),
    worse = sum(significant & x$logrank_direction == "worse", na.rm = TRUE),
    better = sum(significant & x$logrank_direction == "better", na.rm = TRUE),
    not_assigned = sum(significant & x$logrank_direction == "not_assigned", na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

cross_cancer_summary <- lapply(cohorts, function(cohort) {
  cohort_results <- all_results[all_results$cohort == cohort, , drop = FALSE]
  pieces <- lapply(comparison_order, function(comparison) {
    summary_one(cohort_results, comparison)
  })

  data.frame(
    cohort = cohort,
    amp_excluding_ald_tested = pieces[[1]]$tested,
    amp_excluding_ald_significant = pieces[[1]]$significant,
    amp_excluding_ald_worse = pieces[[1]]$worse,
    amp_excluding_ald_better = pieces[[1]]$better,
    amp_all_diploids_tested = pieces[[2]]$tested,
    amp_all_diploids_significant = pieces[[2]]$significant,
    amp_all_diploids_worse = pieces[[2]]$worse,
    amp_all_diploids_better = pieces[[2]]$better,
    amp_all_diploids_change = pieces[[2]]$significant - pieces[[1]]$significant,
    amp_plus_ald_tested = pieces[[3]]$tested,
    amp_plus_ald_significant = pieces[[3]]$significant,
    amp_plus_ald_worse = pieces[[3]]$worse,
    amp_plus_ald_better = pieces[[3]]$better,
    amp_plus_ald_change = pieces[[3]]$significant - pieces[[1]]$significant,
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows() %>%
  arrange(desc(amp_excluding_ald_significant), cohort)

comparison_wide <- all_results %>%
  select(cohort, gene, comparison, FDR, logrank_direction, significant) %>%
  tidyr::pivot_wider(
    names_from = comparison,
    values_from = c(FDR, logrank_direction, significant),
    names_glue = "{comparison}_{.value}"
  )

gene_level_comparison <- comparison_wide %>%
  filter(
    coalesce(AMP_excluding_ALD_significant, FALSE) |
      coalesce(AMP_vs_all_diploids_significant, FALSE) |
      coalesce(AMP_plus_ALD_vs_other_diploids_significant, FALSE)
  ) %>%
  mutate(
    significance_pattern = paste0(
      ifelse(coalesce(AMP_excluding_ALD_significant, FALSE), "AMP_excluding_ALD; ", ""),
      ifelse(coalesce(AMP_vs_all_diploids_significant, FALSE), "AMP_all_diploids; ", ""),
      ifelse(coalesce(AMP_plus_ALD_vs_other_diploids_significant, FALSE), "AMP_plus_ALD; ", "")
    ),
    significance_pattern = sub("; $", "", significance_pattern)
  ) %>%
  transmute(
    cohort,
    gene,
    amp_excluding_ald_FDR = AMP_excluding_ALD_FDR,
    amp_excluding_ald_direction = AMP_excluding_ALD_logrank_direction,
    amp_all_diploids_FDR = AMP_vs_all_diploids_FDR,
    amp_all_diploids_direction = AMP_vs_all_diploids_logrank_direction,
    amp_plus_ald_FDR = AMP_plus_ALD_vs_other_diploids_FDR,
    amp_plus_ald_direction = AMP_plus_ALD_vs_other_diploids_logrank_direction,
    significance_pattern
  ) %>%
  arrange(cohort, pmin(amp_excluding_ald_FDR, amp_all_diploids_FDR, amp_plus_ald_FDR, na.rm = TRUE), gene)

readme <- data.frame(
  section = c(
    "Table overview",
    "Analysis settings",
    "Analysis settings",
    "Analysis settings",
    "Analysis settings",
    "Analysis settings",
    "Analysis settings",
    "Column notes",
    "Sheet contents"
  ),
  item = c(
    "Purpose",
    "Gene universe",
    "AMP excluding ALD",
    "AMP vs all diploids",
    "AMP + ALD",
    "Sample threshold",
    "Multiple testing",
    "Direction",
    "gene_level_comparison"
  ),
  description = c(
    "Sensitivity check of how the AMP reference-group definition affects log-rank survival results.",
    paste0("Default wide_peak amplified-gene lists at prevalence >= ", prev_thr, "."),
    "Default AMP comparison: true AMP vs diploid samples after excluding ALD from the reference group.",
    "New comparison: true AMP vs all diploid samples, including ALD in the reference group.",
    "New comparison: the union of true AMP and ALD vs remaining diploid samples that are not ALD.",
    "At least 10 samples in the tested group and 10 in the reference group after matching to OS data.",
    "Benjamini-Hochberg correction is applied separately within each cohort and comparison; significance is FDR < 0.05.",
    "Worse means observed events exceeded expected events in the tested group; better means fewer events than expected.",
    "Contains genes significant in at least one of the three comparisons; blank AMP-excluding-ALD direction cells indicate nonsignificant default AMP results for which direction was not stored in the original output."
  ),
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

workbook <- openxlsx::createWorkbook()
openxlsx::addWorksheet(workbook, "README")
openxlsx::addWorksheet(workbook, "Settings")
openxlsx::addWorksheet(workbook, "cross_cancer_summary")
openxlsx::addWorksheet(workbook, "gene_level_comparison")

body_wrap <- openxlsx::createStyle(valign = "top", wrapText = TRUE)
fdr_style <- openxlsx::createStyle(numFmt = "0.0000")

openxlsx::writeData(workbook, "README", readme, withFilter = FALSE)
run_settings <- c(expected_table9a_settings, list(cohorts = paste(cohorts, collapse = ", ")))
settings_df <- data.frame(
  setting = names(run_settings),
  value = vapply(run_settings, as.character, character(1)),
  stringsAsFactors = FALSE
)
openxlsx::writeData(workbook, "Settings", settings_df, withFilter = FALSE)
openxlsx::setColWidths(workbook, "Settings", cols = 1:2, widths = c(38, 90))
openxlsx::freezePane(workbook, "Settings", firstRow = TRUE)
openxlsx::writeData(workbook, "cross_cancer_summary", cross_cancer_summary, withFilter = TRUE)
openxlsx::writeData(workbook, "gene_level_comparison", gene_level_comparison, withFilter = TRUE)

sheet_columns <- c(
  README = ncol(readme),
  cross_cancer_summary = ncol(cross_cancer_summary),
  gene_level_comparison = ncol(gene_level_comparison)
)
for (sheet in names(sheet_columns)) {
  n_cols <- sheet_columns[[sheet]]
  openxlsx::freezePane(workbook, sheet, firstRow = TRUE)
}

openxlsx::addStyle(
  workbook, "README", body_wrap,
  rows = 2:(nrow(readme) + 1), cols = 1:3, gridExpand = TRUE
)
openxlsx::setColWidths(workbook, "README", cols = 1, widths = 22)
openxlsx::setColWidths(workbook, "README", cols = 2, widths = 28)
openxlsx::setColWidths(workbook, "README", cols = 3, widths = 105)
openxlsx::setColWidths(workbook, "cross_cancer_summary", cols = 1, widths = 12)
openxlsx::setColWidths(
  workbook, "cross_cancer_summary",
  cols = 2:ncol(cross_cancer_summary), widths = 16
)
openxlsx::setColWidths(workbook, "gene_level_comparison", cols = 1:2, widths = c(12, 18))
openxlsx::setColWidths(
  workbook, "gene_level_comparison",
  cols = 3:ncol(gene_level_comparison), widths = 20
)
openxlsx::addStyle(
  workbook, "gene_level_comparison", fdr_style,
  rows = 2:(nrow(gene_level_comparison) + 1), cols = c(3, 5, 7), gridExpand = TRUE
)

openxlsx::saveWorkbook(workbook, output_file, overwrite = TRUE)

message("[AMP reference sensitivity] Saved: ", output_file)
message("[AMP reference sensitivity] Gene-level rows: ", nrow(gene_level_comparison))
