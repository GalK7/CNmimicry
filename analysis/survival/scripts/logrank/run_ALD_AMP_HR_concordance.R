# Purpose: Evaluate true-AMP effect-size concordance for significant ALD cohort-gene pairs.
# Inputs: Supplementary Table 12 with original ALD directions, core sample groups, and clinical OS.
# Outputs: Workbook under analysis/survival/outputs/logrank, in the configured
#   aldTPMgt<value>/ALD_AMP_HR_concordance subfolder when TPM filtering is enabled,
#   or ALD_AMP_HR_concordance when disabled.

source("config.R")
source(file.path("analysis", "survival", "R", "survival_helpers.R"))
source(file.path("analysis", "survival", "R", "logrank_helpers.R"))

stopifnot(requireNamespace("dplyr", quietly = TRUE))
stopifnot(requireNamespace("openxlsx", quietly = TRUE))
stopifnot(requireNamespace("survival", quietly = TRUE))

library(dplyr)
library(survival)

# ------------------------------------------------------------------------------
# ALD-to-AMP survival effect-size concordance
# ------------------------------------------------------------------------------
cohorts <- names(cancers)
prev_thr <- prevalence_thresholds[1]
ald_tpm <- get_survival_ald_tpm_settings()
min_pos <- 10
min_ref <- 10
fdr_cutoff <- 0.05
hr_thresholds <- c(1.1, 1.2, 1.5, 2.0)

if (!isTRUE(use_gistic_significant) || !identical(gistic_sig_gene_source, "wide_peak")) {
  stop("This analysis must be run with the default wide_peak GISTIC setting.")
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

logrank_output_root <- resolve_logrank_output_root()
input_file <- file.path(
  workflow_root,
  "plots",
  "Figure 7 - ALD DLD survival associations across cancers",
  paste0(
    "Supplementary_Table_12_logrank_",
    gthr_tag, "_", sig_tag, "_", cnv_tag, "_prev", prev_tag,
    "_min", min_pos, "_", min_ref, "_fdr", format(fdr_cutoff, nsmall = 2), ald_tpm$suffix, ".xlsx"
  )
)

output_dir <- file.path(logrank_output_root, "ALD_AMP_HR_concordance")
output_tag <- paste(sig_tag, gthr_tag, cnv_tag, paste0("prev", prev_tag), sep = "_")
output_tag <- paste0(output_tag, ald_tpm$suffix)
output_file <- file.path(output_dir, paste0("ALD_AMP_HR_", output_tag, ".xlsx"))

if (!file.exists(input_file)) {
  stop("Missing matching significant log-rank input file: ", input_file,
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

if (!"Settings" %in% openxlsx::getSheetNames(input_file)) {
  stop("Table 12 has no Settings sheet. Regenerate it with the Figure 7 builder.",
       call. = FALSE)
}
table9a_settings <- openxlsx::read.xlsx(input_file, sheet = "Settings")
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

safe_ratio <- function(numerator, denominator) {
  if (!is.finite(numerator) || !is.finite(denominator) || denominator == 0) {
    return(NA_real_)
  }
  numerator / denominator
}

safe_exp <- function(x) {
  if (!is.finite(x)) return(NA_real_)
  exp(x)
}

classify_hr <- function(ald_direction, amp_hr, threshold) {
  if (is.na(amp_hr) || !is.finite(amp_hr)) return("not_testable")
  if (is.na(ald_direction) || ald_direction == "not_assigned") return("not_assigned")

  lower_threshold <- 1 / threshold

  if (ald_direction == "worse") {
    if (amp_hr >= threshold) return("same_direction")
    if (amp_hr <= lower_threshold) return("opposite_direction")
    return("weak_effect")
  }

  if (ald_direction == "better") {
    if (amp_hr <= lower_threshold) return("same_direction")
    if (amp_hr >= threshold) return("opposite_direction")
    return("weak_effect")
  }

  "not_assigned"
}

run_amp_logrank_with_hr <- function(inputs, gene, min_pos, min_ref) {
  if (!gene %in% rownames(inputs$cnv_use)) {
    return(data.frame(
      gene = gene,
      amp_testable = FALSE,
      amp_status = "gene_missing_from_cnv",
      amp_message = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  gene_cnv <- inputs$cnv_use[gene, , drop = TRUE]

  amp_patients <- names(gene_cnv)[gene_cnv >= inputs$amp_cutoff]
  amp_patients <- unique(substr(amp_patients, 1, 12))

  diploid_patients <- names(gene_cnv)[inputs$is_diploid(gene_cnv)]
  diploid_patients <- unique(substr(diploid_patients, 1, 12))

  ald_patients <- inputs$ald_calls_patient %>%
    filter(.data$Gene == .env$gene) %>%
    pull(.data$patient) %>%
    unique()

  ref_patients <- setdiff(diploid_patients, ald_patients)

  if (length(amp_patients) < min_pos || length(ref_patients) < min_ref) {
    return(data.frame(
      gene = gene,
      amp_testable = FALSE,
      amp_status = "failed_pre_survival_sample_threshold",
      amp_message = NA_character_,
      n_amp = length(amp_patients),
      n_other_diploid = length(ref_patients),
      stringsAsFactors = FALSE
    ))
  }

  dat <- bind_rows(
    data.frame(patient = amp_patients, group = "group", stringsAsFactors = FALSE),
    data.frame(patient = ref_patients, group = "other_diploid", stringsAsFactors = FALSE)
  ) %>%
    distinct(.data$patient, .keep_all = TRUE) %>%
    inner_join(inputs$surv_df, by = "patient") %>%
    filter(!is.na(.data$OS_time), !is.na(.data$OS_event), !is.na(.data$group)) %>%
    mutate(group = factor(.data$group, levels = c("other_diploid", "group")))

  n_amp <- sum(dat$group == "group")
  n_ref <- sum(dat$group == "other_diploid")

  if (n_amp < min_pos || n_ref < min_ref) {
    return(data.frame(
      gene = gene,
      amp_testable = FALSE,
      amp_status = "failed_post_survival_sample_threshold",
      amp_message = NA_character_,
      n_amp = n_amp,
      n_other_diploid = n_ref,
      stringsAsFactors = FALSE
    ))
  }

  logrank_fit <- tryCatch(
    survival::survdiff(survival::Surv(OS_time, OS_event) ~ group, data = dat),
    error = function(e) e
  )

  if (inherits(logrank_fit, "error")) {
    return(data.frame(
      gene = gene,
      amp_testable = TRUE,
      amp_status = "logrank_failed",
      amp_message = conditionMessage(logrank_fit),
      n_amp = n_amp,
      n_other_diploid = n_ref,
      stringsAsFactors = FALSE
    ))
  }

  group_idx <- grep("group=group$", names(logrank_fit$obs))
  ref_idx <- grep("group=other_diploid$", names(logrank_fit$obs))
  if (length(group_idx) != 1) group_idx <- length(logrank_fit$obs)
  if (length(ref_idx) != 1) ref_idx <- setdiff(seq_along(logrank_fit$obs), group_idx)[1]

  observed_group <- unname(logrank_fit$obs[group_idx])
  expected_group <- unname(logrank_fit$exp[group_idx])
  observed_ref <- unname(logrank_fit$obs[ref_idx])
  expected_ref <- unname(logrank_fit$exp[ref_idx])

  oe_group <- safe_ratio(observed_group, expected_group)
  oe_ref <- safe_ratio(observed_ref, expected_ref)
  hr_observed_expected_ratio <- safe_ratio(oe_group, oe_ref)

  observed_minus_expected_group <- observed_group - expected_group
  variance_group <- NA_real_
  if (!is.null(logrank_fit$var)) {
    variance_group <- unname(logrank_fit$var[group_idx, group_idx])
  }
  pike_hr <- safe_exp(observed_minus_expected_group / variance_group)

  amp_direction <- ifelse(
    !is.finite(observed_minus_expected_group) ||
      abs(observed_minus_expected_group) < .Machine$double.eps,
    "not_assigned",
    ifelse(observed_minus_expected_group > 0, "worse", "better")
  )

  chisq <- unname(logrank_fit$chisq)
  p_value <- stats::pchisq(chisq, df = length(logrank_fit$n) - 1, lower.tail = FALSE)

  data.frame(
    gene = gene,
    amp_testable = TRUE,
    amp_status = "ok",
    amp_message = NA_character_,
    n_amp = n_amp,
    n_other_diploid = n_ref,
    n_events_amp = sum(dat$group == "group" & dat$OS_event == 1, na.rm = TRUE),
    n_events_other_diploid = sum(dat$group == "other_diploid" & dat$OS_event == 1, na.rm = TRUE),
    amp_observed = observed_group,
    amp_expected = expected_group,
    amp_observed_expected_ratio = oe_group,
    reference_observed = observed_ref,
    reference_expected = expected_ref,
    reference_observed_expected_ratio = oe_ref,
    amp_HR_observed_expected_ratio = hr_observed_expected_ratio,
    amp_variance = variance_group,
    amp_pike_HR = pike_hr,
    amp_observed_minus_expected = observed_minus_expected_group,
    amp_direction = amp_direction,
    amp_logrank_chisq = chisq,
    amp_logrank_p_value = p_value,
    stringsAsFactors = FALSE
  )
}

message("[ALD-AMP HR concordance] Reading Supplementary Table 12 ALD hits ...")
significant_hits <- openxlsx::read.xlsx(input_file, sheet = "data")
required_cols <- c(
  "cohort", "state", "gene", "FDR", "logrank_p_value",
  "logrank_chisq", "logrank_direction",
  "logrank_observed_group", "logrank_expected_group",
  "logrank_observed_minus_expected_group"
)
missing_cols <- setdiff(required_cols, colnames(significant_hits))
if (length(missing_cols) > 0) {
  stop("Input file is missing required column(s): ", paste(missing_cols, collapse = ", "))
}

ald_hits <- significant_hits %>%
  filter(.data$state == "ALD") %>%
  transmute(
    cohort = as.character(.data$cohort),
    gene = as.character(.data$gene),
    ALD_FDR = as.numeric(.data$FDR),
    ALD_logrank_p_value = as.numeric(.data$logrank_p_value),
    ALD_logrank_chisq = as.numeric(.data$logrank_chisq),
    ALD_logrank_observed_group = as.numeric(.data$logrank_observed_group),
    ALD_logrank_expected_group = as.numeric(.data$logrank_expected_group),
    ALD_logrank_observed_minus_expected_group = as.numeric(.data$logrank_observed_minus_expected_group),
    ALD_direction = as.character(.data$logrank_direction)
  ) %>%
  distinct(.data$cohort, .data$gene, .keep_all = TRUE) %>%
  arrange(.data$cohort, .data$ALD_FDR, .data$gene)

if (nrow(ald_hits) == 0) {
  stop("No significant ALD hits were found in: ", input_file)
}

message("[ALD-AMP HR concordance] Significant ALD cohort/gene pairs: ", nrow(ald_hits))

gene_level_results <- lapply(cohorts, function(cohort) {
  cohort_hits <- ald_hits %>% filter(.data$cohort == .env$cohort)
  if (nrow(cohort_hits) == 0) return(data.frame())

  message("[ALD-AMP HR concordance] Preparing ", cohort, " ...")
  core_inputs <- build_survival_core_inputs(
    cohort = cohort,
    prev_thr = prev_thr,
    force_rebuild_survival_cache = FALSE
  )
  inputs <- attach_amplified_inputs(core_inputs)

  amp_results <- lapply(cohort_hits$gene, function(gene) {
    run_amp_logrank_with_hr(
      inputs = inputs,
      gene = gene,
      min_pos = min_pos,
      min_ref = min_ref
    )
  }) %>%
    bind_rows() %>%
    mutate(cohort = cohort, .before = gene)

  cohort_hits %>%
    left_join(amp_results, by = c("cohort", "gene"))
}) %>%
  bind_rows()

gene_level_results$amp_FDR <- NA_real_
amp_tested <- is.finite(gene_level_results$amp_logrank_p_value)
gene_level_results$amp_FDR[amp_tested] <- stats::p.adjust(
  gene_level_results$amp_logrank_p_value[amp_tested], method = "BH"
)

for (threshold in hr_thresholds) {
  class_col <- paste0("HR_", gsub("\\.", "_", as.character(threshold)), "_class")
  gene_level_results[[class_col]] <- mapply(
    classify_hr,
    gene_level_results$ALD_direction,
    gene_level_results$amp_HR_observed_expected_ratio,
    MoreArgs = list(threshold = threshold),
    USE.NAMES = FALSE
  )
}

summarize_threshold <- function(df, threshold) {
  class_col <- paste0("HR_", gsub("\\.", "_", as.character(threshold)), "_class")
  data.frame(
    threshold = threshold,
    lower_threshold = 1 / threshold,
    input_ALD_hits = nrow(df),
    AMP_testable = sum(df$amp_testable, na.rm = TRUE),
    same_direction = sum(df[[class_col]] == "same_direction", na.rm = TRUE),
    opposite_direction = sum(df[[class_col]] == "opposite_direction", na.rm = TRUE),
    weak_effect = sum(df[[class_col]] == "weak_effect", na.rm = TRUE),
    not_testable = sum(df[[class_col]] == "not_testable", na.rm = TRUE),
    not_assigned = sum(df[[class_col]] == "not_assigned", na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

threshold_summary <- lapply(hr_thresholds, function(threshold) {
  summarize_threshold(gene_level_results, threshold)
}) %>%
  bind_rows()

cohort_summary <- lapply(cohorts, function(cohort) {
  df <- gene_level_results %>% filter(.data$cohort == .env$cohort)
  if (nrow(df) == 0) return(data.frame())
  lapply(hr_thresholds, function(threshold) {
    summarize_threshold(df, threshold) %>% mutate(cohort = cohort, .before = threshold)
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  arrange(.data$cohort, .data$threshold)

readme <- data.frame(
  section = c(
    "Table overview",
    "Input",
    "Comparison",
    "Multiple testing",
    "Sample threshold",
    "HR-like effect size",
    "Secondary HR estimate",
    "Direction matching",
    "Not assigned",
    "Thresholds",
    "Sheet contents",
    "Sheet contents",
    "Sheet contents"
  ),
  item = c(
    "Purpose",
    "ALD hits",
    "AMP follow-up",
    "amp_FDR",
    "10/10 rule",
    "Observed/expected HR",
    "Pike estimator",
    "Same/opposite direction",
    "Not assigned",
    "Threshold grid",
    "threshold_summary",
    "cohort_summary",
    "gene_level_results"
  ),
  description = c(
    "Follow-up testing whether significant ALD survival genes show concordant true-AMP survival effects by HR-like effect-size thresholds.",
    "Default wide_peak ALD log-rank hits at FDR < 0.05 from Figure 7 Supplementary Table 12.",
    "For each ALD hit, the same cohort/gene was tested as true AMP vs diploid reference excluding ALD, matching the default AMP comparison.",
    "amp_FDR is the Benjamini-Hochberg correction of all finite AMP log-rank p-values across cohorts in this follow-up. Untestable pairs retain missing FDR values. This correction does not change effect ratios or same/opposite/weak classifications and is not a test of concordance itself.",
    "At least 10 AMP samples and 10 reference samples are required after matching to OS data.",
    "Primary effect size: (observed_AMP / expected_AMP) divided by (observed_reference / expected_reference). Values >1 indicate worse survival in AMP; values <1 indicate better survival in AMP.",
    "Pike's estimator is included as a secondary column: exp((observed_AMP - expected_AMP) / variance_AMP). The threshold summaries use the observed/expected HR.",
    "Same/opposite direction compares the AMP observed/expected HR direction to the original ALD observed-minus-expected log-rank direction.",
    "Not assigned means the original ALD direction in Supplementary Table 12 could not be assigned or the AMP follow-up was not testable; this should be rare for the original significant ALD hits.",
    "Thresholds tested: 1.1, 1.2, 1.5, and 2.0, with reciprocal lower thresholds for better-survival effects.",
    "Overall counts across all ALD cohort/gene pairs for each HR threshold.",
    "Per-cohort counts for each HR threshold.",
    "One row per significant ALD cohort/gene pair, with ALD statistics, AMP follow-up statistics, HR estimates, and threshold classifications."
  ),
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "README")
openxlsx::addWorksheet(wb, "Settings")
openxlsx::addWorksheet(wb, "threshold_summary")
openxlsx::addWorksheet(wb, "cohort_summary")
openxlsx::addWorksheet(wb, "gene_level_results")

openxlsx::writeData(wb, "README", readme, withFilter = FALSE)
run_settings <- c(expected_table9a_settings, list(
  cohorts = paste(cohorts, collapse = ", "),
  hr_thresholds = paste(hr_thresholds, collapse = ", ")
))
settings_df <- data.frame(
  setting = names(run_settings),
  value = vapply(run_settings, as.character, character(1)),
  stringsAsFactors = FALSE
)
openxlsx::writeData(wb, "Settings", settings_df, withFilter = FALSE)
openxlsx::setColWidths(wb, "Settings", cols = 1:2, widths = c(38, 90))
openxlsx::freezePane(wb, "Settings", firstRow = TRUE)
openxlsx::writeData(wb, "threshold_summary", threshold_summary, withFilter = TRUE)
openxlsx::writeData(wb, "cohort_summary", cohort_summary, withFilter = TRUE)
openxlsx::writeData(wb, "gene_level_results", gene_level_results, withFilter = TRUE)

wrap_style <- openxlsx::createStyle(valign = "top", wrapText = TRUE)
num_style <- openxlsx::createStyle(numFmt = "0.0000")
scientific_style <- openxlsx::createStyle(numFmt = "0.00E+00")

for (sheet in c("README", "threshold_summary", "cohort_summary", "gene_level_results")) {
  openxlsx::freezePane(wb, sheet, firstRow = TRUE)
}

openxlsx::addStyle(
  wb, "README", wrap_style,
  rows = 2:(nrow(readme) + 1), cols = 1:ncol(readme), gridExpand = TRUE
)
openxlsx::setColWidths(wb, "README", cols = 1, widths = 22)
openxlsx::setColWidths(wb, "README", cols = 2, widths = 26)
openxlsx::setColWidths(wb, "README", cols = 3, widths = 110)

openxlsx::setColWidths(wb, "threshold_summary", cols = 1:ncol(threshold_summary), widths = 18)
openxlsx::setColWidths(wb, "cohort_summary", cols = 1:ncol(cohort_summary), widths = 18)
openxlsx::setColWidths(wb, "gene_level_results", cols = 1:ncol(gene_level_results), widths = 18)

numeric_cols <- which(vapply(gene_level_results, is.numeric, logical(1)))
if (length(numeric_cols) > 0) {
  openxlsx::addStyle(
    wb, "gene_level_results", num_style,
    rows = 2:(nrow(gene_level_results) + 1),
    cols = numeric_cols,
    gridExpand = TRUE
  )
}
p_cols <- which(colnames(gene_level_results) %in% c("ALD_logrank_p_value", "amp_logrank_p_value", "amp_FDR"))
if (length(p_cols) > 0) {
  openxlsx::addStyle(
    wb, "gene_level_results", scientific_style,
    rows = 2:(nrow(gene_level_results) + 1),
    cols = p_cols,
    gridExpand = TRUE
  )
}

openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)

message("[ALD-AMP HR concordance] Saved: ", output_file)
message("[ALD-AMP HR concordance] Gene-level rows: ", nrow(gene_level_results))
