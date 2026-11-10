# Purpose: Fit unadjusted and IDH/grade-adjusted Cox models for significant LGG ALD/DLD log-rank hits.
# Inputs: LGG log-rank results, core sample groups, clinical OS, mutation-derived IDH status, and grade.
# Outputs: Cox results under analysis/survival/outputs/cox, with a genes_in_peak_version
#   subfolder for peak-restricted runs and an aldTPMgt<value> subfolder when TPM filtering
#   is enabled, followed by LGG/adjusted_followup_logrank_hits.

source("config.R")
source(file.path("analysis", "survival", "R", "survival_helpers.R"))
source(file.path("analysis", "survival", "R", "logrank_helpers.R"))

stopifnot(requireNamespace("dplyr", quietly = TRUE))
stopifnot(requireNamespace("survival", quietly = TRUE))

library(dplyr)
library(survival)

# ------------------------------------------------------------------------------
# LGG adjusted Cox follow-up for log-rank-selected survival hits
# ------------------------------------------------------------------------------
# The broad survival screen is unadjusted log-rank/Kaplan-Meier. This follow-up
# asks whether LGG ALD/DLD hits from that screen remain associated with survival
# after Cox adjustment for IDH status, tumor grade, or both.

cohort <- "LGG"
prev_thr <- prevalence_thresholds[1]
ald_tpm <- get_survival_ald_tpm_settings()
fdr_cutoff <- 0.05
modes <- c("ALD", "DLD")
adjustment_models <- list(
  unadjusted = character(0),
  IDH_status = c("IDH_status"),
  tumor_grade = c("tumor_grade"),
  IDH_status_tumor_grade = c("IDH_status", "tumor_grade")
)

min_group <- 10
min_ref <- 10
min_events <- 10
sample_filter_tag <- paste0("minPos", min_group, "_minRef", min_ref)
force_rebuild_survival_cache <- FALSE

sig_tag <- if (isTRUE(use_gistic_significant)) {
  paste0("with_SIG_", gistic_sig_gene_source)
} else {
  "no_SIG"
}
gthr_tag <- if (isTRUE(use_gistic_thresholded_intersection)) "with_GISTICthr" else "no_GISTICthr"
cnv_tag <- if (isTRUE(use_cnv_median_normalization)) "cnvnorm" else "rawcnv"
prev_tag <- format(prev_thr, nsmall = 2)

logrank_output_root <- resolve_logrank_output_root()
base_cox_output_root <- file.path(workflow_root, "analysis", "survival", "outputs", "cox")
cox_output_root <- if (isTRUE(use_gistic_significant) && identical(gistic_sig_gene_source, "genes_in_peak")) {
  file.path(base_cox_output_root, "genes_in_peak_version")
} else {
  base_cox_output_root
}
if (ald_tpm$enabled) cox_output_root <- file.path(cox_output_root, ald_tpm$tag)
adjusted_output_dir <- file.path(
  cox_output_root,
  cohort,
  "adjusted_followup_logrank_hits",
  sample_filter_tag
)
dir.create(adjusted_output_dir, recursive = TRUE, showWarnings = FALSE)

adjusted_file_name <- function(adjustment_model = NULL) {
  model_tag <- if (is.null(adjustment_model)) {
    "combined"
  } else if (adjustment_model == "unadjusted") {
    "unadjusted"
  } else if (adjustment_model == "IDH_status") {
    "IDH"
  } else if (adjustment_model == "tumor_grade") {
    "grade"
  } else if (adjustment_model == "IDH_status_tumor_grade") {
    "IDH_grade"
  } else {
    adjustment_model
  }
  
  paste0(
    cohort, "_logrankHits_adjCox_", model_tag, "_",
    gthr_tag, "_", sig_tag, "_", cnv_tag,
    "_", sample_filter_tag,
    "_prev_", prev_tag,
    "_fdr_", format(fdr_cutoff, nsmall = 2),
    ald_tpm$suffix,
    ".csv"
  )
}

combined_output_file <- file.path(adjusted_output_dir, adjusted_file_name())

# ------------------------------------------------------------------------------
# Clinical covariates
# ------------------------------------------------------------------------------

as_missing_na <- function(x) {
  missing_values <- c(
    "",
    "[Not Available]",
    "[Not Applicable]",
    "[Unknown]",
    "[Not Evaluated]",
    "[Discrepancy]"
  )
  x[x %in% missing_values] <- NA
  x
}

find_lgg_patient_file <- function() {
  clinical_root <- file.path(
    workflow_root,
    "analysis", "survival",
    "GDCdata",
    "TCGA-LGG",
    "Clinical",
    "Clinical_Supplement"
  )
  
  hits <- list.files(
    clinical_root,
    pattern = "clinical_patient_lgg\\.txt$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(hits) == 0) {
    stop("Could not find the LGG clinical_patient Biotab file under: ", clinical_root)
  }
  
  hits[1]
}

load_lgg_clinical_covariates <- function() {
  patient_file <- find_lgg_patient_file()
  
  patient_df <- read.delim(
    patient_file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = ""
  )
  
  # BCR Biotab files store a second header row and a CDE row after the true header.
  patient_df <- patient_df[-c(1, 2), , drop = FALSE]
  
  required_cols <- c("bcr_patient_barcode", "tumor_grade")
  missing_cols <- setdiff(required_cols, colnames(patient_df))
  if (length(missing_cols) > 0) {
    stop("LGG clinical patient table is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  patient_df %>%
    transmute(
      patient = bcr_patient_barcode,
      tumor_grade = as_missing_na(tumor_grade)
    ) %>%
    filter(!is.na(patient)) %>%
    distinct(patient, .keep_all = TRUE) %>%
    mutate(tumor_grade = factor(tumor_grade))
}

load_lgg_idh_status <- function() {
  idh_file <- file.path(
    workflow_root,
    "analysis",
    "outputs",
    "idh",
    paste0(
      "LGG_IDH_ALD_per_sample_",
      gthr_tag, "_", sig_tag, "_", cnv_tag,
      "_prev_", prev_tag,
      ald_tpm$suffix,
      ".csv"
    )
  )
  
  if (!file.exists(idh_file)) {
    stop(
      "Missing mutation-derived LGG IDH status file: ", idh_file,
      "\nRun analysis/idh/run_IDH_LGG_ALD_DLD_analysis.R first."
    )
  }
  
  read.csv(idh_file, stringsAsFactors = FALSE) %>%
    transmute(
      patient = patient,
      IDH_status = IDH_status
    ) %>%
    distinct(patient, .keep_all = TRUE) %>%
    mutate(IDH_status = factor(IDH_status, levels = c("IDH-wt", "IDH-mut")))
}

# ------------------------------------------------------------------------------
# Log-rank hit selection and patient sets
# ------------------------------------------------------------------------------

read_logrank_hits <- function(mode) {
  if (mode == "ALD") {
    file_name <- paste0(
      cohort, "_ALD_vs_otherDiploids_OS_logrank_",
      gthr_tag, "_", sig_tag, "_", cnv_tag,
      "_", sample_filter_tag,
      "_prev_", prev_tag,
      ald_tpm$suffix,
      ".csv"
    )
  } else if (mode == "DLD") {
    file_name <- paste0(
      cohort, "_DLD_vs_otherDiploids_OS_logrank_",
      gthr_tag, "_", sig_tag, "_", cnv_tag,
      "_", sample_filter_tag,
      "_prev_", prev_tag,
      ald_tpm$suffix,
      ".csv"
    )
  } else {
    stop("Unsupported mode: ", mode)
  }
  
  path <- file.path(logrank_output_root, cohort, file_name)
  if (!file.exists(path)) stop("Missing log-rank survival file: ", path)
  
  read.csv(path, stringsAsFactors = FALSE) %>%
    filter(
      logrank_status == "ok",
      !is.na(FDR),
      FDR < fdr_cutoff
    ) %>%
    arrange(FDR, logrank_p_value)
}

get_group_and_reference <- function(gene, mode, inputs_deleted, inputs_amplified) {
  if (mode == "ALD") {
    inputs <- inputs_amplified
    
    group_patients <- inputs$ald_calls_patient %>%
      filter(Gene == gene) %>%
      pull(patient) %>%
      unique()
  } else if (mode == "DLD") {
    inputs <- inputs_deleted
    
    group_patients <- inputs$dld_calls_patient %>%
      filter(Gene == gene) %>%
      pull(patient) %>%
      unique()
  } else {
    stop("Unsupported mode: ", mode)
  }
  
  if (!gene %in% rownames(inputs$cnv_use)) {
    return(NULL)
  }
  
  gene_cnv <- inputs$cnv_use[gene, , drop = TRUE]
  diploid_patients <- names(gene_cnv)[inputs$is_diploid(gene_cnv)]
  diploid_patients <- unique(substr(diploid_patients, 1, 12))
  
  list(
    group_patients = group_patients,
    ref_patients = setdiff(diploid_patients, group_patients)
  )
}

# ------------------------------------------------------------------------------
# Adjusted Cox model
# ------------------------------------------------------------------------------

run_adjusted_cox <- function(gene, mode, group_patients, ref_patients,
                             surv_df, covariates_df, adjustment_model,
                             covariate_names) {
  if (length(group_patients) < min_group || length(ref_patients) < min_ref) {
    return(data.frame(
      adjustment_model = adjustment_model,
      n_model = NA_integer_,
      n_group = length(group_patients),
      n_ref = length(ref_patients),
      n_events = NA_integer_,
      adjusted_coef = NA_real_,
      adjusted_HR = NA_real_,
      adjusted_p_value = NA_real_,
      adjusted_lower95 = NA_real_,
      adjusted_upper95 = NA_real_,
      adjusted_cox_status = "too_few_patients_before_covariates",
      adjusted_cox_message = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  
  dat <- bind_rows(
    data.frame(patient = group_patients, group = 1),
    data.frame(patient = ref_patients, group = 0)
  ) %>%
    distinct(patient, .keep_all = TRUE) %>%
    inner_join(surv_df, by = "patient") %>%
    filter(
      !is.na(OS_time),
      !is.na(OS_event),
      !is.na(group)
    ) %>%
    mutate(
      group = factor(group, levels = c(0, 1))
    )

  if (length(covariate_names) > 0) {
    dat <- dat %>%
      inner_join(covariates_df, by = "patient") %>%
      mutate(
        tumor_grade = droplevels(tumor_grade),
        IDH_status = droplevels(IDH_status)
      )
  }
  
  for (covariate in covariate_names) {
    dat <- dat %>% filter(!is.na(.data[[covariate]]))
  }
  
  if (length(covariate_names) > 0) {
    dat <- dat %>%
      mutate(
        tumor_grade = droplevels(tumor_grade),
        IDH_status = droplevels(IDH_status)
      )
  }
  
  n_group <- sum(dat$group == 1)
  n_ref <- sum(dat$group == 0)
  n_events <- sum(dat$OS_event == 1, na.rm = TRUE)
  
  if (n_group < min_group || n_ref < min_ref || n_events < min_events) {
    return(data.frame(
      adjustment_model = adjustment_model,
      n_model = nrow(dat),
      n_group = n_group,
      n_ref = n_ref,
      n_events = n_events,
      adjusted_coef = NA_real_,
      adjusted_HR = NA_real_,
      adjusted_p_value = NA_real_,
      adjusted_lower95 = NA_real_,
      adjusted_upper95 = NA_real_,
      adjusted_cox_status = "too_few_patients_after_covariates",
      adjusted_cox_message = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  
  covariate_levels_ok <- all(vapply(
    covariate_names,
    function(covariate) nlevels(droplevels(dat[[covariate]])) >= 2,
    logical(1)
  ))
  
  if (!covariate_levels_ok) {
    return(data.frame(
      adjustment_model = adjustment_model,
      n_model = nrow(dat),
      n_group = n_group,
      n_ref = n_ref,
      n_events = n_events,
      adjusted_coef = NA_real_,
      adjusted_HR = NA_real_,
      adjusted_p_value = NA_real_,
      adjusted_lower95 = NA_real_,
      adjusted_upper95 = NA_real_,
      adjusted_cox_status = "insufficient_covariate_variation",
      adjusted_cox_message = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  
  cox_warning <- NA_character_
  model_formula <- as.formula(
    paste("Surv(OS_time, OS_event) ~", paste(c("group", covariate_names), collapse = " + "))
  )
  
  fit <- tryCatch(
    withCallingHandlers(
      coxph(model_formula, data = dat),
      warning = function(w) {
        cox_warning <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    return(data.frame(
      adjustment_model = adjustment_model,
      n_model = nrow(dat),
      n_group = n_group,
      n_ref = n_ref,
      n_events = n_events,
      adjusted_coef = NA_real_,
      adjusted_HR = NA_real_,
      adjusted_p_value = NA_real_,
      adjusted_lower95 = NA_real_,
      adjusted_upper95 = NA_real_,
      adjusted_cox_status = "cox_failed",
      adjusted_cox_message = conditionMessage(fit),
      stringsAsFactors = FALSE
    ))
  }
  
  s <- tryCatch(summary(fit), error = function(e) e)
  if (inherits(s, "error")) {
    return(data.frame(
      adjustment_model = adjustment_model,
      n_model = nrow(dat),
      n_group = n_group,
      n_ref = n_ref,
      n_events = n_events,
      adjusted_coef = NA_real_,
      adjusted_HR = NA_real_,
      adjusted_p_value = NA_real_,
      adjusted_lower95 = NA_real_,
      adjusted_upper95 = NA_real_,
      adjusted_cox_status = "summary_failed",
      adjusted_cox_message = conditionMessage(s),
      stringsAsFactors = FALSE
    ))
  }
  
  coef_table <- as.data.frame(s$coefficients)
  ci_table <- as.data.frame(s$conf.int)
  
  group_row <- grep("^group", rownames(coef_table), value = TRUE)[1]
  if (is.na(group_row)) {
    return(data.frame(
      adjustment_model = adjustment_model,
      n_model = nrow(dat),
      n_group = n_group,
      n_ref = n_ref,
      n_events = n_events,
      adjusted_coef = NA_real_,
      adjusted_HR = NA_real_,
      adjusted_p_value = NA_real_,
      adjusted_lower95 = NA_real_,
      adjusted_upper95 = NA_real_,
      adjusted_cox_status = "missing_group_coefficient",
      adjusted_cox_message = cox_warning,
      stringsAsFactors = FALSE
    ))
  }
  
  coef_value <- coef_table[group_row, "coef"]
  hr_value <- coef_table[group_row, "exp(coef)"]
  p_value <- coef_table[group_row, "Pr(>|z|)"]
  lower95 <- ci_table[group_row, "lower .95"]
  upper95 <- ci_table[group_row, "upper .95"]
  
  adjusted_status <- "ok"
  adjusted_message <- cox_warning
  
  if (!all(is.finite(c(coef_value, hr_value, p_value, lower95, upper95)))) {
    adjusted_status <- "non_finite_group_estimate"
    if (is.na(adjusted_message)) {
      adjusted_message <- "Adjusted Cox model produced a non-finite group estimate."
    }
  }
  
  data.frame(
    adjustment_model = adjustment_model,
    n_model = nrow(dat),
    n_group = n_group,
    n_ref = n_ref,
    n_events = n_events,
    adjusted_coef = coef_value,
    adjusted_HR = hr_value,
    adjusted_p_value = p_value,
    adjusted_lower95 = lower95,
    adjusted_upper95 = upper95,
    adjusted_cox_status = adjusted_status,
    adjusted_cox_message = adjusted_message,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# Run follow-up
# ------------------------------------------------------------------------------

message("[LGG log-rank follow-up] Loading LGG core survival/CNV inputs ...")
core_inputs <- build_survival_core_inputs(
  cohort = cohort,
  prev_thr = prev_thr,
  force_rebuild_survival_cache = force_rebuild_survival_cache
)
inputs_deleted <- attach_deleted_inputs(core_inputs)
inputs_amplified <- attach_amplified_inputs(core_inputs)

clinical_covariates <- load_lgg_clinical_covariates()
idh_status_df <- load_lgg_idh_status()
covariates_df <- clinical_covariates %>%
  inner_join(idh_status_df, by = "patient")

message("[LGG log-rank follow-up] Patients with tumor grade and IDH labels: ", nrow(covariates_df))

all_adjusted <- lapply(names(adjustment_models), function(adjustment_model) {
  covariate_names <- adjustment_models[[adjustment_model]]
  message(
    "[LGG log-rank follow-up] Running model ", adjustment_model,
    " with covariates: ", paste(covariate_names, collapse = ", ")
  )
  
  model_res <- lapply(modes, function(mode) {
    message("[LGG log-rank follow-up] Running adjusted Cox follow-up for ", mode, " log-rank hits ...")
    logrank_df <- read_logrank_hits(mode)
    message("[LGG log-rank follow-up] Significant log-rank ", mode, " genes: ", nrow(logrank_df))
    
    mode_res <- lapply(logrank_df$gene, function(gene) {
      patient_sets <- get_group_and_reference(
        gene = gene,
        mode = mode,
        inputs_deleted = inputs_deleted,
        inputs_amplified = inputs_amplified
      )
      if (is.null(patient_sets)) return(NULL)
      
      adjusted <- run_adjusted_cox(
        gene = gene,
        mode = mode,
        group_patients = patient_sets$group_patients,
        ref_patients = patient_sets$ref_patients,
        surv_df = core_inputs$surv_df,
        covariates_df = covariates_df,
        adjustment_model = adjustment_model,
        covariate_names = covariate_names
      )
      
      logrank_gene <- logrank_df %>%
        filter(gene == !!gene) %>%
        slice(1)
      
      data.frame(
        cohort = cohort,
        mode = mode,
        gene = gene,
        logrank_n_group = if (mode == "ALD") logrank_gene$n_ald else logrank_gene$n_dld,
        logrank_n_other_diploid = logrank_gene$n_other_diploid,
        logrank_n_events_group = logrank_gene$n_events_group,
        logrank_n_events_other_diploid = logrank_gene$n_events_other_diploid,
        logrank_km_median_group = logrank_gene$km_median_group,
        logrank_km_median_other_diploid = logrank_gene$km_median_other_diploid,
        logrank_chisq = logrank_gene$logrank_chisq,
        logrank_p_value = logrank_gene$logrank_p_value,
        logrank_FDR = logrank_gene$FDR,
        adjusted,
        stringsAsFactors = FALSE
      )
    }) %>%
      bind_rows()
    
    mode_res
  }) %>%
    bind_rows()
  
  model_res
}) %>%
  bind_rows()

if (nrow(all_adjusted) > 0) {
  all_adjusted <- all_adjusted %>%
    group_by(adjustment_model, mode) %>%
    mutate(adjusted_FDR = p.adjust(adjusted_p_value, method = "BH")) %>%
    ungroup() %>%
    arrange(adjustment_model, mode, adjusted_FDR, adjusted_p_value)
}

write.csv(all_adjusted, combined_output_file, row.names = FALSE)

for (adjustment_model in names(adjustment_models)) {
  model_output_file <- file.path(
    adjusted_output_dir,
    adjusted_file_name(adjustment_model)
  )
  dir.create(dirname(model_output_file), recursive = TRUE, showWarnings = FALSE)
  
  write.csv(
    all_adjusted %>% filter(adjustment_model == !!adjustment_model),
    model_output_file,
    row.names = FALSE
  )
  message("[LGG log-rank follow-up] Saved model-specific file: ", model_output_file)
}

message("[LGG log-rank follow-up] Saved: ", combined_output_file)
message("[LGG log-rank follow-up] Adjusted models written: ", nrow(all_adjusted))
if (nrow(all_adjusted) > 0) {
  adjusted_counts <- all_adjusted %>%
    group_by(adjustment_model, mode) %>%
    summarize(
      fitted_models = n(),
      adjusted_fdr_lt_0p05 = sum(adjusted_FDR < 0.05, na.rm = TRUE),
      .groups = "drop"
    )
  print(adjusted_counts)
}
