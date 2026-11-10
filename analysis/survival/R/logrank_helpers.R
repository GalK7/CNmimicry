# Purpose:
#   Provide log-rank tests, Kaplan-Meier summaries, and result-reading helpers.
#
# Inputs:
#   - Matched sample groups and overall-survival data from survival_helpers.R.
#   - config.R settings and per-cohort log-rank result files for summary helpers.
#
# Outputs:
#   - Gene-level log-rank results, survival directions, and summary counts.
#   - Results are returned to callers, which write them under
#     analysis/survival/outputs/logrank.

source(file.path("analysis", "survival", "R", "survival_helpers.R"))

empty_logrank_df <- function(mode) {
  # Return an empty log-rank result table with the correct columns for this mode.
  n_group_col <- switch(
    mode,
    DLD = "n_dld",
    DEL = "n_del",
    ALD = "n_ald",
    AMP = "n_amp",
    stop("Unsupported mode: ", mode)
  )

  out <- data.frame(
    gene = character(),
    n_other_diploid = integer(),
    n_events_group = integer(),
    n_events_other_diploid = integer(),
    km_median_group = numeric(),
    km_median_other_diploid = numeric(),
    logrank_observed_group = numeric(),
    logrank_expected_group = numeric(),
    logrank_observed_minus_expected_group = numeric(),
    logrank_direction = character(),
    logrank_chisq = numeric(),
    logrank_p_value = numeric(),
    logrank_status = character(),
    logrank_message = character(),
    FDR = numeric(),
    stringsAsFactors = FALSE
  )

  out[[n_group_col]] <- integer()
  out <- out[, c(
    "gene",
    n_group_col,
    "n_other_diploid",
    "n_events_group",
    "n_events_other_diploid",
    "km_median_group",
    "km_median_other_diploid",
    "logrank_observed_group",
    "logrank_expected_group",
    "logrank_observed_minus_expected_group",
    "logrank_direction",
    "logrank_chisq",
    "logrank_p_value",
    "logrank_status",
    "logrank_message",
    "FDR"
  )]

  out
}

resolve_logrank_output_root <- function() {
  base_dir <- file.path(workflow_root, "analysis", "survival", "outputs", "logrank")
  
  if (isTRUE(use_gistic_significant) && identical(gistic_sig_gene_source, "genes_in_peak")) {
    base_dir <- file.path(base_dir, "genes_in_peak_version")
  }
  ald_tpm <- get_survival_ald_tpm_settings()
  if (ald_tpm$enabled) base_dir <- file.path(base_dir, ald_tpm$tag)
  
  base_dir
}

safe_km_median <- function(fit, group_level) {
  # Extract the Kaplan-Meier median survival time for one group from a survfit object.
  fit_summary <- tryCatch(
    summary(fit)$table,
    error = function(e) NULL
  )

  if (is.null(fit_summary)) return(NA_real_)

  if (is.null(dim(fit_summary))) {
    median_value <- fit_summary["median"]
    return(ifelse(is.na(median_value), NA_real_, as.numeric(median_value)))
  }

  row_name <- paste0("group=", group_level)
  if (!row_name %in% rownames(fit_summary)) return(NA_real_)

  median_value <- fit_summary[row_name, "median"]
  ifelse(is.na(median_value), NA_real_, as.numeric(median_value))
}

run_logrank_binary <- function(pos_patients, ref_patients, surv_df,
                               min_pos = 10, min_ref = 10) {
  # Run one two-group log-rank comparison for one gene.
  if (length(pos_patients) < min_pos || length(ref_patients) < min_ref) return(NULL)

  # Build a two-group survival table and match it to available OS data.
  dat <- dplyr::bind_rows(
    data.frame(patient = pos_patients, group = "group"),
    data.frame(patient = ref_patients, group = "other_diploid")
  ) %>%
    dplyr::distinct(patient, .keep_all = TRUE) %>%
    dplyr::inner_join(surv_df, by = "patient") %>%
    dplyr::filter(!is.na(OS_time), !is.na(OS_event), !is.na(group)) %>%
    dplyr::mutate(group = factor(group, levels = c("other_diploid", "group")))

  n_group <- sum(dat$group == "group")
  n_ref <- sum(dat$group == "other_diploid")

  # Re-apply the sample threshold after matching to available OS data.
  if (n_group < min_pos || n_ref < min_ref) return(NULL)

  # Count observed deaths in each group.
  n_events_group <- sum(dat$group == "group" & dat$OS_event == 1, na.rm = TRUE)
  n_events_ref <- sum(dat$group == "other_diploid" & dat$OS_event == 1, na.rm = TRUE)

  # Surv(OS_time, OS_event) defines the survival outcome:
  # OS_time is follow-up time, and OS_event marks death (1) or censoring (0).
  # survfit estimates the Kaplan-Meier survival curve separately for each group.
  km_fit <- tryCatch(
    survival::survfit(survival::Surv(OS_time, OS_event) ~ group, data = dat),
    error = function(e) e
  )

  km_median_group <- NA_real_
  km_median_ref <- NA_real_
  if (!inherits(km_fit, "error")) {
    km_median_group <- safe_km_median(km_fit, "group")
    km_median_ref <- safe_km_median(km_fit, "other_diploid")
  }

  # survdiff performs the log-rank test, asking whether the two KM curves differ.
  logrank_fit <- tryCatch(
    survival::survdiff(survival::Surv(OS_time, OS_event) ~ group, data = dat),
    error = function(e) e
  )

  if (inherits(logrank_fit, "error")) {
    return(data.frame(
      n_group = n_group,
      n_ref = n_ref,
      n_events_group = n_events_group,
      n_events_ref = n_events_ref,
      km_median_group = km_median_group,
      km_median_ref = km_median_ref,
      logrank_observed_group = NA_real_,
      logrank_expected_group = NA_real_,
      logrank_observed_minus_expected_group = NA_real_,
      logrank_direction = "not_assigned",
      logrank_chisq = NA_real_,
      logrank_p_value = NA_real_,
      logrank_status = "logrank_failed",
      logrank_message = conditionMessage(logrank_fit),
      stringsAsFactors = FALSE
    ))
  }

  chisq <- unname(logrank_fit$chisq)
  group_idx <- grep("group=group$", names(logrank_fit$obs))
  if (length(group_idx) != 1) group_idx <- length(logrank_fit$obs)

  observed_group <- unname(logrank_fit$obs[group_idx])
  expected_group <- unname(logrank_fit$exp[group_idx])
  observed_minus_expected_group <- observed_group - expected_group
  logrank_direction <- ifelse(
    !is.finite(observed_minus_expected_group) ||
      abs(observed_minus_expected_group) < .Machine$double.eps,
    "not_assigned",
    ifelse(observed_minus_expected_group > 0, "worse", "better")
  )

  # Convert the survdiff chi-square statistic to a log-rank p-value.
  p_value <- stats::pchisq(chisq, df = length(logrank_fit$n) - 1, lower.tail = FALSE)

  status <- "ok"
  message <- NA_character_
  if (!all(is.finite(c(chisq, p_value)))) {
    status <- "non_finite_estimate"
    message <- "Log-rank test produced a non-finite chi-square statistic or p-value."
  }

  data.frame(
    n_group = n_group,
    n_ref = n_ref,
    n_events_group = n_events_group,
    n_events_ref = n_events_ref,
    km_median_group = km_median_group,
    km_median_ref = km_median_ref,
    logrank_observed_group = observed_group,
    logrank_expected_group = expected_group,
    logrank_observed_minus_expected_group = observed_minus_expected_group,
    logrank_direction = logrank_direction,
    logrank_chisq = chisq,
    logrank_p_value = p_value,
    logrank_status = status,
    logrank_message = message,
    stringsAsFactors = FALSE
  )
}

run_logrank_screen <- function(mode, inputs_deleted, inputs_amplified,
                               min_pos = 10, min_ref = 10) {
  # Run log-rank comparisons across all testable genes for one cohort and one mode.
  if (mode %in% c("DLD", "DEL")) {
    inputs <- inputs_deleted
    genes_to_test <- inputs$deleted_genes_to_test
  } else if (mode %in% c("ALD", "AMP")) {
    inputs <- inputs_amplified
    genes_to_test <- if (mode == "ALD") inputs$ald_eligible_genes_to_test else
      inputs$amplified_genes_to_test
  } else {
    stop("Unsupported mode: ", mode)
  }

  res_df <- lapply(genes_to_test, function(g) {
    if (!g %in% rownames(inputs$cnv_use)) return(NULL)

    gene_cnv <- inputs$cnv_use[g, , drop = TRUE]

    if (mode == "DLD") {
      # Define DLD vs other diploid tumors for this gene.
      group_patients <- inputs$dld_calls_patient %>%
        dplyr::filter(Gene == g) %>%
        dplyr::pull(patient) %>%
        unique()

      diploid_patients <- names(gene_cnv)[inputs$is_diploid(gene_cnv)]
      ref_patients <- setdiff(unique(substr(diploid_patients, 1, 12)), group_patients)
      n_group_col <- "n_dld"
    }

    if (mode == "DEL") {
      # Define true deletion vs diploid tumors excluding DLD tumors for this gene.
      group_patients <- names(gene_cnv)[gene_cnv <= inputs$del_cutoff]
      group_patients <- unique(substr(group_patients, 1, 12))

      diploid_patients <- names(gene_cnv)[inputs$is_diploid(gene_cnv)]
      diploid_patients <- unique(substr(diploid_patients, 1, 12))

      dld_patients <- inputs$dld_calls_patient %>%
        dplyr::filter(Gene == g) %>%
        dplyr::pull(patient) %>%
        unique()

      ref_patients <- setdiff(diploid_patients, dld_patients)
      n_group_col <- "n_del"
    }

    if (mode == "ALD") {
      # Define ALD vs other diploid tumors for this gene.
      group_patients <- inputs$ald_calls_patient %>%
        dplyr::filter(Gene == g) %>%
        dplyr::pull(patient) %>%
        unique()

      diploid_patients <- names(gene_cnv)[inputs$is_diploid(gene_cnv)]
      ref_patients <- setdiff(unique(substr(diploid_patients, 1, 12)), group_patients)
      n_group_col <- "n_ald"
    }

    if (mode == "AMP") {
      # Define true amplification vs diploid tumors excluding ALD tumors for this gene.
      group_patients <- names(gene_cnv)[gene_cnv >= inputs$amp_cutoff]
      group_patients <- unique(substr(group_patients, 1, 12))

      diploid_patients <- names(gene_cnv)[inputs$is_diploid(gene_cnv)]
      diploid_patients <- unique(substr(diploid_patients, 1, 12))

      ald_patients <- inputs$ald_calls_patient %>%
        dplyr::filter(Gene == g) %>%
        dplyr::pull(patient) %>%
        unique()

      ref_patients <- setdiff(diploid_patients, ald_patients)
      n_group_col <- "n_amp"
    }

    # Test this gene; NULL means it failed the sample thresholds after OS matching.
    res <- run_logrank_binary(
      pos_patients = group_patients,
      ref_patients = ref_patients,
      surv_df = inputs$surv_df,
      min_pos = min_pos,
      min_ref = min_ref
    )
    if (is.null(res)) return(NULL)

    out <- data.frame(
      gene = g,
      n_other_diploid = res$n_ref,
      n_events_group = res$n_events_group,
      n_events_other_diploid = res$n_events_ref,
      km_median_group = res$km_median_group,
      km_median_other_diploid = res$km_median_ref,
      logrank_observed_group = res$logrank_observed_group,
      logrank_expected_group = res$logrank_expected_group,
      logrank_observed_minus_expected_group = res$logrank_observed_minus_expected_group,
      logrank_direction = res$logrank_direction,
      logrank_chisq = res$logrank_chisq,
      logrank_p_value = res$logrank_p_value,
      logrank_status = res$logrank_status,
      logrank_message = res$logrank_message,
      stringsAsFactors = FALSE
    )
    out[[n_group_col]] <- res$n_group
    out[, c(
      "gene",
      n_group_col,
      "n_other_diploid",
      "n_events_group",
      "n_events_other_diploid",
      "km_median_group",
      "km_median_other_diploid",
      "logrank_observed_group",
      "logrank_expected_group",
      "logrank_observed_minus_expected_group",
      "logrank_direction",
      "logrank_chisq",
      "logrank_p_value",
      "logrank_status",
      "logrank_message"
    )]
  }) %>%
    dplyr::bind_rows()

  if (nrow(res_df) == 0) return(res_df)

  # Apply BH/FDR correction across genes that were actually tested.
  res_df %>%
    dplyr::mutate(FDR = p.adjust(logrank_p_value, method = "BH")) %>%
    dplyr::arrange(FDR, logrank_p_value)
}

safe_read_logrank_csv <- function(path) {
  # Read an output CSV safely when building cross-cancer summaries.
  if (!file.exists(path)) return(data.frame())

  info <- file.info(path)
  if (is.na(info$size) || info$size == 0) return(data.frame())

  tryCatch(
    read.csv(path),
    error = function(e) data.frame()
  )
}

count_logrank_sig <- function(df, fdr_cutoff) {
  # Count genes passing a selected FDR cutoff.
  if (nrow(df) == 0 || !"FDR" %in% colnames(df)) return(0L)
  sum(df$FDR < fdr_cutoff, na.rm = TRUE)
}

count_logrank_non_ok <- function(df) {
  # Count tested genes whose log-rank model did not finish cleanly.
  if (nrow(df) == 0 || !"logrank_status" %in% colnames(df)) return(0L)
  sum(!is.na(df$logrank_status) & df$logrank_status != "ok")
}
