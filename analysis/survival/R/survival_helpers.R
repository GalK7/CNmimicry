# Purpose:
#   Prepare survival endpoints and sample groups, and provide binary Cox helpers.
#
# Inputs:
#   - config.R settings, core expression/CNV data, ALD/DLD calls, and gene lists.
#   - TCGA clinical data for overall survival.
#
# Outputs:
#   - Matched survival inputs and Cox result tables returned to callers.
#   - Clinical downloads under analysis/survival/GDCdata and processed endpoints
#     under analysis/survival/cache; callers write analysis results.

source(file.path("R", "prepare_data.R"))

pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)][1]
  if (is.na(hit)) return(NA_character_)
  hit
}

get_survival_ald_tpm_settings <- function() {
  enabled <- !is.null(ald_min_amplified_median_tpm)
  if (enabled &&
      (!is.numeric(ald_min_amplified_median_tpm) ||
       length(ald_min_amplified_median_tpm) != 1L ||
       !is.finite(ald_min_amplified_median_tpm) ||
       ald_min_amplified_median_tpm < 0)) {
    stop("ald_min_amplified_median_tpm must be NULL or one finite nonnegative number.")
  }
  tag <- if (enabled) {
    paste0("aldTPMgt", format(ald_min_amplified_median_tpm, trim = TRUE, scientific = FALSE))
  } else {
    "aldTPMoff"
  }
  list(enabled = enabled, tag = tag, suffix = if (enabled) paste0("_", tag) else "")
}

resolve_survival_context <- function(cohort, prev_thr = 0.03) {
  sig_tag <- if (isTRUE(use_gistic_significant)) {
    paste0("with_SIG_", gistic_sig_gene_source)
  } else {
    "no_SIG"
  }
  gthr_tag <- if (isTRUE(use_gistic_thresholded_intersection)) "with_GISTICthr" else "no_GISTICthr"
  cnv_tag  <- if (isTRUE(use_cnv_median_normalization)) "cnvnorm" else "rawcnv"
  
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
  
  ald_tpm <- get_survival_ald_tpm_settings()
  if (ald_tpm$enabled) results_base <- file.path(results_base, ald_tpm$tag)
  res_dir <- file.path(results_base, cohort, paste0("prev_", format(prev_thr, nsmall = 2)))
  
  list(
    cohort = cohort,
    prev_thr = prev_thr,
    sig_tag = sig_tag,
    gthr_tag = gthr_tag,
    cnv_tag = cnv_tag,
    ald_tpm_filter_enabled = ald_tpm$enabled,
    ald_tpm_tag = ald_tpm$tag,
    ald_tpm_suffix = ald_tpm$suffix,
    res_dir = res_dir
  )
}

resolve_survival_clinical_paths <- function(cohort) {
  cache_dir_surv <- file.path(workflow_root, "analysis", "survival", "cache")
  gdc_dir_surv <- file.path(workflow_root, "analysis", "survival", "GDCdata")
  
  list(
    cache_dir = cache_dir_surv,
    gdc_dir = gdc_dir_surv,
    cache_file = file.path(cache_dir_surv, paste0(cohort, "_survival_clinical_OS.rds"))
  )
}

build_survival_clinical_os <- function(cohort, cdef, force_rebuild = FALSE) {
  paths <- resolve_survival_clinical_paths(cohort)
  dir.create(paths$cache_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$gdc_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (file.exists(paths$cache_file) && !isTRUE(force_rebuild)) {
    message("[survival] Using cached OS clinical data: ", paths$cache_file)
    cached <- readRDS(paths$cache_file)
    
    if (is.list(cached) && "surv_df" %in% names(cached)) {
      return(cached)
    }
    
    stop("Survival clinical cache has an unexpected format: ", paths$cache_file)
  }
  
  message("[survival] Preparing OS clinical data for: ", cohort)
  
  query_clin <- TCGAbiolinks::GDCquery(
    project = cdef$gdc_project,
    data.category = "Clinical",
    data.type = "Clinical Supplement",
    data.format = "BCR Biotab"
  )
  
  TCGAbiolinks::GDCdownload(query_clin, directory = paths$gdc_dir)
  clin_list <- TCGAbiolinks::GDCprepare(query_clin, directory = paths$gdc_dir)
  
  patient_df <- clin_list[[paste0("clinical_patient_", tolower(cohort))]]
  if (is.null(patient_df)) {
    stop("Could not find clinical_patient table for ", cohort)
  }
  if (!"bcr_patient_barcode" %in% colnames(patient_df)) {
    stop("Clinical patient table is missing bcr_patient_barcode for ", cohort)
  }
  
  death_col <- pick_col(patient_df, c("death_days_to", "days_to_death"))
  last_contact_col <- pick_col(patient_df, c("last_contact_days_to", "days_to_last_followup", "days_to_last_follow_up"))
  vital_col <- pick_col(patient_df, c("vital_status"))
  
  if (is.na(death_col)) stop("Could not find a death-time column for ", cohort)
  if (is.na(last_contact_col)) stop("Could not find a last-contact column for ", cohort)
  if (is.na(vital_col)) stop("Could not find vital_status column for ", cohort)
  
  surv_df <- patient_df %>%
    dplyr::transmute(
      patient = bcr_patient_barcode,
      vital_status = .data[[vital_col]],
      death_time_raw = suppressWarnings(as.numeric(.data[[death_col]])),
      last_contact_time_raw = suppressWarnings(as.numeric(.data[[last_contact_col]]))
    ) %>%
    dplyr::mutate(
      OS_event = ifelse(vital_status == "Dead", 1, 0),
      OS_time = ifelse(OS_event == 1, death_time_raw, last_contact_time_raw)
    ) %>%
    dplyr::filter(!is.na(patient), !is.na(OS_time), !is.na(OS_event)) %>%
    dplyr::distinct(patient, .keep_all = TRUE) %>%
    dplyr::filter(OS_time >= 0)
  
  cached <- list(
    surv_df = surv_df,
    metadata = list(
      cohort = cohort,
      gdc_project = cdef$gdc_project,
      data_category = "Clinical",
      data_type = "Clinical Supplement",
      data_format = "BCR Biotab",
      source_gdc_dir = paths$gdc_dir,
      cache_file = paths$cache_file,
      death_time_column = death_col,
      last_contact_column = last_contact_col,
      vital_status_column = vital_col,
      n_patients = nrow(surv_df),
      n_events = sum(surv_df$OS_event == 1, na.rm = TRUE),
      created_at = as.character(Sys.time())
    )
  )
  
  saveRDS(cached, paths$cache_file)
  message("[survival] Saved OS clinical cache: ", paths$cache_file)
  
  cached
}

build_survival_core_inputs <- function(cohort, prev_thr = 0.03,
                                       force_rebuild_survival_cache = FALSE) {
  ctx <- resolve_survival_context(cohort, prev_thr)
  cdef <- cancers[[cohort]]
  if (is.null(cdef)) stop("Cohort is not defined in config.R: ", cohort)
  
  clinical_inputs <- build_survival_clinical_os(
    cohort = cohort,
    cdef = cdef,
    force_rebuild = force_rebuild_survival_cache
  )
  surv_df <- clinical_inputs$surv_df
  
  prepared <- prepare_data(
    cohort = cohort,
    gdc_project = cdef$gdc_project,
    gdc_dir = cdef$gdc_dir,
    cache_dir = cache_dir,
    er_pos_only = isTRUE(cdef$er_pos_only)
  )
  
  cnv_mat <- prepared$cnv_mat
  
  if (isTRUE(use_cnv_median_normalization)) {
    denom <- apply(cnv_mat, 2, median, na.rm = TRUE)
    denom[is.na(denom) | denom == 0] <- NA_real_
    cnv_use <- sweep(cnv_mat, 2, denom, "/")
    
    amp_cutoff <- cnv_amp_cutoff_norm
    del_cutoff <- cnv_del_cutoff_norm
    dip_center <- cnv_diploid_norm
    dip_tol <- cnv_diploid_tol_norm
  } else {
    cnv_use <- cnv_mat
    amp_cutoff <- cnv_amp_cutoff
    del_cutoff <- cnv_del_cutoff
    dip_center <- cnv_diploid
    dip_tol <- 0
  }
  
  is_diploid <- function(x) {
    if (dip_tol == 0) {
      x == dip_center
    } else {
      abs(x - dip_center) <= dip_tol
    }
  }
  
  list(
    ctx = ctx,
    surv_df = surv_df,
    survival_clinical_metadata = clinical_inputs$metadata,
    cnv_use = cnv_use,
    amp_cutoff = amp_cutoff,
    del_cutoff = del_cutoff,
    is_diploid = is_diploid
  )
}

attach_deleted_inputs <- function(core) {
  dld_calls_path <- file.path(core$ctx$res_dir, "dld_calls_long.rds")
  del_genes_path <- file.path(core$ctx$res_dir, "deleted_genes.rds")
  
  if (!file.exists(dld_calls_path)) stop("Missing: ", dld_calls_path)
  if (!file.exists(del_genes_path)) stop("Missing: ", del_genes_path)
  
  dld_calls_long <- readRDS(dld_calls_path)
  deleted_genes <- readRDS(del_genes_path)
  
  dld_calls_patient <- dld_calls_long %>%
    dplyr::transmute(
      patient = substr(Sample, 1, 12),
      Gene = Gene
    ) %>%
    dplyr::distinct()
  
  deleted_genes_to_test <- intersect(deleted_genes, rownames(core$cnv_use))
  
  c(
    core,
    list(
      dld_calls_patient = dld_calls_patient,
      deleted_genes_to_test = deleted_genes_to_test
    )
  )
}

attach_amplified_inputs <- function(core) {
  ald_calls_path <- file.path(core$ctx$res_dir, "ald_calls_long.rds")
  amp_genes_path <- file.path(core$ctx$res_dir, "amplified_genes.rds")
  ald_genes_path <- if (core$ctx$ald_tpm_filter_enabled) {
    file.path(core$ctx$res_dir, "ald_eligible_genes.rds")
  } else {
    amp_genes_path
  }
  
  if (!file.exists(ald_calls_path)) stop("Missing: ", ald_calls_path)
  if (!file.exists(amp_genes_path)) stop("Missing: ", amp_genes_path)
  if (!file.exists(ald_genes_path)) stop("Missing: ", ald_genes_path)
  
  ald_calls_long <- readRDS(ald_calls_path)
  amplified_genes <- readRDS(amp_genes_path)
  ald_eligible_genes <- readRDS(ald_genes_path)
  
  ald_calls_patient <- ald_calls_long %>%
    dplyr::transmute(
      patient = substr(Sample, 1, 12),
      Gene = Gene
    ) %>%
    dplyr::distinct()
  
  amplified_genes_to_test <- intersect(amplified_genes, rownames(core$cnv_use))
  
  c(
    core,
    list(
      ald_calls_patient = ald_calls_patient,
      ald_eligible_genes_to_test = intersect(ald_eligible_genes, rownames(core$cnv_use)),
      amplified_genes_to_test = amplified_genes_to_test
    )
  )
}

run_cox_binary <- function(pos_patients, ref_patients, surv_df, min_pos = 10, min_ref = 10) {
  if (length(pos_patients) < min_pos || length(ref_patients) < min_ref) return(NULL)
  
  dat <- dplyr::bind_rows(
    data.frame(patient = pos_patients, group = 1),
    data.frame(patient = ref_patients, group = 0)
  ) %>%
    dplyr::distinct(patient, .keep_all = TRUE) %>%
    dplyr::inner_join(surv_df, by = "patient")
  
  n_pos <- sum(dat$group == 1)
  n_ref <- sum(dat$group == 0)
  
  if (n_pos < min_pos || n_ref < min_ref) return(NULL)
  
  cox_warning <- NA_character_
  fit <- tryCatch(
    withCallingHandlers(
      survival::coxph(survival::Surv(OS_time, OS_event) ~ group, data = dat),
      warning = function(w) {
        cox_warning <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    return(data.frame(
      n_pos = n_pos,
      n_ref = n_ref,
      coef = NA_real_,
      HR = NA_real_,
      p_value = NA_real_,
      lower95 = NA_real_,
      upper95 = NA_real_,
      cox_status = "cox_failed",
      cox_message = conditionMessage(fit),
      stringsAsFactors = FALSE
    ))
  }
  
  s <- tryCatch(
    summary(fit),
    error = function(e) e
  )
  
  if (inherits(s, "error")) {
    return(data.frame(
      n_pos = n_pos,
      n_ref = n_ref,
      coef = NA_real_,
      HR = NA_real_,
      p_value = NA_real_,
      lower95 = NA_real_,
      upper95 = NA_real_,
      cox_status = "summary_failed",
      cox_message = conditionMessage(s),
      stringsAsFactors = FALSE
    ))
  }
  
  coef_value <- s$coefficients[1, "coef"]
  hr_value <- s$coefficients[1, "exp(coef)"]
  p_value <- s$coefficients[1, "Pr(>|z|)"]
  lower95 <- s$conf.int[1, "lower .95"]
  upper95 <- s$conf.int[1, "upper .95"]
  
  cox_status <- "ok"
  cox_message <- cox_warning
  
  if (!all(is.finite(c(coef_value, hr_value, p_value, lower95, upper95)))) {
    cox_status <- "non_finite_estimate"
    if (is.na(cox_message)) {
      cox_message <- "Cox model produced a non-finite coefficient, hazard ratio, confidence interval, or p-value."
    }
  }
  
  data.frame(
    n_pos = n_pos,
    n_ref = n_ref,
    coef = coef_value,
    HR = hr_value,
    p_value = p_value,
    lower95 = lower95,
    upper95 = upper95,
    cox_status = cox_status,
    cox_message = cox_message,
    stringsAsFactors = FALSE
  )
}

run_survival_screen_DLD <- function(inputs, min_dld = 10, min_ref = 10) {
  res_df <- lapply(inputs$deleted_genes_to_test, function(g) {
    dld_patients <- inputs$dld_calls_patient %>%
      dplyr::filter(Gene == g) %>%
      dplyr::pull(patient) %>%
      unique()
    
    if (!g %in% rownames(inputs$cnv_use)) return(NULL)
    
    gene_cnv <- inputs$cnv_use[g, , drop = TRUE]
    dip_patients <- names(gene_cnv)[inputs$is_diploid(gene_cnv)]
    dip_patients <- unique(substr(dip_patients, 1, 12))
    
    ref_patients <- setdiff(dip_patients, dld_patients)
    
    res <- run_cox_binary(dld_patients, ref_patients, inputs$surv_df,
                          min_pos = min_dld, min_ref = min_ref)
    if (is.null(res)) return(NULL)
    
    cbind(
      data.frame(
        gene = g,
        n_dld = res$n_pos,
        n_other_diploid = res$n_ref,
        stringsAsFactors = FALSE
      ),
      res[, c("coef", "HR", "p_value", "lower95", "upper95", "cox_status", "cox_message")]
    )
  }) %>%
    dplyr::bind_rows()
  
  if (nrow(res_df) == 0) {
    return(res_df)
  }
  
  res_df %>%
    dplyr::mutate(FDR = p.adjust(p_value, method = "BH")) %>%
    dplyr::arrange(FDR, p_value)
}

run_survival_screen_DEL <- function(inputs, min_del = 10, min_ref = 10) {
  res_df <- lapply(inputs$deleted_genes_to_test, function(g) {
    if (!g %in% rownames(inputs$cnv_use)) return(NULL)
    
    gene_cnv <- inputs$cnv_use[g, , drop = TRUE]
    
    del_patients <- names(gene_cnv)[gene_cnv <= inputs$del_cutoff]
    del_patients <- unique(substr(del_patients, 1, 12))
    
    dip_patients <- names(gene_cnv)[inputs$is_diploid(gene_cnv)]
    dip_patients <- unique(substr(dip_patients, 1, 12))
    
    dld_patients <- inputs$dld_calls_patient %>%
      dplyr::filter(Gene == g) %>%
      dplyr::pull(patient) %>%
      unique()
    
    ref_patients <- setdiff(dip_patients, dld_patients)
    
    res <- run_cox_binary(del_patients, ref_patients, inputs$surv_df,
                          min_pos = min_del, min_ref = min_ref)
    if (is.null(res)) return(NULL)
    
    cbind(
      data.frame(
        gene = g,
        n_del = res$n_pos,
        n_other_diploid = res$n_ref,
        stringsAsFactors = FALSE
      ),
      res[, c("coef", "HR", "p_value", "lower95", "upper95", "cox_status", "cox_message")]
    )
  }) %>%
    dplyr::bind_rows()
  
  if (nrow(res_df) == 0) {
    return(res_df)
  }
  
  res_df %>%
    dplyr::mutate(FDR = p.adjust(p_value, method = "BH")) %>%
    dplyr::arrange(FDR, p_value)
}

run_survival_screen_DEL_all_others <- function(inputs, min_del = 10, min_ref = 10) {
  res_df <- lapply(inputs$deleted_genes_to_test, function(g) {
    if (!g %in% rownames(inputs$cnv_use)) return(NULL)
    
    gene_cnv <- inputs$cnv_use[g, , drop = TRUE]
    evaluable_patients <- unique(substr(names(gene_cnv)[!is.na(gene_cnv)], 1, 12))
    
    del_patients <- names(gene_cnv)[gene_cnv <= inputs$del_cutoff]
    del_patients <- unique(substr(del_patients, 1, 12))
    
    ref_patients <- setdiff(evaluable_patients, del_patients)
    
    res <- run_cox_binary(del_patients, ref_patients, inputs$surv_df,
                          min_pos = min_del, min_ref = min_ref)
    if (is.null(res)) return(NULL)
    
    cbind(
      data.frame(
        gene = g,
        n_del = res$n_pos,
        n_other = res$n_ref,
        stringsAsFactors = FALSE
      ),
      res[, c("coef", "HR", "p_value", "lower95", "upper95", "cox_status", "cox_message")]
    )
  }) %>%
    dplyr::bind_rows()
  
  if (nrow(res_df) == 0) {
    return(res_df)
  }
  
  res_df %>%
    dplyr::mutate(FDR = p.adjust(p_value, method = "BH")) %>%
    dplyr::arrange(FDR, p_value)
}

run_survival_screen_ALD <- function(inputs, min_ald = 10, min_ref = 10) {
  res_df <- lapply(inputs$ald_eligible_genes_to_test, function(g) {
    ald_patients <- inputs$ald_calls_patient %>%
      dplyr::filter(Gene == g) %>%
      dplyr::pull(patient) %>%
      unique()
    
    if (!g %in% rownames(inputs$cnv_use)) return(NULL)
    
    gene_cnv <- inputs$cnv_use[g, , drop = TRUE]
    dip_patients <- names(gene_cnv)[inputs$is_diploid(gene_cnv)]
    dip_patients <- unique(substr(dip_patients, 1, 12))
    
    ref_patients <- setdiff(dip_patients, ald_patients)
    
    res <- run_cox_binary(ald_patients, ref_patients, inputs$surv_df,
                          min_pos = min_ald, min_ref = min_ref)
    if (is.null(res)) return(NULL)
    
    cbind(
      data.frame(
        gene = g,
        n_ald = res$n_pos,
        n_other_diploid = res$n_ref,
        stringsAsFactors = FALSE
      ),
      res[, c("coef", "HR", "p_value", "lower95", "upper95", "cox_status", "cox_message")]
    )
  }) %>%
    dplyr::bind_rows()
  
  if (nrow(res_df) == 0) {
    return(res_df)
  }
  
  res_df %>%
    dplyr::mutate(FDR = p.adjust(p_value, method = "BH")) %>%
    dplyr::arrange(FDR, p_value)
}

run_survival_screen_AMP <- function(inputs, min_amp = 10, min_ref = 10) {
  res_df <- lapply(inputs$amplified_genes_to_test, function(g) {
    if (!g %in% rownames(inputs$cnv_use)) return(NULL)
    
    gene_cnv <- inputs$cnv_use[g, , drop = TRUE]
    
    amp_patients <- names(gene_cnv)[gene_cnv >= inputs$amp_cutoff]
    amp_patients <- unique(substr(amp_patients, 1, 12))
    
    dip_patients <- names(gene_cnv)[inputs$is_diploid(gene_cnv)]
    dip_patients <- unique(substr(dip_patients, 1, 12))
    
    ald_patients <- inputs$ald_calls_patient %>%
      dplyr::filter(Gene == g) %>%
      dplyr::pull(patient) %>%
      unique()
    
    ref_patients <- setdiff(dip_patients, ald_patients)
    
    res <- run_cox_binary(amp_patients, ref_patients, inputs$surv_df,
                          min_pos = min_amp, min_ref = min_ref)
    if (is.null(res)) return(NULL)
    
    cbind(
      data.frame(
        gene = g,
        n_amp = res$n_pos,
        n_other_diploid = res$n_ref,
        stringsAsFactors = FALSE
      ),
      res[, c("coef", "HR", "p_value", "lower95", "upper95", "cox_status", "cox_message")]
    )
  }) %>%
    dplyr::bind_rows()
  
  if (nrow(res_df) == 0) {
    return(res_df)
  }
  
  res_df %>%
    dplyr::mutate(FDR = p.adjust(p_value, method = "BH")) %>%
    dplyr::arrange(FDR, p_value)
}

run_survival_screen_AMP_all_others <- function(inputs, min_amp = 10, min_ref = 10) {
  res_df <- lapply(inputs$amplified_genes_to_test, function(g) {
    if (!g %in% rownames(inputs$cnv_use)) return(NULL)
    
    gene_cnv <- inputs$cnv_use[g, , drop = TRUE]
    evaluable_patients <- unique(substr(names(gene_cnv)[!is.na(gene_cnv)], 1, 12))
    
    amp_patients <- names(gene_cnv)[gene_cnv >= inputs$amp_cutoff]
    amp_patients <- unique(substr(amp_patients, 1, 12))
    
    ref_patients <- setdiff(evaluable_patients, amp_patients)
    
    res <- run_cox_binary(amp_patients, ref_patients, inputs$surv_df,
                          min_pos = min_amp, min_ref = min_ref)
    if (is.null(res)) return(NULL)
    
    cbind(
      data.frame(
        gene = g,
        n_amp = res$n_pos,
        n_other = res$n_ref,
        stringsAsFactors = FALSE
      ),
      res[, c("coef", "HR", "p_value", "lower95", "upper95", "cox_status", "cox_message")]
    )
  }) %>%
    dplyr::bind_rows()
  
  if (nrow(res_df) == 0) {
    return(res_df)
  }
  
  res_df %>%
    dplyr::mutate(FDR = p.adjust(p_value, method = "BH")) %>%
    dplyr::arrange(FDR, p_value)
}
