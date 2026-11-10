# Purpose:
#   Recalculate ALD/DLD within predefined TCGA subtype groups and compare
#   subtype-level ALD/DLD retention and survival results.
#
# Inputs:
#   - Prepared expression and copy-number data for the selected cohort.
#   - TCGA subtype annotations defined in this script.
#   - Whole-cohort survival results and cached overall-survival clinical data.
#   - config.R ALD amplified-median TPM eligibility setting, reapplied per subtype.
#
# Outputs:
#   - Subtype retention summaries, gene-level retention tables, and survival
#     comparison CSV files written to analysis/outputs/subtype.
#   - Enabled TPM filtering adds a separate aldTPMgt<value> output tag.

source("config.R")
source(file.path("R", "prepare_data.R"))
source(file.path("R", "get_gistic_significant_sets.R"))
source(file.path("analysis", "survival", "R", "survival_helpers.R"))
source(file.path("analysis", "survival", "R", "logrank_helpers.R"))

stopifnot(requireNamespace("TCGAbiolinks", quietly = TRUE))
stopifnot(requireNamespace("dplyr", quietly = TRUE))
stopifnot(requireNamespace("tidyr", quietly = TRUE))
stopifnot(requireNamespace("tibble", quietly = TRUE))
stopifnot(requireNamespace("survival", quietly = TRUE))

library(dplyr)
library(tidyr)
library(tibble)
library(survival)

# ------------------------------------------------------------------------------
# Subtype ALD/DLD pilot settings
# ------------------------------------------------------------------------------
pilot_cohort <- "UCEC"
prev_thr <- prevalence_thresholds[1]
ald_dld_percent_cutoff <- 20
min_pos <- 10
min_ref <- 10

pilot_configs <- list(
  LGG = list(
    subtype_column = "IDH.codel.subtype",
    subtype_map = list(
      `IDHmut-codel` = "IDHmut-codel",
      `IDHmut-non-codel` = "IDHmut-non-codel",
      IDHwt = "IDHwt"
    )
  ),
  SARC = list(
    subtype_column = "short histo",
    subtype_map = list(
      SARC_DDLPS = "DDLPS",
      SARC_LMS = c("STLMS", "ULMS"),
      SARC_MFS_UPS = c("MFS", "UPS"),
      SARC_Other = c("SS", "MPNST")
    )
  ),
  UCEC = list(
    subtype_column = "IntegrativeCluster",
    subtype_map = list(
      UCEC_CN_high = "CN high",
      UCEC_CN_low = "CN low",
      UCEC_MSI = "MSI",
      UCEC_POLE = "POLE"
    )
  )
)

if (!pilot_cohort %in% names(pilot_configs)) {
  stop("No subtype pilot config is defined for: ", pilot_cohort)
}

pilot_config <- pilot_configs[[pilot_cohort]]
sample_filter_tag <- paste0("minPos", min_pos, "_minRef", min_ref)
gthr_tag <- if (isTRUE(use_gistic_thresholded_intersection)) "with_GISTICthr" else "no_GISTICthr"
sig_tag <- if (isTRUE(use_gistic_significant)) {
  paste0("with_SIG_", gistic_sig_gene_source)
} else {
  "no_SIG"
}
cnv_tag <- if (isTRUE(use_cnv_median_normalization)) "cnvnorm" else "rawcnv"
prev_tag <- format(prev_thr, nsmall = 2)
ald_tpm <- get_survival_ald_tpm_settings()
whole_label <- paste0("Whole_", pilot_cohort)

output_dir <- file.path(
  workflow_root,
  "analysis",
  "outputs",
  "subtype",
  pilot_cohort
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_tag <- paste(sig_tag, gthr_tag, cnv_tag, paste0("prev", prev_tag), sep = "_")
output_tag <- paste0(output_tag, ald_tpm$suffix)
summary_file <- file.path(output_dir, paste0(pilot_cohort, "_retention_summary_", output_tag, ".csv"))
gene_file <- file.path(output_dir, paste0(pilot_cohort, "_retention_genes_", output_tag, ".csv"))
survival_file <- file.path(output_dir, paste0(pilot_cohort, "_survival_comparison_", output_tag, ".csv"))

clean_subtype <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "Notassigned", "Not assigned", "Not Assigned")] <- NA_character_
  x
}

as_logical_flag <- function(x) {
  if (is.logical(x)) return(x)
  ifelse(is.na(x), NA, tolower(as.character(x)) == "true")
}

load_pilot_subtypes <- function(cohort, config) {
  subtype_df <- TCGAbiolinks::TCGAquery_subtype(cohort)
  required_cols <- c("patient", config$subtype_column)
  missing_cols <- setdiff(required_cols, colnames(subtype_df))
  if (length(missing_cols) > 0) {
    stop("Subtype table is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  raw <- subtype_df %>%
    transmute(
      patient = patient,
      subtype_raw = clean_subtype(.data[[config$subtype_column]])
    ) %>%
    filter(!is.na(patient), !is.na(subtype_raw)) %>%
    distinct(patient, .keep_all = TRUE)

  mapped <- bind_rows(lapply(names(config$subtype_map), function(label) {
    values <- config$subtype_map[[label]]
    raw %>%
      filter(subtype_raw %in% values) %>%
      mutate(subtype = label)
  }))

  mapped %>%
    select(patient, subtype, subtype_raw) %>%
    distinct(patient, .keep_all = TRUE)
}

prepare_pilot_inputs <- function(cohort) {
  cdef <- cancers[[cohort]]
  if (is.null(cdef)) stop("Missing cohort in config.R: ", cohort)

  prepared <- prepare_data(
    cohort = cohort,
    gdc_project = cdef$gdc_project,
    gdc_dir = cdef$gdc_dir,
    cache_dir = cache_dir,
    er_pos_only = isTRUE(cdef$er_pos_only)
  )

  expr_mat <- prepared$expr_mat
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

  gistic_sig <- NULL
  if (isTRUE(use_gistic_significant)) {
    gistic_sig <- get_gistic_significant_sets(
      cohort = cohort,
      cache_dir = cache_dir,
      gdac_resources_dir = gdac_resources_dir,
      gdac_run_id = gdac_run_id,
      gistic_sig_q_cutoff = gistic_sig_q_cutoff,
      gistic_sig_gene_source = gistic_sig_gene_source
    )
  }

  list(
    expr_mat = expr_mat,
    cnv_use = cnv_use,
    amp_cutoff = amp_cutoff,
    del_cutoff = del_cutoff,
    is_diploid = is_diploid,
    gistic_sig = gistic_sig
  )
}

calc_ald_dld_for_samples <- function(label, sample_ids, inputs) {
  sample_ids <- intersect(sample_ids, colnames(inputs$expr_mat))
  sample_ids <- intersect(sample_ids, colnames(inputs$cnv_use))
  if (length(sample_ids) == 0) stop("No matching samples for label: ", label)

  expr_mat <- inputs$expr_mat[, sample_ids, drop = FALSE]
  cnv_use <- inputs$cnv_use[, sample_ids, drop = FALSE]

  binary_amp <- cnv_use >= inputs$amp_cutoff
  amp_pct <- rowSums(binary_amp, na.rm = TRUE) / ncol(binary_amp)
  amplified_genes <- names(amp_pct)[amp_pct > prev_thr]
  amplified_genes <- intersect(amplified_genes, rownames(expr_mat))

  if (isTRUE(use_gistic_significant) && !is.null(inputs$gistic_sig)) {
    amplified_genes <- intersect(amplified_genes, inputs$gistic_sig$sig_amp_genes)
  }

  binary_del <- cnv_use <= inputs$del_cutoff
  del_pct <- rowSums(binary_del, na.rm = TRUE) / ncol(binary_del)
  deleted_genes <- names(del_pct)[del_pct > prev_thr]
  deleted_genes <- intersect(deleted_genes, rownames(expr_mat))

  if (isTRUE(use_gistic_significant) && !is.null(inputs$gistic_sig)) {
    deleted_genes <- intersect(deleted_genes, inputs$gistic_sig$sig_del_genes)
  }

  if (length(amplified_genes) > 0) {
    amp_status_df <- as.data.frame(t(binary_amp[amplified_genes, , drop = FALSE])) %>%
      rownames_to_column("Sample") %>%
      pivot_longer(-Sample, names_to = "Gene", values_to = "Amplification_Status")

    expr_amp_df <- as.data.frame(t(expr_mat[amplified_genes, , drop = FALSE])) %>%
      rownames_to_column("Sample") %>%
      pivot_longer(-Sample, names_to = "Gene", values_to = "Expression")

    cnv_amp_df <- as.data.frame(t(cnv_use[amplified_genes, , drop = FALSE])) %>%
      rownames_to_column("Sample") %>%
      pivot_longer(-Sample, names_to = "Gene", values_to = "cnv_value")

    combined_amp <- expr_amp_df %>%
      inner_join(amp_status_df, by = c("Sample", "Gene")) %>%
      inner_join(cnv_amp_df, by = c("Sample", "Gene"))

    median_amp <- combined_amp %>%
      filter(Amplification_Status == TRUE) %>%
      group_by(Gene) %>%
      summarise(median_amplified = median(Expression, na.rm = TRUE), .groups = "drop")

    ald_gene_table <- combined_amp %>%
      left_join(median_amp, by = "Gene") %>%
      filter(
        if (ald_tpm$enabled) {
          is.finite(median_amplified) & median_amplified > ald_min_amplified_median_tpm
        } else TRUE
      ) %>%
      group_by(Gene) %>%
      summarise(
        outliers_count = sum(
          Amplification_Status == FALSE &
            inputs$is_diploid(cnv_value) &
            Expression > median_amplified,
          na.rm = TRUE
        ),
        total_diploid_reference = sum(Amplification_Status == FALSE & inputs$is_diploid(cnv_value), na.rm = TRUE),
        percentage = ifelse(total_diploid_reference > 0, 100 * outliers_count / total_diploid_reference, NA_real_),
        .groups = "drop"
      ) %>%
      mutate(analysis_label = label, state = "ALD", eligible_gene = TRUE, hit = percentage > ald_dld_percent_cutoff)
  } else {
    ald_gene_table <- tibble(
      analysis_label = character(),
      state = character(),
      Gene = character(),
      outliers_count = integer(),
      total_diploid_reference = integer(),
      percentage = numeric(),
      eligible_gene = logical(),
      hit = logical()
    )
  }

  if (length(deleted_genes) > 0) {
    del_status_df <- as.data.frame(t(binary_del[deleted_genes, , drop = FALSE])) %>%
      rownames_to_column("Sample") %>%
      pivot_longer(-Sample, names_to = "Gene", values_to = "Deletion_Status")

    expr_del_df <- as.data.frame(t(expr_mat[deleted_genes, , drop = FALSE])) %>%
      rownames_to_column("Sample") %>%
      pivot_longer(-Sample, names_to = "Gene", values_to = "Expression")

    cnv_del_df <- as.data.frame(t(cnv_use[deleted_genes, , drop = FALSE])) %>%
      rownames_to_column("Sample") %>%
      pivot_longer(-Sample, names_to = "Gene", values_to = "cnv_value")

    combined_del <- expr_del_df %>%
      inner_join(del_status_df, by = c("Sample", "Gene")) %>%
      inner_join(cnv_del_df, by = c("Sample", "Gene"))

    median_del <- combined_del %>%
      filter(Deletion_Status == TRUE) %>%
      group_by(Gene) %>%
      summarise(median_deleted = median(Expression, na.rm = TRUE), .groups = "drop")

    dld_gene_table <- combined_del %>%
      left_join(median_del, by = "Gene") %>%
      group_by(Gene) %>%
      summarise(
        outliers_count = sum(
          Deletion_Status == FALSE &
            inputs$is_diploid(cnv_value) &
            Expression < median_deleted,
          na.rm = TRUE
        ),
        total_diploid_reference = sum(Deletion_Status == FALSE & inputs$is_diploid(cnv_value), na.rm = TRUE),
        percentage = ifelse(total_diploid_reference > 0, 100 * outliers_count / total_diploid_reference, NA_real_),
        .groups = "drop"
      ) %>%
      mutate(analysis_label = label, state = "DLD", eligible_gene = TRUE, hit = percentage > ald_dld_percent_cutoff)
  } else {
    dld_gene_table <- tibble(
      analysis_label = character(),
      state = character(),
      Gene = character(),
      outliers_count = integer(),
      total_diploid_reference = integer(),
      percentage = numeric(),
      eligible_gene = logical(),
      hit = logical()
    )
  }

  gene_table <- bind_rows(ald_gene_table, dld_gene_table) %>%
    mutate(
      n_samples = length(sample_ids),
      prevalence_threshold = prev_thr,
      ald_dld_percent_cutoff = ald_dld_percent_cutoff
    ) %>%
    select(
      analysis_label,
      state,
      Gene,
      n_samples,
      eligible_gene,
      hit,
      outliers_count,
      total_diploid_reference,
      percentage,
      prevalence_threshold,
      ald_dld_percent_cutoff
    )

  summary_table <- gene_table %>%
    group_by(analysis_label, state) %>%
    summarise(
      n_samples = first(n_samples),
      eligible_genes = n(),
      hit_genes = sum(hit, na.rm = TRUE),
      cancer_score = ifelse(eligible_genes > 0, hit_genes / eligible_genes, NA_real_),
      .groups = "drop"
    )

  # Preserve an explicit ALD summary when the TPM rule excludes every gene.
  if (ald_tpm$enabled && !any(summary_table$state == "ALD")) {
    summary_table <- bind_rows(summary_table, tibble(
      analysis_label = label, state = "ALD", n_samples = length(sample_ids),
      eligible_genes = 0L, hit_genes = 0L, cancer_score = NA_real_
    ))
  }

  list(gene_table = gene_table, summary_table = summary_table)
}

get_group_reference <- function(gene, state, sample_ids, inputs) {
  sample_ids <- intersect(sample_ids, colnames(inputs$expr_mat))
  sample_ids <- intersect(sample_ids, colnames(inputs$cnv_use))
  if (!gene %in% rownames(inputs$expr_mat) || !gene %in% rownames(inputs$cnv_use)) return(NULL)
  if (length(sample_ids) == 0) return(NULL)

  expr <- inputs$expr_mat[gene, sample_ids, drop = TRUE]
  cnv <- inputs$cnv_use[gene, sample_ids, drop = TRUE]
  patients <- substr(sample_ids, 1, 12)

  if (state == "ALD") {
    altered <- cnv >= inputs$amp_cutoff
    altered_median <- median(expr[altered], na.rm = TRUE)
    if (!is.finite(altered_median)) return(NULL)
    if (ald_tpm$enabled && altered_median <= ald_min_amplified_median_tpm) return(NULL)

    group_samples <- !altered & inputs$is_diploid(cnv) & expr > altered_median
    ref_samples <- !altered & inputs$is_diploid(cnv) & !group_samples
  } else if (state == "DLD") {
    altered <- cnv <= inputs$del_cutoff
    altered_median <- median(expr[altered], na.rm = TRUE)
    if (!is.finite(altered_median)) return(NULL)

    group_samples <- !altered & inputs$is_diploid(cnv) & expr < altered_median
    ref_samples <- !altered & inputs$is_diploid(cnv) & !group_samples
  } else {
    stop("Unsupported state: ", state)
  }

  list(
    group_patients = unique(patients[group_samples]),
    ref_patients = unique(patients[ref_samples])
  )
}

run_subtype_logrank <- function(gene, state, subtype, sample_sets, inputs, surv_df) {
  groups <- get_group_reference(
    gene = gene,
    state = state,
    sample_ids = sample_sets[[subtype]],
    inputs = inputs
  )

  if (is.null(groups)) {
    return(data.frame(
      subtype_n_group = NA_integer_,
      subtype_n_reference = NA_integer_,
      subtype_n_events_group = NA_integer_,
      subtype_n_events_reference = NA_integer_,
      subtype_logrank_chisq = NA_real_,
      subtype_logrank_p_value = NA_real_,
      subtype_logrank_direction = NA_character_,
      subtype_logrank_status = "not_testable_no_group_definition",
      subtype_logrank_message = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  if (length(groups$group_patients) < min_pos || length(groups$ref_patients) < min_ref) {
    return(data.frame(
      subtype_n_group = length(groups$group_patients),
      subtype_n_reference = length(groups$ref_patients),
      subtype_n_events_group = NA_integer_,
      subtype_n_events_reference = NA_integer_,
      subtype_logrank_chisq = NA_real_,
      subtype_logrank_p_value = NA_real_,
      subtype_logrank_direction = NA_character_,
      subtype_logrank_status = "not_testable_low_samples",
      subtype_logrank_message = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  res <- run_logrank_binary(
    pos_patients = groups$group_patients,
    ref_patients = groups$ref_patients,
    surv_df = surv_df,
    min_pos = min_pos,
    min_ref = min_ref
  )

  if (is.null(res)) {
    return(data.frame(
      subtype_n_group = length(groups$group_patients),
      subtype_n_reference = length(groups$ref_patients),
      subtype_n_events_group = NA_integer_,
      subtype_n_events_reference = NA_integer_,
      subtype_logrank_chisq = NA_real_,
      subtype_logrank_p_value = NA_real_,
      subtype_logrank_direction = NA_character_,
      subtype_logrank_status = "not_testable_low_samples_after_survival_match",
      subtype_logrank_message = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    subtype_n_group = res$n_group,
    subtype_n_reference = res$n_ref,
    subtype_n_events_group = res$n_events_group,
    subtype_n_events_reference = res$n_events_ref,
    subtype_logrank_chisq = res$logrank_chisq,
    subtype_logrank_p_value = res$logrank_p_value,
    subtype_logrank_direction = res$logrank_direction,
    subtype_logrank_status = res$logrank_status,
    subtype_logrank_message = res$logrank_message,
    stringsAsFactors = FALSE
  )
}

load_whole_logrank <- function(cohort, state, inputs, surv_df, all_samples) {
  comparison_label <- switch(
    state,
    DLD = "DLD_vs_otherDiploids",
    ALD = "ALD_vs_otherDiploids",
    stop("Unsupported state: ", state)
  )
  state_group_col <- switch(state, DLD = "n_dld", ALD = "n_ald")

  whole_file <- file.path(
    resolve_logrank_output_root(),
    cohort,
    paste0(
      cohort, "_", comparison_label, "_OS_logrank_",
      gthr_tag, "_", sig_tag, "_", cnv_tag,
      "_", sample_filter_tag,
      "_prev_", prev_tag,
      ald_tpm$suffix,
      ".csv"
    )
  )

  if (!file.exists(whole_file)) {
    stop("Missing whole-cohort log-rank output: ", whole_file)
  }

  whole_df <- read.csv(whole_file, check.names = FALSE)

  whole_out <- whole_df %>%
    transmute(
      state = state,
      Gene = gene,
      whole_n_group = .data[[state_group_col]],
      whole_n_reference = n_other_diploid,
      whole_n_events_group = n_events_group,
      whole_n_events_reference = n_events_other_diploid,
      whole_logrank_chisq = logrank_chisq,
      whole_logrank_p_value = logrank_p_value,
      whole_logrank_FDR = FDR,
      whole_logrank_status = logrank_status
    )

  if ("logrank_direction" %in% colnames(whole_df)) {
    whole_out$whole_logrank_direction <- whole_df$logrank_direction
  } else {
    whole_out$whole_logrank_direction <- vapply(whole_out$Gene, function(g) {
      groups <- get_group_reference(
        gene = g,
        state = state,
        sample_ids = all_samples,
        inputs = inputs
      )
      if (is.null(groups)) return(NA_character_)

      res <- run_logrank_binary(
        pos_patients = groups$group_patients,
        ref_patients = groups$ref_patients,
        surv_df = surv_df,
        min_pos = min_pos,
        min_ref = min_ref
      )
      if (is.null(res)) return(NA_character_)
      res$logrank_direction
    }, character(1))
  }

  whole_out %>%
    select(
      state,
      Gene,
      whole_n_group,
      whole_n_reference,
      whole_n_events_group,
      whole_n_events_reference,
      whole_logrank_chisq,
      whole_logrank_p_value,
      whole_logrank_FDR,
      whole_logrank_direction,
      whole_logrank_status
    )
}

message("[subtype pilot] Cohort: ", pilot_cohort)
message("[subtype pilot] Loading subtype labels ...")
subtypes <- load_pilot_subtypes(pilot_cohort, pilot_config)

message("[subtype pilot] Loading expression/CNV inputs ...")
inputs <- prepare_pilot_inputs(pilot_cohort)
all_samples <- colnames(inputs$expr_mat)
all_patients <- substr(all_samples, 1, 12)

sample_sets <- list()
sample_sets[[whole_label]] <- all_samples

for (subtype_name in names(pilot_config$subtype_map)) {
  subtype_patients <- subtypes %>%
    filter(subtype == subtype_name) %>%
    pull(patient)
  sample_sets[[subtype_name]] <- all_samples[all_patients %in% subtype_patients]
}

message("[subtype pilot] Samples by group:")
print(data.frame(group = names(sample_sets), n_samples = vapply(sample_sets, length, integer(1))))

message("[subtype pilot] Recalculating ALD/DLD for whole cohort and subtypes ...")
all_calc <- lapply(names(sample_sets), function(label) {
  message("[subtype pilot] ", label, " | samples: ", length(sample_sets[[label]]))
  calc_ald_dld_for_samples(label, sample_sets[[label]], inputs)
})
names(all_calc) <- names(sample_sets)

gene_long <- bind_rows(lapply(all_calc, `[[`, "gene_table"))
summary_basic <- bind_rows(lapply(all_calc, `[[`, "summary_table"))

whole_hits <- gene_long %>%
  filter(analysis_label == whole_label, hit) %>%
  select(state, Gene) %>%
  distinct()

subtype_labels <- setdiff(names(sample_sets), whole_label)

retention_summary <- bind_rows(lapply(subtype_labels, function(label) {
  bind_rows(lapply(c("ALD", "DLD"), function(state_name) {
    whole_state_hits <- whole_hits %>%
      filter(state == state_name) %>%
      pull(Gene)
    subtype_hits <- gene_long %>%
      filter(analysis_label == label, state == state_name, hit) %>%
      pull(Gene)

    data.frame(
      subtype = label,
      state = state_name,
      n_samples = length(sample_sets[[label]]),
      whole_hit_genes = length(whole_state_hits),
      subtype_hit_genes = length(subtype_hits),
      retained_whole_hits = length(intersect(whole_state_hits, subtype_hits)),
      lost_whole_hits = length(setdiff(whole_state_hits, subtype_hits)),
      new_subtype_hits = length(setdiff(subtype_hits, whole_state_hits)),
      retained_pct_of_whole_hits = ifelse(
        length(whole_state_hits) > 0,
        100 * length(intersect(whole_state_hits, subtype_hits)) / length(whole_state_hits),
        NA_real_
      ),
      stringsAsFactors = FALSE
    )
  }))
}))

summary_out <- summary_basic %>%
  mutate(subtype = ifelse(analysis_label == whole_label, whole_label, analysis_label)) %>%
  select(subtype, state, n_samples, eligible_genes, hit_genes, cancer_score) %>%
  left_join(retention_summary, by = c("subtype", "state", "n_samples"))

gene_wide <- gene_long %>%
  mutate(
    percent_col = paste0(analysis_label, "_percentage"),
    hit_col = paste0(analysis_label, "_hit"),
    eligible_col = paste0(analysis_label, "_eligible")
  )

gene_percent <- gene_wide %>%
  select(state, Gene, percent_col, percentage) %>%
  pivot_wider(names_from = percent_col, values_from = percentage)

gene_hit <- gene_wide %>%
  select(state, Gene, hit_col, hit) %>%
  pivot_wider(names_from = hit_col, values_from = hit)

gene_eligible <- gene_wide %>%
  select(state, Gene, eligible_col, eligible_gene) %>%
  pivot_wider(names_from = eligible_col, values_from = eligible_gene)

gene_out <- gene_percent %>%
  full_join(gene_hit, by = c("state", "Gene")) %>%
  full_join(gene_eligible, by = c("state", "Gene")) %>%
  arrange(state, Gene)

write.csv(summary_out, summary_file, row.names = FALSE)
write.csv(gene_out, gene_file, row.names = FALSE)

message("[subtype pilot] Saved retention summary: ", summary_file)
message("[subtype pilot] Saved retention gene-level table: ", gene_file)

message("[subtype pilot] Running subtype log-rank comparisons ...")
cdef <- cancers[[pilot_cohort]]
clinical_inputs <- build_survival_clinical_os(
  cohort = pilot_cohort,
  cdef = cdef,
  force_rebuild = FALSE
)
surv_df <- clinical_inputs$surv_df

whole_logrank <- bind_rows(
  load_whole_logrank(pilot_cohort, "ALD", inputs, surv_df, all_samples),
  load_whole_logrank(pilot_cohort, "DLD", inputs, surv_df, all_samples)
)

survival_rows <- bind_rows(lapply(subtype_labels, function(subtype_name) {
  bind_rows(lapply(c("ALD", "DLD"), function(state_name) {
    subtype_hit_col <- paste0(subtype_name, "_hit")
    subtype_pct_col <- paste0(subtype_name, "_percentage")
    whole_hit_col <- paste0(whole_label, "_hit")
    whole_pct_col <- paste0(whole_label, "_percentage")

    state_retention <- gene_out %>%
      filter(state == state_name) %>%
      transmute(
        subtype = subtype_name,
        state = state_name,
        Gene = Gene,
        whole_ald_dld_percent = .data[[whole_pct_col]],
        subtype_ald_dld_percent = .data[[subtype_pct_col]],
        whole_ald_dld_hit = as_logical_flag(.data[[whole_hit_col]]),
        subtype_ald_dld_hit = as_logical_flag(.data[[subtype_hit_col]])
      ) %>%
      filter(whole_ald_dld_hit %in% TRUE | subtype_ald_dld_hit %in% TRUE) %>%
      mutate(
        ald_dld_hit_comparison_status = case_when(
          whole_ald_dld_hit %in% TRUE & subtype_ald_dld_hit %in% TRUE ~ "retained_ALD_DLD_hit",
          whole_ald_dld_hit %in% TRUE & !(subtype_ald_dld_hit %in% TRUE) ~ "lost_whole_ALD_DLD_hit",
          !(whole_ald_dld_hit %in% TRUE) & subtype_ald_dld_hit %in% TRUE ~ "new_subtype_ALD_DLD_hit",
          TRUE ~ "not_compared"
        )
      )

    if (nrow(state_retention) == 0) return(state_retention)

    subtype_survival <- lapply(seq_len(nrow(state_retention)), function(i) {
      row <- state_retention[i, ]

      if (!isTRUE(row$subtype_ald_dld_hit)) {
        return(data.frame(
          subtype_n_group = NA_integer_,
          subtype_n_reference = NA_integer_,
          subtype_n_events_group = NA_integer_,
          subtype_n_events_reference = NA_integer_,
          subtype_logrank_chisq = NA_real_,
          subtype_logrank_p_value = NA_real_,
          subtype_logrank_direction = NA_character_,
          subtype_logrank_status = "not_tested_not_subtype_ALD_DLD_hit",
          subtype_logrank_message = NA_character_,
          stringsAsFactors = FALSE
        ))
      }

      run_subtype_logrank(
        gene = row$Gene,
        state = row$state,
        subtype = row$subtype,
        sample_sets = sample_sets,
        inputs = inputs,
        surv_df = surv_df
      )
    }) %>%
      bind_rows()

    bind_cols(state_retention, subtype_survival)
  }))
}))

survival_out <- survival_rows %>%
  group_by(subtype, state) %>%
  mutate(
    subtype_logrank_FDR = ifelse(
      subtype_logrank_status == "ok" & !is.na(subtype_logrank_p_value),
      p.adjust(subtype_logrank_p_value, method = "BH"),
      NA_real_
    )
  ) %>%
  ungroup() %>%
  left_join(whole_logrank, by = c("state", "Gene")) %>%
  select(
    subtype,
    state,
    Gene,
    ald_dld_hit_comparison_status,
    whole_ald_dld_percent,
    subtype_ald_dld_percent,
    whole_ald_dld_hit,
    subtype_ald_dld_hit,
    subtype_n_group,
    subtype_n_reference,
    subtype_n_events_group,
    subtype_n_events_reference,
    subtype_logrank_chisq,
    subtype_logrank_p_value,
    subtype_logrank_FDR,
    subtype_logrank_direction,
    subtype_logrank_status,
    subtype_logrank_message,
    whole_n_group,
    whole_n_reference,
    whole_n_events_group,
    whole_n_events_reference,
    whole_logrank_chisq,
    whole_logrank_p_value,
    whole_logrank_FDR,
    whole_logrank_direction,
    whole_logrank_status
  ) %>%
  arrange(subtype, state, ald_dld_hit_comparison_status, subtype_logrank_FDR, whole_logrank_FDR, Gene)

write.csv(survival_out, survival_file, row.names = FALSE)

message("[subtype pilot] Saved survival comparison: ", survival_file)
message("[subtype pilot] Survival rows: ", nrow(survival_out))
message("[subtype pilot] Subtype log-rank FDR < 0.05: ", sum(survival_out$subtype_logrank_FDR < 0.05, na.rm = TRUE))
message("[subtype pilot] Done.")
