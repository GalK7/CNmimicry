# Purpose:
#   Classify ALD/DLD gene-sample pairs and summarize their burdens for one
#   cohort and prevalence threshold, using the supplied CNV, GISTIC and
#   amplified-group median-TPM eligibility settings.
#
# Inputs:
#   - Cohort, paths and analysis settings supplied by run_all.R from config.R.
#   - Prepared expression/CNV matrices and optional GISTIC gene sets loaded
#     through the data-preparation helpers.
#
# Outputs:
#   - RDS files containing gene sets, ALD expression eligibility, ALD/DLD calls,
#     sample-level burdens and diploid counts, plus summary.csv.
#   - Files organized by cohort and prevalence within the selected results
#     directory, with a separate aldTPMgt<value> branch when the filter is enabled.
#   - The cohort summary returned invisibly to the caller.
#
# Usage:
#   Defines run_core(), called by run_all.R for each configured cohort/threshold.

run_core <- function(cohort, gdc_project, gdc_dir,
                     prev_thr,
                     cache_dir, results_dir,
                     cnv_amp_cutoff = 7,
                     cnv_del_cutoff = 1,
                     cnv_diploid = 2,
                     sample_fraction_cutoff = 0.05,
                     er_pos_only = FALSE,
                     use_cnv_median_normalization = FALSE,
                     cnv_norm_stat = "median",
                     cnv_amp_cutoff_norm = 2,
                     cnv_del_cutoff_norm = 0.67,
                     cnv_diploid_norm = 1,
                     cnv_diploid_tol_norm = 0.05,
                     use_gistic_thresholded_intersection = TRUE,
                     use_gistic_significant = FALSE,
                     gdac_run_id = "analyses__2016_01_28",
                     gistic_sig_q_cutoff = 0.25,
                     gistic_sig_gene_source = "wide_peak",
                     gdac_resources_dir = NULL,
                     ald_min_amplified_median_tpm = NULL) {
  
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
  
  # Packages
  req_pkgs <- c("dplyr", "tidyr", "tibble", "readr")
  for (p in req_pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) stop("Package required: ", p)
  }
  
  # ---- Load prepared data + optional thresholded GISTIC sets (cached/reused) ----
  prepared <- prepare_data(
    cohort = cohort, gdc_project = gdc_project, gdc_dir = gdc_dir,
    cache_dir = cache_dir,
    er_pos_only = er_pos_only
  )
  
  gistic <- NULL
  if (isTRUE(use_gistic_thresholded_intersection)) {
    gistic <- get_gistic_sets(
      cohort = cohort, cache_dir = cache_dir, prev_thr = prev_thr
    )
  }
  
  expr_mat <- prepared$expr_mat
  cnv_mat  <- prepared$cnv_mat
  
  # ---- Optional CNV normalization (per-sample) ----
  cnv_use <- cnv_mat
  
  if (isTRUE(use_cnv_median_normalization)) {
    denom <- apply(cnv_mat, 2, median, na.rm = TRUE)
    denom[is.na(denom) | denom == 0] <- NA_real_
    cnv_use <- sweep(cnv_mat, 2, denom, "/")
  }
  
  amp_cutoff_use <- if (isTRUE(use_cnv_median_normalization)) cnv_amp_cutoff_norm else cnv_amp_cutoff
  del_cutoff_use <- if (isTRUE(use_cnv_median_normalization)) cnv_del_cutoff_norm else cnv_del_cutoff
  
  is_diploid <- function(x) {
    if (isTRUE(use_cnv_median_normalization)) {
      abs(x - cnv_diploid_norm) <= cnv_diploid_tol_norm   # [0.95, 1.05]
    } else {
      x == cnv_diploid
    }
  }
  
  # ---- Load GISTIC significant data ----
  gistic_sig <- NULL
  if (isTRUE(use_gistic_significant)) {
    if (is.null(gdac_resources_dir)) {
      stop("use_gistic_significant=TRUE but gdac_resources_dir is NULL")
    }
    gistic_sig <- get_gistic_significant_sets(
      cohort = cohort,
      cache_dir = cache_dir,
      gdac_resources_dir = gdac_resources_dir,
      gdac_run_id = gdac_run_id,
      gistic_sig_q_cutoff = gistic_sig_q_cutoff,
      gistic_sig_gene_source = gistic_sig_gene_source
    )
  }
  
  # ---- Output directory ----
  thr_tag <- format(prev_thr, nsmall = 2)
  if (ald_tpm_filter_enabled) {
    results_dir <- file.path(results_dir, ald_tpm_tag)
  }
  out_dir <- file.path(results_dir, cohort, paste0("prev_", thr_tag))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  message("[run_core] ", ald_tpm_tag, " | Output: ", out_dir)
  
  # =========================================================
  # ALD: Amplified (CN>=7) vs Diploid (CN==2) + outliers
  # =========================================================
  binary_amp <- cnv_use >= amp_cutoff_use
  amp_pct <- rowSums(binary_amp, na.rm = TRUE) / ncol(binary_amp)
  amplified_genes <- names(amp_pct)[amp_pct > prev_thr]
  
  # Mandatory GISTIC intersection (amps)
  if (isTRUE(use_gistic_thresholded_intersection)) {
    amplified_genes <- intersect(amplified_genes, gistic$gistic_amp_genes)
  }
  
  # Keep only amplified genes present in expression matrix (double check for safety)
  amplified_genes <- intersect(amplified_genes, rownames(expr_mat))
  
  # Optional: significant GISTIC genes intersection (amps)
  if (isTRUE(use_gistic_significant) && !is.null(gistic_sig)) {
    amplified_genes <- intersect(amplified_genes, gistic_sig$sig_amp_genes)
  }
  
  
  # True amplification eligibility remains independent of the ALD expression filter.
  ald_eligible_genes <- amplified_genes
  amplified_expression_eligibility <- tibble::tibble(
    Gene = character(), median_amplified = numeric(), ald_eligible = logical()
  )
  
  if (length(amplified_genes) == 0) {
    
    # Empty ALD outputs (prevents pivot_longer crash)
    outliers_percentage_ald <- dplyr::tibble(
      Gene = character(),
      outliers_count = integer(),
      total_non_amplified = integer(),
      percentage = numeric()
    )
    special_samples_filtered <- character(0)
    
    ald_per_sample <- tibble::tibble(
      Sample = colnames(expr_mat),
      n_ald_genes = 0L,
      ald_frac = NA_real_   # no ALD-eligible gene denominator
    )
    
    ald_calls_long <- tibble::tibble(Sample = character(), Gene = character())
    
    
  } else {
    binary_amp <- binary_amp[amplified_genes, , drop = FALSE]
    
    # Build long combined_df_ald for ALD
    amp_status_df <- as.data.frame(t(binary_amp)) |>
      tibble::rownames_to_column("Sample") |>
      tidyr::pivot_longer(-Sample, names_to = "Gene", values_to = "Amplification_Status")
    
    expr_amp_df <- as.data.frame(t(expr_mat[amplified_genes, , drop = FALSE])) |>
      tibble::rownames_to_column("Sample") |>
      tidyr::pivot_longer(-Sample, names_to = "Gene", values_to = "Expression")
    
    cnv_amp_df <- as.data.frame(t(cnv_use[amplified_genes, , drop = FALSE])) |>
      tibble::rownames_to_column("Sample") |>
      tidyr::pivot_longer(-Sample, names_to = "Gene", values_to = "cnv_value")
    
    combined_df_ald <- expr_amp_df |>
      dplyr::inner_join(amp_status_df, by = c("Sample", "Gene")) |>
      dplyr::inner_join(cnv_amp_df,     by = c("Sample", "Gene"))
    
    gene_stats_ald <- combined_df_ald |>
      dplyr::filter(Amplification_Status == TRUE) |>
      dplyr::group_by(Gene) |>
      dplyr::summarise(median_amplified = median(Expression, na.rm = TRUE), .groups = "drop")
    
    amplified_expression_eligibility <- tibble::tibble(Gene = amplified_genes) |>
      dplyr::left_join(gene_stats_ald, by = "Gene") |>
      dplyr::mutate(
        ald_eligible = if (ald_tpm_filter_enabled) {
          is.finite(median_amplified) & median_amplified > ald_min_amplified_median_tpm
        } else {
          TRUE
        }
      )
    ald_eligible_genes <- amplified_expression_eligibility$Gene[
      amplified_expression_eligibility$ald_eligible
    ]
    
    combined_df_ald <- combined_df_ald |>
      dplyr::left_join(gene_stats_ald, by = "Gene")
    if (ald_tpm_filter_enabled) {
      combined_df_ald <- combined_df_ald |>
        dplyr::filter(Gene %in% ald_eligible_genes)
    }
    
    non_amp_outliers <- combined_df_ald |>
      dplyr::filter(Amplification_Status == FALSE,
                    is_diploid(cnv_value),
                    Expression > median_amplified)
    
    ald_calls_long <- non_amp_outliers |>
      dplyr::select(Sample, Gene) |>
      dplyr::distinct()
    
    total_non_amp <- combined_df_ald |>
      dplyr::filter(Amplification_Status == FALSE, is_diploid(cnv_value)) |>
      dplyr::group_by(Gene) |>
      dplyr::summarise(total_non_amplified = dplyr::n(), .groups = "drop")
    
    outliers_percentage_ald <- non_amp_outliers |>
      dplyr::group_by(Gene) |>
      dplyr::summarise(outliers_count = dplyr::n(), .groups = "drop") |>
      dplyr::left_join(total_non_amp, by = "Gene") |>
      dplyr::mutate(percentage = 100 * outliers_count / total_non_amplified) |>
      dplyr::arrange(dplyr::desc(percentage))

    
    # Sample-level ALD uses the TPM-eligible gene denominator when the filter is enabled.
    cutoff_ald <- sample_fraction_cutoff * length(ald_eligible_genes)
    
    special_matrix <- non_amp_outliers |>
      dplyr::distinct(Sample, Gene) |>
      dplyr::mutate(value = 1) |>
      tidyr::pivot_wider(names_from = Gene, values_from = value, values_fill = 0) |>
      tibble::column_to_rownames("Sample")
    
    # Per-sample ALD burden: keep all samples and count zero-ALD samples explicitly.
    ald_per_sample <- tibble::tibble(Sample = colnames(expr_mat)) |>
      dplyr::left_join(
        non_amp_outliers |>
          dplyr::group_by(Sample) |>
          dplyr::summarise(n_ald_genes = dplyr::n_distinct(Gene), .groups = "drop"),
        by = "Sample"
      ) |>
      dplyr::mutate(
        n_ald_genes = dplyr::coalesce(n_ald_genes, 0L),
        ald_frac = if (length(ald_eligible_genes) > 0) n_ald_genes / length(ald_eligible_genes) else NA_real_
      )
    
    
    if (nrow(special_matrix) == 0) {
      special_samples_filtered <- character(0)
    } else {
      al_count <- rowSums(special_matrix)
      special_samples_filtered <- names(al_count)[al_count > cutoff_ald]  # STRICT >
    }
  }
  
  # =========================================================
  # DLD: Deleted (CN<=1) vs Diploid (CN==2) + outliers
  # =========================================================
  binary_del <- cnv_use <= del_cutoff_use
  del_pct <- rowSums(binary_del, na.rm = TRUE) / ncol(binary_del)
  deleted_genes <- names(del_pct)[del_pct > prev_thr]
  
  # Mandatory GISTIC intersection (dels)
  if (isTRUE(use_gistic_thresholded_intersection)) {
    deleted_genes <- intersect(deleted_genes, gistic$gistic_del_genes)
  }
  
  # Keep only deleted genes present in expression matrix (double check for safety)
  deleted_genes <- intersect(deleted_genes, rownames(expr_mat))
  
  # Optional: significant GISTIC genes intersection (dels)
  if (isTRUE(use_gistic_significant) && !is.null(gistic_sig)) {
    deleted_genes <- intersect(deleted_genes, gistic_sig$sig_del_genes)
  }
  
  
  if (length(deleted_genes) == 0) {
    
    # Empty DLD outputs (prevents pivot_longer crash)
    outliers_percentage_dld <- dplyr::tibble(
      Gene = character(),
      total_diploid = integer(),
      outliers_count = integer(),
      percentage = numeric()
    )
    
    dld_per_sample <- dplyr::tibble(
      Sample = character(),
      n_dld_genes = integer(),
      dld_frac = numeric()
    )
    
    dld_calls_long <- tibble::tibble(Sample = character(), Gene = character())
    
    n_samples_gt5pct_dld <- 0L
    
  } else {
    binary_del <- binary_del[deleted_genes, , drop = FALSE]
    
    del_status_df <- as.data.frame(t(binary_del)) |>
      tibble::rownames_to_column("Sample") |>
      tidyr::pivot_longer(-Sample, names_to = "Gene", values_to = "Deletion_Status")
    
    expr_del_df <- as.data.frame(t(expr_mat[deleted_genes, , drop = FALSE])) |>
      tibble::rownames_to_column("Sample") |>
      tidyr::pivot_longer(-Sample, names_to = "Gene", values_to = "Expression")
    
    cnv_del_df <- as.data.frame(t(cnv_use[deleted_genes, , drop = FALSE])) |>
      tibble::rownames_to_column("Sample") |>
      tidyr::pivot_longer(-Sample, names_to = "Gene", values_to = "cnv_value")
    
    combined_df_dld <- expr_del_df |>
      dplyr::inner_join(del_status_df, by = c("Sample", "Gene")) |>
      dplyr::inner_join(cnv_del_df,    by = c("Sample", "Gene"))
    
    gene_stats_dld <- combined_df_dld |>
      dplyr::filter(Deletion_Status == TRUE) |>
      dplyr::group_by(Gene) |>
      dplyr::summarise(median_deleted = median(Expression, na.rm = TRUE), .groups = "drop")
    
    combined_df_dld <- combined_df_dld |>
      dplyr::left_join(gene_stats_dld, by = "Gene")
    
    diploid_low_outliers_dld <- combined_df_dld |>
      dplyr::filter(Deletion_Status == FALSE,
                    is_diploid(cnv_value),
                    Expression < median_deleted)
    
    dld_calls_long <- diploid_low_outliers_dld |>
      dplyr::select(Sample, Gene) |>
      dplyr::distinct()

    
    outliers_percentage_dld <- combined_df_dld |>
      dplyr::group_by(Gene) |>
      dplyr::summarise(
        total_diploid  = sum(Deletion_Status == FALSE & is_diploid(cnv_value), na.rm = TRUE),
        outliers_count = sum(Deletion_Status == FALSE & is_diploid(cnv_value) & Expression < median_deleted, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(percentage = dplyr::if_else(total_diploid > 0, 100 * outliers_count / total_diploid, NA_real_)) |>
      dplyr::arrange(dplyr::desc(percentage))
    
    # Sample-level DLD (>5% of deleted genes), STRICT >
    n_deleted_genes <- length(deleted_genes)
    
    dld_per_sample <- tibble::tibble(Sample = colnames(expr_mat)) |>
      dplyr::left_join(
        diploid_low_outliers_dld |>
          dplyr::group_by(Sample) |>
          dplyr::summarise(n_dld_genes = dplyr::n_distinct(Gene), .groups = "drop"),
        by = "Sample"
      ) |>
      dplyr::mutate(
        n_dld_genes = dplyr::coalesce(n_dld_genes, 0L),
        dld_frac = if (n_deleted_genes > 0) n_dld_genes / n_deleted_genes else NA_real_
      )
    
    n_samples_gt5pct_dld <- sum(dld_per_sample$dld_frac > sample_fraction_cutoff, na.rm = TRUE)
  }
  
  
  # =========================================================
  # Sample-level diploid classification (genome-wide; ASCAT3 CN only)
  # =========================================================
  diploid_mask <- is_diploid(cnv_use)  # genes x samples logical
  diploid_frac <- colMeans(diploid_mask, na.rm = TRUE)
  
  diploid_sample_cutoff <- 0.90
  
  diploid_per_sample <- tibble::tibble(
    Sample = names(diploid_frac),
    diploid_frac_genome = as.numeric(diploid_frac),
    is_diploid_sample = diploid_frac >= diploid_sample_cutoff
  )
  
  # =========================================================
  # Cohort-level summary
  # =========================================================
  n_amp_genes <- length(amplified_genes)
  n_ald_eligible_genes <- length(ald_eligible_genes)
  n_del_genes <- length(deleted_genes)
  
  n_genes_gt20_ald <- sum(outliers_percentage_ald$percentage > 20, na.rm = TRUE)
  n_genes_gt50_ald <- sum(outliers_percentage_ald$percentage > 50, na.rm = TRUE)
  n_genes_gt20_dld <- sum(outliers_percentage_dld$percentage > 20, na.rm = TRUE)
  n_genes_gt50_dld <- sum(outliers_percentage_dld$percentage > 50, na.rm = TRUE)
  
  top10_ald <- paste(head(outliers_percentage_ald$Gene, 10), collapse = ", ")
  top10_dld <- paste(head(outliers_percentage_dld$Gene, 10), collapse = ", ")
  
  # Cancer score = (# genes above cutoff) / (# event genes)
  score20_ald <- if (n_ald_eligible_genes > 0) n_genes_gt20_ald / n_ald_eligible_genes else NA_real_
  score50_ald <- if (n_ald_eligible_genes > 0) n_genes_gt50_ald / n_ald_eligible_genes else NA_real_
  score20_dld <- if (n_del_genes > 0) n_genes_gt20_dld / n_del_genes else NA_real_
  score50_dld <- if (n_del_genes > 0) n_genes_gt50_dld / n_del_genes else NA_real_
  
  summary_df <- data.frame(
    cancer = cohort,
    threshold = prev_thr,
    n_amp_genes = n_amp_genes,
    n_del_genes = n_del_genes,
    n_samples_gt5pct_ald = length(special_samples_filtered),
    n_samples_gt5pct_dld = n_samples_gt5pct_dld,
    n_genes_gt20_ald = n_genes_gt20_ald,
    n_genes_gt20_dld = n_genes_gt20_dld,
    top10_genes_ald = top10_ald,
    top10_genes_dld = top10_dld,
    cancer_score20_ald = score20_ald,
    cancer_score50_ald = score50_ald,
    cancer_score20_dld = score20_dld,
    cancer_score50_dld = score50_dld,
    use_gistic_significant = isTRUE(use_gistic_significant),
    use_cnv_median_normalization = isTRUE(use_cnv_median_normalization),
    cnv_norm_stat = if (isTRUE(use_cnv_median_normalization)) cnv_norm_stat else NA_character_,
    gdac_run_id = if (isTRUE(use_gistic_significant)) gdac_run_id else NA_character_,
    gistic_sig_q_cutoff = if (isTRUE(use_gistic_significant)) gistic_sig_q_cutoff else NA_real_,
    gistic_sig_gene_source = if (isTRUE(use_gistic_significant)) gistic_sig_gene_source else NA_character_,
    ald_tpm_filter_enabled = ald_tpm_filter_enabled,
    ald_min_amplified_median_tpm = if (ald_tpm_filter_enabled) ald_min_amplified_median_tpm else NA_real_,
    ald_tpm_tag = ald_tpm_tag,
    n_ald_eligible_genes = n_ald_eligible_genes,
    n_amp_genes_excluded_from_ald = n_amp_genes - n_ald_eligible_genes,
    ald_gene_denominator = "ald_eligible_genes",
    stringsAsFactors = FALSE
  )
  
  # ---- Save outputs ----
  saveRDS(amplified_genes,     file.path(out_dir, "amplified_genes.rds"))
  saveRDS(ald_eligible_genes, file.path(out_dir, "ald_eligible_genes.rds"))
  saveRDS(amplified_expression_eligibility, file.path(out_dir, "amplified_expression_eligibility.rds"))
  saveRDS(deleted_genes,    file.path(out_dir, "deleted_genes.rds"))
  saveRDS(outliers_percentage_ald, file.path(out_dir, "outliers_percentage_ald.rds"))
  saveRDS(outliers_percentage_dld, file.path(out_dir, "outliers_percentage_dld.rds"))
  saveRDS(ald_per_sample, file.path(out_dir, "ald_per_sample.rds"))
  saveRDS(ald_calls_long, file.path(out_dir, "ald_calls_long.rds"))
  saveRDS(dld_per_sample, file.path(out_dir, "dld_per_sample.rds"))
  saveRDS(dld_calls_long, file.path(out_dir, "dld_calls_long.rds"))
  saveRDS(diploid_per_sample, file.path(out_dir, "diploid_per_sample.rds"))

  readr::write_csv(summary_df, file.path(out_dir, "summary.csv"))
  
  invisible(summary_df)
}
