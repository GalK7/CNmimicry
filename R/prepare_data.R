# Purpose:
#   Load or prepare aligned TCGA expression and copy-number matrices for
#   the core workflow, applying the upstream expression filter and optional
#   BRCA ER-positive sample restriction when preparing data.
#
# Inputs:
#   - Cohort, GDC project, data/cache paths and filtering arguments.
#   - An existing prepared-data cache, or GDC TPM expression, ASCAT3 gene-level
#     copy-number data and clinical data when ER-positive filtering is requested.
#
# Outputs:
#   - A prepared-data list containing aligned matrices and sample/gene metadata.
#   - A cohort-specific RDS cache under cache_dir, with an ERpos suffix when used.
#
# Usage:
#   Defines prepare_data(), called by run_core(); existing caches are reused.

prepare_data <- function(cohort, gdc_project, gdc_dir,
                         cache_dir,
                         er_pos_only = FALSE,
                         expr_min_tpm = 1,
                         expr_min_frac = 0.05) {
  
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    stop("Package 'TCGAbiolinks' is required.")
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' is required.")
  }
  
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(gdc_dir, showWarnings = FALSE, recursive = TRUE)
  
  gdc_dir_norm <- normalizePath(gdc_dir, winslash = "/", mustWork = TRUE)
  
  # IMPORTANT: cache must distinguish BRCA ER+ vs all-samples
  suffix <- if (isTRUE(er_pos_only)) "_ERpos" else ""
  cache_file <- file.path(cache_dir, paste0(cohort, suffix, "_prepared.rds"))
  
  if (file.exists(cache_file)) {
    message("[prepare_data] Using cached prepared data: ", cache_file)
    return(readRDS(cache_file))
  }
  
  message("[prepare_data] Preparing data for: ", cohort, if (isTRUE(er_pos_only)) " (ER+ only)" else "")
  message("[prepare_data] GDC directory: ", gdc_dir_norm)
  
  # -----------------------------
  # Optional BRCA ER-positive sample selection
  # -----------------------------
  er_positive_samples <- NULL
  
  if (isTRUE(er_pos_only)) {
    q_clin <- TCGAbiolinks::GDCquery(
      project = gdc_project,
      data.category = "Clinical",
      data.type = "Clinical Supplement",
      data.format = "BCR Biotab"
    )
    TCGAbiolinks::GDCdownload(q_clin, directory = gdc_dir)
    clinical_data <- TCGAbiolinks::GDCprepare(q_clin, directory = gdc_dir)
    
    clinical_df <- clinical_data[["clinical_patient_brca"]]
    if (is.null(clinical_df)) {
      stop("ER+ requested but clinical_patient_brca was not found in clinical_data.")
    }
    if (!all(c("er_status_by_ihc", "bcr_patient_barcode") %in% colnames(clinical_df))) {
      stop("ER+ requested but clinical_patient_brca is missing expected columns er_status_by_ihc / bcr_patient_barcode.")
    }
    
    er_positive_samples <- clinical_df[clinical_df$er_status_by_ihc == "Positive", "bcr_patient_barcode"]$bcr_patient_barcode
    er_positive_samples <- unique(as.character(er_positive_samples))
    er_positive_samples <- er_positive_samples[!is.na(er_positive_samples) & nzchar(er_positive_samples)]
    
    if (length(er_positive_samples) == 0) {
      stop("ER+ filtering produced 0 samples. Check clinical er_status_by_ihc values.")
    }
  }
  
  # -----------------------------
  # Expression (TPM)
  # -----------------------------
  if (isTRUE(er_pos_only)) {
    q_expr <- TCGAbiolinks::GDCquery(
      project = gdc_project,
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      data.format = "TPM",
      barcode = er_positive_samples
    )
  } else {
    q_expr <- TCGAbiolinks::GDCquery(
      project = gdc_project,
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      data.format = "TPM"
    )
  }
  TCGAbiolinks::GDCdownload(q_expr, directory = gdc_dir)
  se_expr <- TCGAbiolinks::GDCprepare(q_expr, directory = gdc_dir)
  
  expr_mat <- SummarizedExperiment::assay(se_expr, "tpm_unstrand")
  rownames(expr_mat) <- SummarizedExperiment::rowRanges(se_expr)$gene_name
  colnames(expr_mat) <- se_expr$patient
  
  # -----------------------------
  # CNV (ASCAT3 gene-level copy number)
  # -----------------------------
  if (isTRUE(er_pos_only)) {
    q_cnv <- TCGAbiolinks::GDCquery(
      project = gdc_project,
      data.category = "Copy Number Variation",
      data.type = "Gene Level Copy Number",
      access = "open",
      workflow.type = "ASCAT3",
      barcode = er_positive_samples
    )
  } else {
    q_cnv <- TCGAbiolinks::GDCquery(
      project = gdc_project,
      data.category = "Copy Number Variation",
      data.type = "Gene Level Copy Number",
      access = "open",
      workflow.type = "ASCAT3"
    )
  }
  TCGAbiolinks::GDCdownload(q_cnv, directory = gdc_dir)
  se_cnv <- TCGAbiolinks::GDCprepare(q_cnv, directory = gdc_dir)
  
  cnv_mat <- SummarizedExperiment::assay(se_cnv, "copy_number")
  rownames(cnv_mat) <- SummarizedExperiment::rowRanges(se_cnv)$gene_name
  colnames(cnv_mat) <- se_cnv$patient
  
  # -----------------------------
  # Align samples and genes
  # -----------------------------
  common_samples <- intersect(colnames(expr_mat), colnames(cnv_mat))
  expr_mat <- expr_mat[, common_samples, drop = FALSE]
  cnv_mat  <- cnv_mat[, common_samples, drop = FALSE]
  
  common_genes <- intersect(rownames(expr_mat), rownames(cnv_mat))
  expr_mat <- expr_mat[common_genes, , drop = FALSE]
  cnv_mat  <- cnv_mat[common_genes, , drop = FALSE]
  
  # -----------------------------
  # Valid-gene filter: TPM > expr_min_tpm in > expr_min_frac samples
  # -----------------------------
  sample_threshold <- expr_min_frac * ncol(expr_mat)
  samples_with_expr <- rowSums(expr_mat > expr_min_tpm, na.rm = TRUE)
  valid_genes <- rownames(expr_mat)[samples_with_expr > sample_threshold]
  
  expr_mat <- expr_mat[valid_genes, , drop = FALSE]
  cnv_mat  <- cnv_mat[valid_genes, , drop = FALSE]
  
  prepared <- list(
    cohort = cohort,
    gdc_project = gdc_project,
    gdc_dir = gdc_dir_norm,
    er_pos_only = er_pos_only,
    expr_mat = expr_mat,
    cnv_mat = cnv_mat,
    common_samples = common_samples,
    valid_genes = valid_genes,
    prepared_at = as.character(Sys.time())
  )
  
  saveRDS(prepared, cache_file)
  message("[prepare_data] Saved cache: ", cache_file)
  
  prepared
}
