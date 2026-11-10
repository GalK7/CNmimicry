# Purpose:
#   Prepare promoter methylation data and provide ALD/DLD screening functions.
#
# Inputs:
#   - config.R settings, core expression/CNV data, ALD/DLD calls, and gene lists.
#   - TCGA methylation data and probe annotations.
#
# Outputs:
#   - Promoter matrices, matched sample groups, Wilcoxon results, and input summaries.
#   - Downloaded data under analysis/methylation/GDCdata and promoter matrix
#     caches under analysis/methylation/outputs; callers write screening results.

source(file.path("R", "prepare_data.R"))

# ------------------------------------------------------------------------------
# Resolve output paths
# ------------------------------------------------------------------------------

resolve_methylation_output_dir <- function(cohort) {
  out_dir <- file.path(workflow_root, "analysis", "methylation", "outputs", cohort)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir
}

# ------------------------------------------------------------------------------
# Resolve CNmimicry result directory
# ------------------------------------------------------------------------------

get_methylation_ald_tpm_settings <- function() {
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

resolve_methylation_context <- function(cohort, prev_thr = 0.03,
                                        apply_ald_tpm_filter = FALSE) {
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
  
  # Both screens read the configured core run; only ALD outputs carry the TPM tag.
  core_tpm <- get_methylation_ald_tpm_settings()
  ald_tpm <- if (apply_ald_tpm_filter) core_tpm else
    list(enabled = FALSE, tag = "aldTPMoff", suffix = "")
  if (core_tpm$enabled) results_base <- file.path(results_base, core_tpm$tag)
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

# ------------------------------------------------------------------------------
# Load methylation beta data
# ------------------------------------------------------------------------------

load_methylation_se <- function(cohort,
                                platform = "Illumina Human Methylation 450",
                                sample_type = "Primary Tumor") {
  
  cdef <- cancers[[cohort]]
  
  methylation_gdc_dir <- file.path(
    workflow_root,
    "analysis", "methylation",
    "GDCdata",
    cohort
  )
  
  dir.create(methylation_gdc_dir, recursive = TRUE, showWarnings = FALSE)
  
  query_met <- tryCatch(
    TCGAbiolinks::GDCquery(
      project = cdef$gdc_project,
      data.category = "DNA Methylation",
      platform = platform,
      data.type = "Methylation Beta Value",
      sample.type = sample_type
    ),
    error = function(e) {
      stop(
        "No matching methylation data for ",
        cohort,
        " using platform='", platform,
        "', sample_type='", sample_type,
        "'. Original error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  TCGAbiolinks::GDCdownload(
    query_met,
    directory = methylation_gdc_dir,
    method = "api",
    files.per.chunk = 50
  )
  
  met_se <- TCGAbiolinks::GDCprepare(
    query_met,
    directory = methylation_gdc_dir
  )
  
  met_se
}

# ------------------------------------------------------------------------------
# Inspect methylation object annotation
# ------------------------------------------------------------------------------

inspect_methylation_annotation <- function(met_se) {
  anno <- as.data.frame(SummarizedExperiment::rowRanges(met_se))
  
  cat("Annotation columns:\n")
  print(colnames(anno))
  
  cat("\nFirst rows:\n")
  print(head(anno))
  
  invisible(anno)
}

# ------------------------------------------------------------------------------
# Collapse duplicate patient columns
# ------------------------------------------------------------------------------

collapse_methylation_by_patient <- function(beta_mat) {
  patient_ids <- substr(colnames(beta_mat), 1, 12)
  patient_groups <- split(seq_along(patient_ids), patient_ids)
  
  collapsed <- sapply(patient_groups, function(idx) {
    if (length(idx) == 1) {
      beta_mat[, idx]
    } else {
      apply(beta_mat[, idx, drop = FALSE], 1, median, na.rm = TRUE)
    }
  })
  
  collapsed <- as.matrix(collapsed)
  rownames(collapsed) <- rownames(beta_mat)
  collapsed
}

# ------------------------------------------------------------------------------
# Probe annotation helpers
# ------------------------------------------------------------------------------

pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)][1]
  if (is.na(hit)) return(NA_character_)
  hit
}

split_annotation_field <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  strsplit(x, ";", fixed = TRUE)[[1]]
}

make_probe_gene_promoter_map <- function(met_se,
                                         promoter_regions = c("TSS200", "TSS1500", "5'UTR")) {
  
  if (!requireNamespace("minfi", quietly = TRUE)) {
    stop("Package 'minfi' is required. Install it with BiocManager::install('minfi')")
  }
  
  if (!requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
    stop("Package 'IlluminaHumanMethylation450kanno.ilmn12.hg19' is required.")
  }
  
  suppressPackageStartupMessages(
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  )
  
  data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  
  anno <- as.data.frame(
    minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  )
  
  anno$probe_id <- rownames(anno)
  
  required_cols <- c("UCSC_RefGene_Name", "UCSC_RefGene_Group")
  missing_cols <- setdiff(required_cols, colnames(anno))
  
  if (length(missing_cols) > 0) {
    stop(
      "The 450K annotation is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      "\nAvailable columns:\n",
      paste(colnames(anno), collapse = ", ")
    )
  }
  
  anno_use <- anno[
    anno$probe_id %in% rownames(SummarizedExperiment::assay(met_se)) &
      !is.na(anno$UCSC_RefGene_Name) &
      anno$UCSC_RefGene_Name != "" &
      !is.na(anno$UCSC_RefGene_Group) &
      anno$UCSC_RefGene_Group != "",
  ]
  
  map_list <- lapply(seq_len(nrow(anno_use)), function(i) {
    
    genes <- split_annotation_field(anno_use$UCSC_RefGene_Name[i])
    regions <- split_annotation_field(anno_use$UCSC_RefGene_Group[i])
    
    genes <- genes[genes != ""]
    regions <- regions[regions != ""]
    
    if (length(genes) == 0 || length(regions) == 0) return(NULL)
    
    if (length(genes) == length(regions)) {
      tmp <- data.frame(
        probe_id = anno_use$probe_id[i],
        gene = genes,
        region = regions,
        stringsAsFactors = FALSE
      )
    } else {
      tmp <- expand.grid(
        probe_id = anno_use$probe_id[i],
        gene = genes,
        region = regions,
        stringsAsFactors = FALSE
      )
    }
    
    tmp[tmp$region %in% promoter_regions, ]
  })
  
  probe_gene_map <- dplyr::bind_rows(map_list) %>%
    dplyr::filter(!is.na(gene), gene != "") %>%
    dplyr::distinct(probe_id, gene, region)
  
  probe_gene_map
}
# ------------------------------------------------------------------------------
# Build or load gene-level promoter methylation matrix
# ------------------------------------------------------------------------------

make_cache_tag <- function(x) {
  tag <- paste(x, collapse = "_")
  tag <- gsub("[^A-Za-z0-9]+", "_", tag)
  tag <- gsub("_+", "_", tag)
  tag <- gsub("^_|_$", "", tag)
  
  if (!nzchar(tag)) {
    tag <- "unspecified"
  }
  
  tag
}

make_promoter_region_tag <- function(promoter_regions) {
  make_cache_tag(promoter_regions)
}

build_promoter_methylation_matrix <- function(cohort,
                                              promoter_regions = c("TSS200", "TSS1500", "5'UTR"),
                                              platform = "Illumina Human Methylation 450",
                                              sample_type = "Primary Tumor",
                                              force_rebuild = FALSE) {
  
  out_dir <- resolve_methylation_output_dir(cohort)
  promoter_tag <- make_promoter_region_tag(promoter_regions)
  platform_tag <- make_cache_tag(platform)
  sample_type_tag <- make_cache_tag(sample_type)
  
  cache_file <- file.path(
    out_dir,
    paste0(
      cohort,
      "_promoter_methylation_matrix_",
      promoter_tag,
      "_",
      platform_tag,
      "_",
      sample_type_tag,
      ".rds"
    )
  )
  
  if (file.exists(cache_file) && !isTRUE(force_rebuild)) {
    message("[methylation] Using cached promoter methylation matrix: ", cache_file)
    return(readRDS(cache_file))
  }
  
  message("[methylation] Loading methylation data for ", cohort, " ...")
  met_se <- load_methylation_se(
    cohort = cohort,
    platform = platform,
    sample_type = sample_type
  )
  
  beta_mat <- SummarizedExperiment::assay(met_se)
  beta_mat <- as.matrix(beta_mat)
  
  message("[methylation] Collapsing methylation samples to patient barcodes ...")
  beta_mat <- collapse_methylation_by_patient(beta_mat)
  
  message("[methylation] Mapping promoter probes to genes ...")
  probe_gene_map <- make_probe_gene_promoter_map(
    met_se = met_se,
    promoter_regions = promoter_regions
  )
  
  probe_gene_map <- probe_gene_map %>%
    dplyr::filter(probe_id %in% rownames(beta_mat))
  
  message("[methylation] Promoter probes retained: ", length(unique(probe_gene_map$probe_id)))
  message("[methylation] Genes with promoter probes: ", length(unique(probe_gene_map$gene)))
  
  probes_by_gene <- split(probe_gene_map$probe_id, probe_gene_map$gene)
  probes_by_gene <- lapply(probes_by_gene, unique)
  
  message("[methylation] Summarizing promoter methylation by gene ...")
  
  promoter_mat <- t(sapply(names(probes_by_gene), function(g) {
    probes <- probes_by_gene[[g]]
    apply(beta_mat[probes, , drop = FALSE], 2, median, na.rm = TRUE)
  }))
  
  promoter_mat <- as.matrix(promoter_mat)
  
  out <- list(
    promoter_methylation_matrix = promoter_mat,
    probe_gene_map = probe_gene_map,
    promoter_regions = promoter_regions,
    platform = platform,
    sample_type = sample_type,
    cohort = cohort
  )
  
  saveRDS(out, cache_file)
  message("[methylation] Saved promoter methylation matrix: ", cache_file)
  
  out
}

# ------------------------------------------------------------------------------
# Load DLD/CNV inputs
# ------------------------------------------------------------------------------

load_dld_methylation_inputs <- function(cohort, prev_thr = 0.03) {
  ctx <- resolve_methylation_context(cohort, prev_thr)
  cdef <- cancers[[cohort]]
  
  dld_calls_path <- file.path(ctx$res_dir, "dld_calls_long.rds")
  del_genes_path <- file.path(ctx$res_dir, "deleted_genes.rds")
  
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
    
    dip_center <- cnv_diploid_norm
    dip_tol <- cnv_diploid_tol_norm
  } else {
    cnv_use <- cnv_mat
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
  
  genes_to_test <- intersect(deleted_genes, rownames(cnv_use))
  
  list(
    ctx = ctx,
    dld_calls_patient = dld_calls_patient,
    cnv_use = cnv_use,
    is_diploid = is_diploid,
    genes_to_test = genes_to_test
  )
}

# ------------------------------------------------------------------------------
# One-gene methylation test: DLD vs other diploids
# ------------------------------------------------------------------------------

run_one_gene_dld_methylation <- function(g,
                                         promoter_mat,
                                         dld_inputs,
                                         min_dld = 5,
                                         min_ref = 10) {
  
  if (!g %in% rownames(promoter_mat)) return(NULL)
  if (!g %in% rownames(dld_inputs$cnv_use)) return(NULL)
  
  dld_patients <- dld_inputs$dld_calls_patient %>%
    dplyr::filter(Gene == g) %>%
    dplyr::pull(patient) %>%
    unique()
  
  gene_cnv <- dld_inputs$cnv_use[g, , drop = TRUE]
  
  dip_patients <- names(gene_cnv)[dld_inputs$is_diploid(gene_cnv)]
  dip_patients <- unique(substr(dip_patients, 1, 12))
  
  ref_patients <- setdiff(dip_patients, dld_patients)
  
  available_patients <- colnames(promoter_mat)
  
  dld_patients <- intersect(dld_patients, available_patients)
  ref_patients <- intersect(ref_patients, available_patients)
  
  if (length(dld_patients) < min_dld || length(ref_patients) < min_ref) {
    return(NULL)
  }
  
  dld_values <- as.numeric(promoter_mat[g, dld_patients, drop = TRUE])
  ref_values <- as.numeric(promoter_mat[g, ref_patients, drop = TRUE])
  
  dld_values <- dld_values[!is.na(dld_values)]
  ref_values <- ref_values[!is.na(ref_values)]
  
  if (length(dld_values) < min_dld || length(ref_values) < min_ref) {
    return(NULL)
  }
  
  wt <- tryCatch(
    wilcox.test(dld_values, ref_values, exact = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(wt)) return(NULL)
  
  data.frame(
    gene = g,
    n_dld = length(dld_values),
    n_other_diploid = length(ref_values),
    median_methylation_dld = median(dld_values, na.rm = TRUE),
    median_methylation_other = median(ref_values, na.rm = TRUE),
    delta_methylation = median(dld_values, na.rm = TRUE) - median(ref_values, na.rm = TRUE),
    p_value = wt$p.value,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# Full DLD promoter methylation screen
# ------------------------------------------------------------------------------

run_dld_promoter_methylation_screen <- function(promoter_mat,
                                                dld_inputs,
                                                min_dld = 5,
                                                min_ref = 10) {
  
  genes_to_test <- intersect(dld_inputs$genes_to_test, rownames(promoter_mat))
  
  res_df <- lapply(genes_to_test, function(g) {
    run_one_gene_dld_methylation(
      g = g,
      promoter_mat = promoter_mat,
      dld_inputs = dld_inputs,
      min_dld = min_dld,
      min_ref = min_ref
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

# ------------------------------------------------------------------------------
# Input-overlap summary
# ------------------------------------------------------------------------------

summarize_dld_methylation_inputs <- function(promoter_mat,
                                             dld_inputs,
                                             res_df = NULL,
                                             promoter_regions = NULL,
                                             platform = NA_character_,
                                             sample_type = NA_character_,
                                             min_dld = NA_integer_,
                                             min_ref = NA_integer_,
                                             delta_methylation_cutoff = NA_real_) {
  
  methylation_patients <- colnames(promoter_mat)
  cnv_patients <- unique(substr(colnames(dld_inputs$cnv_use), 1, 12))
  dld_call_patients <- unique(dld_inputs$dld_calls_patient$patient)
  
  genes_with_promoter <- rownames(promoter_mat)
  genes_to_test <- dld_inputs$genes_to_test
  genes_testable_by_promoter <- intersect(genes_to_test, genes_with_promoter)
  
  data.frame(
    cohort = dld_inputs$ctx$cohort,
    prev_thr = dld_inputs$ctx$prev_thr,
    gistic_thresholded_intersection = dld_inputs$ctx$gthr_tag,
    gistic_significant = dld_inputs$ctx$sig_tag,
    cnv_mode = dld_inputs$ctx$cnv_tag,
    promoter_regions = paste(promoter_regions, collapse = ";"),
    methylation_platform = platform,
    methylation_sample_type = sample_type,
    min_dld = min_dld,
    min_ref = min_ref,
    delta_methylation_cutoff = delta_methylation_cutoff,
    n_methylation_patients = length(unique(methylation_patients)),
    n_cnv_patients = length(unique(cnv_patients)),
    n_dld_call_patients = length(unique(dld_call_patients)),
    n_methylation_cnv_overlap_patients = length(intersect(methylation_patients, cnv_patients)),
    n_methylation_dld_call_overlap_patients = length(intersect(methylation_patients, dld_call_patients)),
    n_deleted_genes_in_cnmimicry_results = length(genes_to_test),
    n_genes_with_promoter_methylation = length(genes_with_promoter),
    n_deleted_genes_with_promoter_methylation = length(genes_testable_by_promoter),
    n_genes_tested_after_min_sample_filters = if (is.null(res_df)) NA_integer_ else nrow(res_df),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# Load ALD/CNV inputs
# ------------------------------------------------------------------------------

load_ald_methylation_inputs <- function(cohort, prev_thr = 0.03) {
  ctx <- resolve_methylation_context(cohort, prev_thr, apply_ald_tpm_filter = TRUE)
  cdef <- cancers[[cohort]]
  
  ald_calls_path <- file.path(ctx$res_dir, "ald_calls_long.rds")
  amp_genes_path <- file.path(ctx$res_dir, "amplified_genes.rds")
  ald_genes_path <- if (ctx$ald_tpm_filter_enabled) {
    file.path(ctx$res_dir, "ald_eligible_genes.rds")
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
    
    dip_center <- cnv_diploid_norm
    dip_tol <- cnv_diploid_tol_norm
  } else {
    cnv_use <- cnv_mat
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
  
  genes_to_test <- intersect(ald_eligible_genes, rownames(cnv_use))
  
  list(
    ctx = ctx,
    ald_calls_patient = ald_calls_patient,
    amplified_genes = amplified_genes,
    cnv_use = cnv_use,
    is_diploid = is_diploid,
    genes_to_test = genes_to_test
  )
}

# ------------------------------------------------------------------------------
# One-gene methylation test: ALD vs other diploids
# ------------------------------------------------------------------------------

run_one_gene_ald_methylation <- function(g,
                                         promoter_mat,
                                         ald_inputs,
                                         min_ald = 10,
                                         min_ref = 10) {
  
  if (!g %in% rownames(promoter_mat)) return(NULL)
  if (!g %in% rownames(ald_inputs$cnv_use)) return(NULL)
  
  ald_patients <- ald_inputs$ald_calls_patient %>%
    dplyr::filter(Gene == g) %>%
    dplyr::pull(patient) %>%
    unique()
  
  gene_cnv <- ald_inputs$cnv_use[g, , drop = TRUE]
  
  dip_patients <- names(gene_cnv)[ald_inputs$is_diploid(gene_cnv)]
  dip_patients <- unique(substr(dip_patients, 1, 12))
  
  ref_patients <- setdiff(dip_patients, ald_patients)
  
  available_patients <- colnames(promoter_mat)
  
  ald_patients <- intersect(ald_patients, available_patients)
  ref_patients <- intersect(ref_patients, available_patients)
  
  if (length(ald_patients) < min_ald || length(ref_patients) < min_ref) {
    return(NULL)
  }
  
  ald_values <- as.numeric(promoter_mat[g, ald_patients, drop = TRUE])
  ref_values <- as.numeric(promoter_mat[g, ref_patients, drop = TRUE])
  
  ald_values <- ald_values[!is.na(ald_values)]
  ref_values <- ref_values[!is.na(ref_values)]
  
  if (length(ald_values) < min_ald || length(ref_values) < min_ref) {
    return(NULL)
  }
  
  wt <- tryCatch(
    wilcox.test(ald_values, ref_values, exact = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(wt)) return(NULL)
  
  data.frame(
    gene = g,
    n_ald = length(ald_values),
    n_other_diploid = length(ref_values),
    median_methylation_ald = median(ald_values, na.rm = TRUE),
    median_methylation_other = median(ref_values, na.rm = TRUE),
    delta_methylation = median(ald_values, na.rm = TRUE) - median(ref_values, na.rm = TRUE),
    p_value = wt$p.value,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# Full ALD promoter methylation screen
# ------------------------------------------------------------------------------

run_ald_promoter_methylation_screen <- function(promoter_mat,
                                                ald_inputs,
                                                min_ald = 10,
                                                min_ref = 10) {
  
  genes_to_test <- intersect(ald_inputs$genes_to_test, rownames(promoter_mat))
  
  res_df <- lapply(genes_to_test, function(g) {
    run_one_gene_ald_methylation(
      g = g,
      promoter_mat = promoter_mat,
      ald_inputs = ald_inputs,
      min_ald = min_ald,
      min_ref = min_ref
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

# ------------------------------------------------------------------------------
# ALD input-overlap summary
# ------------------------------------------------------------------------------

summarize_ald_methylation_inputs <- function(promoter_mat,
                                             ald_inputs,
                                             res_df = NULL,
                                             promoter_regions = NULL,
                                             platform = NA_character_,
                                             sample_type = NA_character_,
                                             min_ald = NA_integer_,
                                             min_ref = NA_integer_,
                                             delta_methylation_cutoff = NA_real_) {
  
  methylation_patients <- colnames(promoter_mat)
  cnv_patients <- unique(substr(colnames(ald_inputs$cnv_use), 1, 12))
  ald_call_patients <- unique(ald_inputs$ald_calls_patient$patient)
  
  genes_with_promoter <- rownames(promoter_mat)
  genes_to_test <- ald_inputs$genes_to_test
  genes_testable_by_promoter <- intersect(genes_to_test, genes_with_promoter)
  
  data.frame(
    cohort = ald_inputs$ctx$cohort,
    prev_thr = ald_inputs$ctx$prev_thr,
    gistic_thresholded_intersection = ald_inputs$ctx$gthr_tag,
    gistic_significant = ald_inputs$ctx$sig_tag,
    cnv_mode = ald_inputs$ctx$cnv_tag,
    ald_tpm_filter_enabled = ald_inputs$ctx$ald_tpm_filter_enabled,
    ald_min_amplified_median_tpm = if (ald_inputs$ctx$ald_tpm_filter_enabled) ald_min_amplified_median_tpm else NA_real_,
    ald_tpm_tag = ald_inputs$ctx$ald_tpm_tag,
    promoter_regions = paste(promoter_regions, collapse = ";"),
    methylation_platform = platform,
    methylation_sample_type = sample_type,
    min_ald = min_ald,
    min_ref = min_ref,
    delta_methylation_cutoff = delta_methylation_cutoff,
    n_methylation_patients = length(unique(methylation_patients)),
    n_cnv_patients = length(unique(cnv_patients)),
    n_ald_call_patients = length(unique(ald_call_patients)),
    n_methylation_cnv_overlap_patients = length(intersect(methylation_patients, cnv_patients)),
    n_methylation_ald_call_overlap_patients = length(intersect(methylation_patients, ald_call_patients)),
    n_amplified_genes_in_cnmimicry_results = length(ald_inputs$amplified_genes),
    n_ald_eligible_genes_in_cnmimicry_results = length(genes_to_test),
    n_genes_with_promoter_methylation = length(genes_with_promoter),
    n_amplified_genes_with_promoter_methylation = length(intersect(ald_inputs$amplified_genes, genes_with_promoter)),
    n_ald_eligible_genes_with_promoter_methylation = length(genes_testable_by_promoter),
    n_genes_tested_after_min_sample_filters = if (is.null(res_df)) NA_integer_ else nrow(res_df),
    stringsAsFactors = FALSE
  )
}
