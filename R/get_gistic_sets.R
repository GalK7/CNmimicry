# Purpose:
#   Build amplification/deletion gene sets from thresholded GISTIC calls for
#   the optional intersection filter, disabled in the default configuration.
#
# Inputs:
#   - Cohort, cache path, amplification prevalence and deletion thresholds.
#   - Cached gene sets or thresholded calls retrieved with TCGAbiolinks::getGistic().
#
# Outputs:
#   - A list of amplification/deletion gene sets and their threshold settings.
#   - Parsed gene-set caches under cache/gistic_sets and downloaded resources
#     under cache/gistic, relative to the supplied cache_dir.
#
# Usage:
#   Defines get_gistic_sets(), called by run_core() when the filter is enabled.

get_gistic_sets <- function(cohort,
                            cache_dir,
                            prev_thr,
                            del_base_cutoff = 0.10,
                            del_focal_cutoff = 0.01) {
  # cohort: "GBM", "BRCA", etc. (NO "TCGA-")
  # prev_thr: 0.03 / 0.05 / 0.10  (MUST match ASCAT3 prevalence threshold for this run)
  # cache_dir: project-local cache directory, e.g. file.path(workflow_root, "cache")
  
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    stop("Package 'TCGAbiolinks' is required.")
  }
  
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  
  thr_tag <- format(prev_thr, nsmall = 2)
  gistic_sets_cache_dir <- file.path(cache_dir, "gistic_sets")
  dir.create(gistic_sets_cache_dir, showWarnings = FALSE, recursive = TRUE)
  cache_file <- file.path(gistic_sets_cache_dir, paste0(cohort, "_gistic_sets_prev_", thr_tag, ".rds"))
  
  if (file.exists(cache_file)) {
    message("[get_gistic_sets] Using cached GISTIC sets: ", cache_file)
    return(readRDS(cache_file))
  }
  
  message("[get_gistic_sets] Building GISTIC sets for: ", cohort, " | prev_thr=", prev_thr)
  
  # ---- Run getGistic() inside cache/gistic to avoid cluttering Workflow/ ----
  gistic_cache_dir <- file.path(cache_dir, "gistic")
  dir.create(gistic_cache_dir, showWarnings = FALSE, recursive = TRUE)
  
  old_wd <- getwd()
  setwd(gistic_cache_dir)
  on.exit(setwd(old_wd), add = TRUE)
  
  # ---- getGistic() ----
  gistic_dir <- TCGAbiolinks::getGistic(cohort)
  
  # Sample columns (TCGA barcodes)
  gistic_samp <- grep("^TCGA-", colnames(gistic_dir))
  if (length(gistic_samp) == 0) stop("No TCGA sample columns found in getGistic() output.")
  
  # Numeric GISTIC matrix {-2,-1,0,1,2}
  gistic_mat <- as.matrix(
    suppressWarnings(apply(gistic_dir[, gistic_samp, drop = FALSE], 2, as.numeric))
  )
  rownames(gistic_mat) <- seq_len(nrow(gistic_dir))
  
  # Gene symbol column
  gistic_gene_col <- "Gene Symbol"
  if (!gistic_gene_col %in% colnames(gistic_dir)) {
    stop("Expected gene column 'Gene Symbol' not found in getGistic() output.")
  }
  
  # ---- Thresholded sets (match ASCAT3 prevalence threshold) ----
  # Amplification set: fraction of samples with GISTIC == 2 > prev_thr
  amp_mask <- rowMeans(gistic_mat == 2, na.rm = TRUE) > prev_thr
  
  # Deletion gene-set criteria:
  # base: (<= -1) in >=10% of samples
  # focal: (-2) in >=1% of samples
  del_base  <- rowMeans(gistic_mat <= -1, na.rm = TRUE) >= del_base_cutoff
  del_focal <- rowMeans(gistic_mat == -2, na.rm = TRUE) >= del_focal_cutoff
  del_mask  <- del_base & del_focal
  
  gistic_amp_genes <- unique(gistic_dir[[gistic_gene_col]][amp_mask])
  gistic_del_genes <- unique(gistic_dir[[gistic_gene_col]][del_mask])
  
  # Clean
  gistic_amp_genes <- gistic_amp_genes[!is.na(gistic_amp_genes) & nzchar(gistic_amp_genes)]
  gistic_del_genes <- gistic_del_genes[!is.na(gistic_del_genes) & nzchar(gistic_del_genes)]
  
  out <- list(
    cohort = cohort,
    prev_thr = prev_thr,
    del_base_cutoff = del_base_cutoff,
    del_focal_cutoff = del_focal_cutoff,
    gistic_amp_genes = gistic_amp_genes,
    gistic_del_genes = gistic_del_genes
  )
  
  saveRDS(out, cache_file)
  message("[get_gistic_sets] Saved cache: ", cache_file)
  
  out
}
