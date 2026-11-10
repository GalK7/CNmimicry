# Purpose:
#   Define shared settings for the core workflow, downstream analyses and
#   figure/table builders.
#
# Settings:
#   - Project-relative paths, cohort selection and BRCA ER-positive restriction.
#   - GISTIC filters, resource version and significant-peak gene source.
#   - CNV normalization, copy-number thresholds and prevalence thresholds.
#   - Minimum amplified-group median TPM for ALD classification and the
#     sample-level burden cutoff.
#
# Usage:
#   Sourced by workflow scripts; this file does not generate analysis outputs.

# ---- Project paths ----
get_workflow_root <- function() {
  frame_files <- vapply(
    sys.frames(),
    function(frame) {
      if (is.null(frame$ofile)) {
        NA_character_
      } else {
        frame$ofile
      }
    },
    character(1)
  )
  
  config_candidates <- frame_files[
    !is.na(frame_files) & grepl("(^|[/\\\\])config\\.R$", frame_files)
  ]
  
  if (length(config_candidates) > 0) {
    config_file <- normalizePath(tail(config_candidates, 1), winslash = "/", mustWork = TRUE)
    return(dirname(config_file))
  }
  
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

workflow_root <- get_workflow_root()
cache_dir     <- file.path(workflow_root, "cache")
gdc_data_dir  <- file.path(workflow_root, "GDCdata")
results_dir   <- file.path(workflow_root, "results")
results_dir_sig  <- file.path(workflow_root, "results_sig")

# ---- Results paths for no-thresholded-GISTIC-intersection test ----
results_dir_no_gistic_thr <- file.path(workflow_root, "results_noGISTICthr")
results_dir_sig_no_gistic_thr <- file.path(workflow_root, "results_sig_noGISTICthr")
results_dir_cnvnorm_no_gistic_thr <- file.path(workflow_root, "results_cnvnorm_noGISTICthr")

# ---- Optional intersection with thresholded GISTIC2 calls ----
use_gistic_thresholded_intersection <- FALSE


# ---- Additional optional filter: GDAC GISTIC2 significant genes ----
# Parses GDAC "CopyNumber_Gistic2.Level_4" outputs (amp_genes.conf_99.txt / del_genes.conf_99.txt)
# and extracts genes from significant peaks using either the current wide-peak definition
# or the stricter genes_in_peak definition from table_amp/table_del.

use_gistic_significant <- TRUE         # toggle: TRUE/FALSE
gdac_run_id <- "analyses__2016_01_28"    # reproducible GDAC run
gistic_sig_q_cutoff <- 0.25             # peak q-value cutoff
gistic_sig_gene_source <- "wide_peak"    # "wide_peak" or "genes_in_peak"
gdac_resources_dir <- file.path(workflow_root, "resources", "gdac_gistic_level4")


# ---- Optional test: CNV median-normalization (ASCAT3) ----
use_cnv_median_normalization <- TRUE   # toggle: TRUE/FALSE
cnv_norm_stat <- "median"               # keep as string for future ("mean")
results_dir_cnvnorm <- file.path(workflow_root, "results_cnvnorm")


cnv_amp_cutoff_norm <- 2
cnv_del_cutoff_norm <- 0.67
cnv_diploid_norm <- 1
cnv_diploid_tol_norm <- 0.05   # diploid window: 1 ± 0.05  => [0.95, 1.05]


# ---- Cancers to run (EXACT Excel order) ----
make_cancer <- function(cohort, er_pos_only = FALSE) {
  out <- list(
    cohort = cohort,
    gdc_project = paste0("TCGA-", cohort),
    gdc_dir = gdc_data_dir
  )
  
  if (isTRUE(er_pos_only)) {
    out$er_pos_only <- TRUE
  }
  
  out
}

cancers <- list(
  BRCA = make_cancer("BRCA", er_pos_only = TRUE),
  LIHC = make_cancer("LIHC"),
  LUAD = make_cancer("LUAD"),
  LUSC = make_cancer("LUSC"),
  SKCM = make_cancer("SKCM"),
  COAD = make_cancer("COAD"),
  BLCA = make_cancer("BLCA"),
  SARC = make_cancer("SARC"),
  OV   = make_cancer("OV"),
  ACC  = make_cancer("ACC"),
  CESC = make_cancer("CESC"),
  CHOL = make_cancer("CHOL"),
  ESCA = make_cancer("ESCA"),
  GBM  = make_cancer("GBM"),
  HNSC = make_cancer("HNSC"),
  STAD = make_cancer("STAD"),
  LGG  = make_cancer("LGG"),
  UCEC = make_cancer("UCEC"),
  PCPG = make_cancer("PCPG"),
  LAML = make_cancer("LAML")
)

# ---- Thresholds ----
prevalence_thresholds <- c(0.03)

# ---- ALD amplified-expression eligibility ----
# NULL disables this filter; otherwise require median amplified TPM strictly > this value.
# Filtered core results are written separately under aldTPMgt<value>.
ald_min_amplified_median_tpm <- 1

# ---- Raw ASCAT3 CNV definitions (used only when use_cnv_median_normalization = FALSE) ----
cnv_amp_cutoff <- 7
cnv_del_cutoff <- 1
cnv_diploid    <- 2

# ---- Sample-level cutoff ----
sample_fraction_cutoff <- 0.05   # strictly >
