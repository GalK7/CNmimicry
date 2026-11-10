# Purpose:
#   Run the core ALD/DLD workflow for the cohorts and prevalence thresholds
#   selected in config.R.
#
# Inputs:
#   - Project settings from config.R and helper functions under R/.
#   - Prepared expression/CNV data and GISTIC resources, cached or retrieved
#     as required by the selected settings.
#
# Outputs:
#   - Gene sets, ALD/DLD calls, sample-level burdens and cohort summaries
#     in the configured results directories.
#
# Usage:
#   From the Workflow root: source("run_all.R")

# 1) Load config (paths, cancers, thresholds)
source("config.R")

# 2) Load functions
source(file.path("R", "prepare_data.R"))
source(file.path("R", "get_gistic_sets.R"))
source(file.path("R", "run_core.R"))
source(file.path("R", "get_gistic_significant_sets.R"))

# 3) Make sure output dirs exist
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

if (isTRUE(use_gistic_thresholded_intersection)) {
  results_dir_base <- results_dir
  results_dir_sig_base <- results_dir_sig
  results_dir_cnvnorm_base <- results_dir_cnvnorm
} else {
  results_dir_base <- results_dir_no_gistic_thr
  results_dir_sig_base <- results_dir_sig_no_gistic_thr
  results_dir_cnvnorm_base <- results_dir_cnvnorm_no_gistic_thr
}

gistic_sig_folder_tag <- if (isTRUE(use_gistic_significant)) {
  paste0("with_SIG_", gistic_sig_gene_source)
} else {
  "no_SIG"
}

results_dir_use <- if (isTRUE(use_cnv_median_normalization)) {
  file.path(results_dir_cnvnorm_base, gistic_sig_folder_tag)
} else {
  if (isTRUE(use_gistic_significant)) results_dir_sig_base else results_dir_base
}

dir.create(results_dir_use, showWarnings = FALSE, recursive = TRUE)

dir.create(gdac_resources_dir, showWarnings = FALSE, recursive = TRUE)

# 4) Run all cancers x thresholds
all_summaries <- list()

for (cancer_name in names(cancers)) {
  cdef <- cancers[[cancer_name]]
  
  for (thr in prevalence_thresholds) {
    message("\n==============================")
    message("Running: ", cdef$cohort, " | prev_thr = ", thr)
    message("==============================\n")
    
    s <- run_core(
      cohort = cdef$cohort,
      gdc_project = cdef$gdc_project,
      gdc_dir = cdef$gdc_dir,
      prev_thr = thr,
      cache_dir = cache_dir,
      results_dir = results_dir_use,
      cnv_amp_cutoff = cnv_amp_cutoff,
      cnv_del_cutoff = cnv_del_cutoff,
      cnv_diploid = cnv_diploid,
      sample_fraction_cutoff = sample_fraction_cutoff,
      er_pos_only = isTRUE(cdef$er_pos_only),
      use_cnv_median_normalization = isTRUE(use_cnv_median_normalization),
      cnv_norm_stat = cnv_norm_stat,
      cnv_amp_cutoff_norm = cnv_amp_cutoff_norm,
      cnv_del_cutoff_norm = cnv_del_cutoff_norm,
      cnv_diploid_norm = cnv_diploid_norm,
      cnv_diploid_tol_norm = cnv_diploid_tol_norm,
      use_gistic_thresholded_intersection = isTRUE(use_gistic_thresholded_intersection),
      use_gistic_significant = isTRUE(use_gistic_significant),
      gdac_run_id = gdac_run_id,
      gistic_sig_q_cutoff = gistic_sig_q_cutoff,
      gistic_sig_gene_source = gistic_sig_gene_source,
      gdac_resources_dir = gdac_resources_dir,
      ald_min_amplified_median_tpm = ald_min_amplified_median_tpm
    )
    
    all_summaries[[paste0(cdef$cohort, "_", format(thr, nsmall = 2))]] <- s
  }
}

message("\nDONE.")
