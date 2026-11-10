################################################################################
# Purpose
# Compare observed ALD/DLD cancer scores with simulated null distributions.
# Scores measure the fraction of eligible genes with mimicry calls in >20% of
# diploid samples. Simulated expression is generated from diploid expression
# means/SDs and multiplied by the CNV matrix used for scoring.
#
# Inputs
# - config.R: cohorts, prevalence_thresholds[1], GISTIC filters/source, CNV
#   normalization/cutoffs, and the optional amplified-median TPM threshold.
# - Matching core summary.csv files for observed cancer_score20_ald/dld scores.
# - Prepared expression/CNV caches and GISTIC gene sets through workflow helpers.
#
# Outputs
# - A permutation summary CSV under analysis/outputs/permutations, containing observed
#   scores, null medians/ranges, empirical p-values, and TPM-filter settings.
# - The filename records CNV/GISTIC settings, prevalence and permutation count;
#   enabled TPM filtering adds aldTPMgt<value> to preserve baseline outputs.
#
# Run settings
# - 1,000 simulations per cohort; random seed 123; score cutoff fixed at >20%.
# - The ALD TPM eligibility rule is reapplied within each simulation; DLD scoring
#   is independent of this filter. Missing/mismatched observed TPM-filter inputs
#   stop the script before the simulation loop.
################################################################################

source("config.R")
source(file.path("R", "get_gistic_significant_sets.R"))

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
ald_tpm_suffix <- if (ald_tpm_filter_enabled) paste0("_", ald_tpm_tag) else ""

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
if (ald_tpm_filter_enabled) results_dir_use <- file.path(results_dir_use, ald_tpm_tag)


if (isTRUE(use_cnv_median_normalization)) {
  cnv_amp_cutoff_use <- cnv_amp_cutoff_norm
  cnv_del_cutoff_use <- cnv_del_cutoff_norm
  dip_center_use     <- cnv_diploid_norm
  dip_tol_use        <- cnv_diploid_tol_norm
} else {
  cnv_amp_cutoff_use <- cnv_amp_cutoff
  cnv_del_cutoff_use <- cnv_del_cutoff
  dip_center_use     <- cnv_diploid
  dip_tol_use        <- 0
}


# ---- packages ----
if (!requireNamespace("readr", quietly = TRUE)) stop("Install readr: install.packages('readr')")
if (!requireNamespace("dplyr", quietly = TRUE)) stop("Install dplyr: install.packages('dplyr')")

library(readr)
library(dplyr)

# ---- load workflow helpers (reuse existing pipeline building blocks) ----
source(file.path("R", "prepare_data.R"))
source(file.path("R", "get_gistic_sets.R"))

# ---- Permutation settings ----
prev_thr <- prevalence_thresholds[1]
score_cutoff_pct <- 20
n_perm <- 1000
set.seed(123)

# ---- helper: read observed score from existing summary.csv ----
read_observed <- function(results_dir, cohort, prev_thr) {
  thr_tag <- format(prev_thr, nsmall = 2)
  f <- file.path(results_dir, cohort, paste0("prev_", thr_tag), "summary.csv")
  if (!file.exists(f)) stop("Missing summary.csv: ", f)
  s <- readr::read_csv(f, show_col_types = FALSE)
  if (ald_tpm_filter_enabled) {
    required <- c("ald_tpm_filter_enabled", "ald_min_amplified_median_tpm", "ald_tpm_tag")
    if (!all(required %in% names(s)) || nrow(s) != 1L ||
        !isTRUE(s$ald_tpm_filter_enabled[1]) ||
        !isTRUE(s$ald_min_amplified_median_tpm[1] == ald_min_amplified_median_tpm) ||
        !isTRUE(s$ald_tpm_tag[1] == ald_tpm_tag)) {
      stop("Observed core summary does not match the configured ALD TPM filter: ", f)
    }
  }
  list(
    obs_ald = s$cancer_score20_ald[1],
    obs_dld = s$cancer_score20_dld[1]
  )
}

#  ---- helper: helper function for normalized perm run ----
make_dip_mask <- function(cnv_mat, center, tol) {
  if (tol == 0) {
    cnv_mat == center
  } else {
    cnv_mat >= (center - tol) & cnv_mat <= (center + tol)
  }
}

# ---- Define eligible gene sets ----
get_event_genes <- function(expr_mat, cnv_mat, gistic_sets, prev_thr,
                            cnv_amp_cutoff, cnv_del_cutoff,
                            use_gistic_thresholded_intersection = TRUE) {
  
  amp_pct <- rowSums(cnv_mat >= cnv_amp_cutoff, na.rm = TRUE) / ncol(cnv_mat)
  amp_genes <- names(amp_pct)[amp_pct > prev_thr]
  if (isTRUE(use_gistic_thresholded_intersection)) {
    if (is.null(gistic_sets)) {
      stop("use_gistic_thresholded_intersection=TRUE but gistic_sets is NULL")
    }
    amp_genes <- intersect(amp_genes, gistic_sets$gistic_amp_genes)
  }
  amp_genes <- intersect(amp_genes, rownames(expr_mat))
  
  del_pct <- rowSums(cnv_mat <= cnv_del_cutoff, na.rm = TRUE) / ncol(cnv_mat)
  del_genes <- names(del_pct)[del_pct > prev_thr]
  if (isTRUE(use_gistic_thresholded_intersection)) {
    if (is.null(gistic_sets)) {
      stop("use_gistic_thresholded_intersection=TRUE but gistic_sets is NULL")
    }
    del_genes <- intersect(del_genes, gistic_sets$gistic_del_genes)
  }
  del_genes <- intersect(del_genes, rownames(expr_mat))
  
  list(amp_genes = amp_genes, del_genes = del_genes)
}

# ---- ALD/DLD cancer-score functions ----
score20_ald <- function(expr_mat, cnv_mat, amp_genes,
                       cnv_amp_cutoff, cnv_diploid, score_cutoff_pct,
                       ald_min_amplified_median_tpm = NULL) {
  if (length(amp_genes) == 0) return(NA_real_)
  
  expr <- expr_mat[amp_genes, , drop = FALSE]
  cnv  <- cnv_mat[amp_genes, , drop = FALSE]
  if (nrow(expr) == 0) return(NA_real_)
  
  # median expression among amplified samples (per gene)
  med_amp <- vapply(seq_len(nrow(expr)), function(i) {
    idx <- which(cnv[i, ] >= cnv_amp_cutoff)
    if (length(idx) == 0) return(NA_real_)
    median(expr[i, idx], na.rm = TRUE)
  }, numeric(1))

  # Reapply the core eligibility rule to each simulated expression matrix.
  if (!is.null(ald_min_amplified_median_tpm)) {
    eligible <- is.finite(med_amp) & med_amp > ald_min_amplified_median_tpm
    expr <- expr[eligible, , drop = FALSE]
    cnv <- cnv[eligible, , drop = FALSE]
    med_amp <- med_amp[eligible]
    if (nrow(expr) == 0) return(NA_real_)
  }
  
  dip_mask <- make_dip_mask(cnv, dip_center_use, dip_tol_use)
  
  # percent of diploid samples above amplified median
  out_count <- vapply(seq_len(nrow(expr)), function(i) {
    if (!is.finite(med_amp[i])) return(0)
    sum(dip_mask[i, ] & (expr[i, ] > med_amp[i]), na.rm = TRUE)
  }, numeric(1))
  
  total_dip <- rowSums(dip_mask, na.rm = TRUE)
  
  pct <- ifelse(total_dip > 0, 100 * out_count / total_dip, NA_real_)
  n_gt <- sum(pct > score_cutoff_pct, na.rm = TRUE)
  
  n_gt / nrow(expr)
}

score20_dld <- function(expr_mat, cnv_mat, del_genes,
                       cnv_del_cutoff, cnv_diploid, score_cutoff_pct) {
  if (length(del_genes) == 0) return(NA_real_)
  
  expr <- expr_mat[del_genes, , drop = FALSE]
  cnv  <- cnv_mat[del_genes, , drop = FALSE]
  if (nrow(expr) == 0) return(NA_real_)
  
  # median expression among deleted samples (per gene)
  med_del <- vapply(seq_len(nrow(expr)), function(i) {
    idx <- which(cnv[i, ] <= cnv_del_cutoff)
    if (length(idx) == 0) return(NA_real_)
    median(expr[i, idx], na.rm = TRUE)
  }, numeric(1))
  
  dip_mask <- make_dip_mask(cnv, dip_center_use, dip_tol_use)
  
  # percent of diploid samples below deleted median
  out_count <- vapply(seq_len(nrow(expr)), function(i) {
    if (!is.finite(med_del[i])) return(0)
    sum(dip_mask[i, ] & (expr[i, ] < med_del[i]), na.rm = TRUE)
  }, numeric(1))
  
  total_dip <- rowSums(dip_mask, na.rm = TRUE)
  
  pct <- ifelse(total_dip > 0, 100 * out_count / total_dip, NA_real_)
  n_gt <- sum(pct > score_cutoff_pct, na.rm = TRUE)
  
  n_gt / length(del_genes)
}


# ---- permutation (simulate diploid expression, then multiply by raw CNV) ----
# For each gene:
#   - estimate mean/sd from DIPLOID samples (CN == cnv_diploid)
#   - simulate expression for all samples from N(mean_dip, sd_dip)
#   - multiply by RAW CNV 

compute_diploid_params <- function(expr_mat, cnv_mat, cnv_diploid) {
  dip_mask <- make_dip_mask(cnv_mat, dip_center_use, dip_tol_use)
  
  n_dip <- rowSums(dip_mask, na.rm = TRUE)
  
  # mean of diploid expression per gene
  sum_dip <- rowSums(expr_mat * dip_mask, na.rm = TRUE)
  mean_dip <- ifelse(n_dip > 0, sum_dip / n_dip, NA_real_)
  
  # variance via E[x^2] - (E[x])^2, using only diploid samples
  sumsq_dip <- rowSums((expr_mat^2) * dip_mask, na.rm = TRUE)
  # sample variance: (sumsq - n*mean^2) / (n-1)
  var_dip <- ifelse(n_dip > 1, (sumsq_dip - n_dip * (mean_dip^2)) / (n_dip - 1), 0)
  
  # numerical guard
  var_dip[!is.finite(var_dip) | var_dip < 0] <- 0
  sd_dip <- sqrt(var_dip)
  
  # fallback for genes with no diploid samples: use overall mean, sd=0
  overall_mean <- rowMeans(expr_mat, na.rm = TRUE)
  mean_dip[!is.finite(mean_dip)] <- overall_mean[!is.finite(mean_dip)]
  sd_dip[!is.finite(sd_dip)] <- 0
  
  list(mean_dip = mean_dip, sd_dip = sd_dip)
}

simulate_expr_cnv <- function(mean_dip, sd_dip, cnv_mat) {
  n_genes <- length(mean_dip)
  n_samp  <- ncol(cnv_mat)
  
  # simulate N(mean_dip, sd_dip) per gene across samples
  z <- matrix(rnorm(n_genes * n_samp), nrow = n_genes, ncol = n_samp)
  sim_expr <- sweep(z, 1, sd_dip, `*`)
  sim_expr <- sweep(sim_expr, 1, mean_dip, `+`)
  
  # multiply by CNV matrix used for scoring (raw or normalized)
  sim_expr * cnv_mat
}


# ---- output ----
# Validate every observed input before starting the long permutation run.
observed_scores <- lapply(cancers, function(cdef) {
  read_observed(results_dir_use, cdef$cohort, prev_thr)
})

perm_dir <- file.path(workflow_root, "analysis", "outputs", "permutations")
dir.create(perm_dir, showWarnings = FALSE, recursive = TRUE)

out_rows <- list()

for (cancer_name in names(cancers)) {
  cdef <- cancers[[cancer_name]]
  cohort <- cdef$cohort
  
  message("\n--- ", cohort, " ---")
  
  obs <- observed_scores[[cancer_name]]
  
  prepared <- prepare_data(
    cohort = cohort,
    gdc_project = cdef$gdc_project,
    gdc_dir = cdef$gdc_dir,
    cache_dir = cache_dir,
    er_pos_only = isTRUE(cdef$er_pos_only)
  )
  
  expr_mat <- prepared$expr_mat
  cnv_mat  <- prepared$cnv_mat
  
  if (isTRUE(use_cnv_median_normalization)) {
    denom <- apply(cnv_mat, 2, median, na.rm = TRUE)
    cnv_mat <- sweep(cnv_mat, 2, denom, "/")
  }
  
  gistic_sets <- NULL
  if (isTRUE(use_gistic_thresholded_intersection)) {
    gistic_sets <- get_gistic_sets(cohort = cohort, cache_dir = cache_dir, prev_thr = prev_thr)
  }
  
  genes <- get_event_genes(
    expr_mat = expr_mat,
    cnv_mat = cnv_mat,
    gistic_sets = gistic_sets,
    prev_thr = prev_thr,
    cnv_amp_cutoff = cnv_amp_cutoff_use,
    cnv_del_cutoff = cnv_del_cutoff_use,
    use_gistic_thresholded_intersection = isTRUE(use_gistic_thresholded_intersection)
  )
  
  if (isTRUE(use_gistic_significant)) {
    gistic_sig <- get_gistic_significant_sets(
      cohort = cohort,
      cache_dir = cache_dir,
      gdac_resources_dir = gdac_resources_dir,
      gdac_run_id = gdac_run_id,
      gistic_sig_q_cutoff = gistic_sig_q_cutoff,
      gistic_sig_gene_source = gistic_sig_gene_source
    )
    genes$amp_genes <- intersect(genes$amp_genes, gistic_sig$sig_amp_genes)
    genes$del_genes <- intersect(genes$del_genes, gistic_sig$sig_del_genes)
  }
  
  
  dip_params <- compute_diploid_params(
    expr_mat = expr_mat,
    cnv_mat  = cnv_mat,
    cnv_diploid = cnv_diploid
  )
  
  perm_ald <- numeric(n_perm)
  perm_dld <- numeric(n_perm)
  
  for (b in seq_len(n_perm)) {
    expr_perm <- simulate_expr_cnv(
      mean_dip = dip_params$mean_dip,
      sd_dip   = dip_params$sd_dip,
      cnv_mat  = cnv_mat
    )
    
    perm_ald[b] <- score20_ald(expr_perm, cnv_mat, genes$amp_genes,
                             cnv_amp_cutoff_use, cnv_diploid, score_cutoff_pct,
                             ald_min_amplified_median_tpm = ald_min_amplified_median_tpm)
    
    perm_dld[b] <- score20_dld(expr_perm, cnv_mat, genes$del_genes,
                             cnv_del_cutoff_use, cnv_diploid, score_cutoff_pct)
  }
  
  # An unavailable observed score or empty null distribution is not a test.
  has_ald <- any(!is.na(perm_ald))
  has_dld <- any(!is.na(perm_dld))
  p_ald <- if (!is.na(obs$obs_ald) && has_ald) {
    (sum(perm_ald >= obs$obs_ald, na.rm = TRUE) + 1) / (sum(!is.na(perm_ald)) + 1)
  } else NA_real_
  p_dld <- if (!is.na(obs$obs_dld) && has_dld) {
    (sum(perm_dld >= obs$obs_dld, na.rm = TRUE) + 1) / (sum(!is.na(perm_dld)) + 1)
  } else NA_real_
  
  out_rows[[cohort]] <- data.frame(
    cancer = cohort,
    prev_thr = prev_thr,
    n_perm = n_perm,
    ald_tpm_filter_enabled = ald_tpm_filter_enabled,
    ald_min_amplified_median_tpm = if (ald_tpm_filter_enabled) ald_min_amplified_median_tpm else NA_real_,
    ald_tpm_tag = ald_tpm_tag,
    
    observed_score20_ald = obs$obs_ald,
    perm_median_score20_ald = median(perm_ald, na.rm = TRUE),
    perm_min_score20_ald = if (has_ald) min(perm_ald, na.rm = TRUE) else NA_real_,
    perm_max_score20_ald = if (has_ald) max(perm_ald, na.rm = TRUE) else NA_real_,
    p_value_ald = p_ald,
    
    observed_score20_dld = obs$obs_dld,
    perm_median_score20_dld = median(perm_dld, na.rm = TRUE),
    perm_min_score20_dld = if (has_dld) min(perm_dld, na.rm = TRUE) else NA_real_,
    perm_max_score20_dld = if (has_dld) max(perm_dld, na.rm = TRUE) else NA_real_,
    p_value_dld = p_dld,
    
    min_p_value = if (all(is.na(c(p_ald, p_dld)))) NA_real_ else min(p_ald, p_dld, na.rm = TRUE),
    
    stringsAsFactors = FALSE
  )
}

perm_summary <- dplyr::bind_rows(out_rows)

gthr_tag <- if (isTRUE(use_gistic_thresholded_intersection)) "withGISTICthr" else "noGISTICthr"
sig_tag  <- if (isTRUE(use_gistic_significant)) {
  paste0("withSIG_", gistic_sig_gene_source)
} else {
  "noSIG"
}
cnv_tag  <- if (isTRUE(use_cnv_median_normalization)) "cnvnorm" else "rawcnv"

tag <- paste(cnv_tag, gthr_tag, sig_tag, sep = "_")
out_file <- file.path(perm_dir, paste0("perm_summary_", tag, "_prev_", format(prev_thr, nsmall = 2), "_n", n_perm, ald_tpm_suffix, ".csv"))
readr::write_csv(perm_summary, out_file)

message("\nSaved: ", out_file)
