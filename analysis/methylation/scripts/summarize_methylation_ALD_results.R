################################################################################
# Purpose
# Summarize ALD promoter methylation results across cohorts.
#
# Inputs
# - Per-cohort ALD promoter methylation Wilcoxon result files produced by
#   analysis/methylation/scripts/run_methylation_ALD_screen.R.
#
# Outputs
# - Cross-cancer ALD promoter methylation summary files.
# - Files are written under analysis/methylation/outputs/fdr_0.05.
################################################################################

source("config.R")
source(file.path("analysis", "methylation", "R", "methylation_helpers.R"))

stopifnot(requireNamespace("dplyr", quietly = TRUE))
stopifnot(requireNamespace("readr", quietly = TRUE))

library(dplyr)
library(readr)

# ------------------------------------------------------------------------------
# Filename tag helpers
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

# ------------------------------------------------------------------------------
# Safe reader
# ------------------------------------------------------------------------------
safe_read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing methylation result: ", path, call. = FALSE)
  
  info <- file.info(path)
  if (is.na(info$size) || info$size == 0) {
    stop("Methylation result is unreadable or has no CSV header: ", path, call. = FALSE)
  }
  
  df <- tryCatch(
    read.csv(path),
    error = function(e) stop("Cannot read methylation result ", path,
                             ": ", conditionMessage(e), call. = FALSE)
  )
  missing_cols <- setdiff(c("gene", "FDR", "delta_methylation"), names(df))
  if (length(missing_cols) > 0) {
    stop("Methylation result is missing columns: ", paste(missing_cols, collapse = ", "),
         " | ", path, call. = FALSE)
  }
  df
}

# ------------------------------------------------------------------------------
# Small helpers
# ------------------------------------------------------------------------------
count_sig <- function(df, fdr_cutoff = 0.05, delta_cutoff = 0) {
  if (nrow(df) == 0 || !"FDR" %in% colnames(df) || !"delta_methylation" %in% colnames(df)) return(0L)
  sum(df$FDR < fdr_cutoff & abs(df$delta_methylation) >= delta_cutoff, na.rm = TRUE)
}

count_sig_hyper <- function(df, fdr_cutoff = 0.05, delta_cutoff = 0) {
  if (nrow(df) == 0 || !"FDR" %in% colnames(df) || !"delta_methylation" %in% colnames(df)) return(0L)
  sum(df$FDR < fdr_cutoff & df$delta_methylation >= delta_cutoff, na.rm = TRUE)
}

count_sig_hypo <- function(df, fdr_cutoff = 0.05, delta_cutoff = 0) {
  if (nrow(df) == 0 || !"FDR" %in% colnames(df) || !"delta_methylation" %in% colnames(df)) return(0L)
  sum(df$FDR < fdr_cutoff & df$delta_methylation <= -delta_cutoff, na.rm = TRUE)
}

get_top_hypo_gene <- function(df, fdr_cutoff = 0.05, delta_cutoff = 0) {
  if (nrow(df) == 0 || !"FDR" %in% colnames(df) || !"delta_methylation" %in% colnames(df)) return(NA_character_)
  
  tmp <- df %>%
    filter(FDR < fdr_cutoff, delta_methylation <= -delta_cutoff) %>%
    arrange(FDR, delta_methylation)
  
  if (nrow(tmp) == 0) return(NA_character_)
  tmp$gene[1]
}

get_top_hypo_delta <- function(df, fdr_cutoff = 0.05, delta_cutoff = 0) {
  if (nrow(df) == 0 || !"FDR" %in% colnames(df) || !"delta_methylation" %in% colnames(df)) return(NA_real_)
  
  tmp <- df %>%
    filter(FDR < fdr_cutoff, delta_methylation <= -delta_cutoff) %>%
    arrange(FDR, delta_methylation)
  
  if (nrow(tmp) == 0) return(NA_real_)
  tmp$delta_methylation[1]
}

# ------------------------------------------------------------------------------
# Run settings
# ------------------------------------------------------------------------------
cohorts <- names(cancers)
prev_thr <- prevalence_thresholds[1]
ald_tpm <- get_methylation_ald_tpm_settings()
fdr_cutoff <- 0.05
promoter_regions <- c("TSS200", "TSS1500", "5'UTR")
methylation_platform <- "Illumina Human Methylation 450"
methylation_sample_type <- "Primary Tumor"
min_ald <- 10
min_ref <- 10
delta_methylation_cutoff <- 0.20

sig_tag <- if (isTRUE(use_gistic_significant)) {
  paste0("with_SIG_", gistic_sig_gene_source)
} else {
  "no_SIG"
}
gthr_tag <- if (isTRUE(use_gistic_thresholded_intersection)) "with_GISTICthr" else "no_GISTICthr"
cnv_tag  <- if (isTRUE(use_cnv_median_normalization)) "cnvnorm" else "rawcnv"
prev_tag <- format(prev_thr, nsmall = 2)
fdr_tag <- paste0("fdr_", format(fdr_cutoff, nsmall = 2))
sample_filter_tag <- paste0("minALD", min_ald, "_minRef", min_ref)
delta_filter_tag <- paste0(
  "delta",
  gsub("\\.", "p", format(delta_methylation_cutoff, nsmall = 2))
)
methylation_setting_tag <- paste(
  make_promoter_region_tag(promoter_regions),
  "HM450",
  "PrimaryTumor",
  sep = "_"
)

# ------------------------------------------------------------------------------
# Summary loop
# ------------------------------------------------------------------------------
summary_list <- lapply(cohorts, function(cohort) {
  
  out_dir <- file.path(workflow_root, "analysis", "methylation", "outputs", cohort)
  
  result_file <- file.path(
    out_dir,
    paste0(
      cohort,
      "_ALD_promoter_methylation_wilcox_",
      gthr_tag, "_",
      sig_tag, "_",
      cnv_tag,
      "_methylation_", methylation_setting_tag,
      "_", sample_filter_tag,
      "_", delta_filter_tag,
      "_prev_", prev_tag,
      ald_tpm$suffix,
      ".csv"
    )
  )
  
  df <- safe_read_csv(result_file)
  
  sig_total <- count_sig(df, fdr_cutoff, delta_methylation_cutoff)
  sig_hyper <- count_sig_hyper(df, fdr_cutoff, delta_methylation_cutoff)
  sig_hypo <- count_sig_hypo(df, fdr_cutoff, delta_methylation_cutoff)
  
  data.frame(
    cohort = cohort,
    genes_tested = nrow(df),
    fdr_cutoff = fdr_cutoff,
    delta_methylation_cutoff = delta_methylation_cutoff,
    sig_genes = sig_total,
    sig_hypermethylated_ALD = sig_hyper,
    sig_hypomethylated_ALD = sig_hypo,
    percent_sig_hypermethylated = ifelse(sig_total > 0, 100 * sig_hyper / sig_total, NA_real_),
    percent_sig_hypomethylated = ifelse(sig_total > 0, 100 * sig_hypo / sig_total, NA_real_),
    top_hypomethylated_gene = get_top_hypo_gene(df, fdr_cutoff, delta_methylation_cutoff),
    top_hypomethylated_delta = get_top_hypo_delta(df, fdr_cutoff, delta_methylation_cutoff),
    stringsAsFactors = FALSE
  )
})

summary_df <- bind_rows(summary_list) %>%
  arrange(desc(sig_genes), cohort)

print(summary_df)

# ------------------------------------------------------------------------------
# Save summary
# ------------------------------------------------------------------------------
summary_out_dir <- file.path(workflow_root, "analysis", "methylation", "outputs", fdr_tag)
dir.create(summary_out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(
  summary_out_dir,
  paste0(
    "ALD_promoter_methylation_summary_",
    gthr_tag, "_",
    sig_tag, "_",
    cnv_tag,
    "_methylation_", methylation_setting_tag,
    "_", sample_filter_tag,
    "_", delta_filter_tag,
    "_prev_", prev_tag,
    ald_tpm$suffix,
    ".csv"
  )
)

write.csv(summary_df, out_file, row.names = FALSE)
message("Saved summary: ", out_file)
