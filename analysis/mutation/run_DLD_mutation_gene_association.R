# Purpose:
#   Test gene-level associations between DLD calls and somatic mutation status
#   across cohorts.
#
# Inputs:
#   - Core workflow DLD call files from result directories defined in config.R.
#   - Cached or downloaded TCGA somatic mutation data.
#
# Outputs:
#   - Mutation association workbooks written to analysis/outputs/mutation.

source("config.R")

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

stopifnot(requireNamespace("dplyr", quietly = TRUE))
stopifnot(requireNamespace("readr", quietly = TRUE))
stopifnot(requireNamespace("purrr", quietly = TRUE))
stopifnot(requireNamespace("tibble", quietly = TRUE))
stopifnot(requireNamespace("TCGAbiolinks", quietly = TRUE))
stopifnot(requireNamespace("openxlsx", quietly = TRUE))


library(dplyr)


# ----------------------------
# Mutation classes used for the optional mutation-set filters.
# ----------------------------
nonsilent_classes <- c(
  "Missense_Mutation",
  "Nonsense_Mutation",
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "Splice_Site",
  "In_Frame_Del",
  "In_Frame_Ins"
)

silent_classes <- c("Silent")
# ----------------------------
# MAF caching (Masked Somatic Mutation + Ensemble Merging/Masking)
# ----------------------------
get_maf_cached_ensemble <- function(gdc_project, cache_dir, gdc_dir = NULL) {
  maf_dir <- file.path(workflow_root, "analysis", "cache", "maf_ensemble_masked")
  maf_gdc_dir <- file.path(workflow_root, "analysis", "GDCdata")
  dir.create(maf_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(maf_gdc_dir, recursive = TRUE, showWarnings = FALSE)
  
  maf_path <- file.path(maf_dir, paste0(gdc_project, "_ensemble_masked_maf.rds"))
  if (file.exists(maf_path)) return(readRDS(maf_path))
  
  query_maf <- TCGAbiolinks::GDCquery(
    project       = gdc_project,
    data.category = "Simple Nucleotide Variation",
    data.type     = "Masked Somatic Mutation",
    access        = "open",
    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
  )
  
  TCGAbiolinks::GDCdownload(query_maf, directory = maf_gdc_dir)
  maf_data <- TCGAbiolinks::GDCprepare(query_maf, directory = maf_gdc_dir)
  
  saveRDS(maf_data, maf_path)
  maf_data
}

# ----------------------------
# Fisher tests per gene: DLD(gene) vs MUT(gene)
# ----------------------------
empty_mutation_association <- function() {
  tibble::tibble(
    Gene = character(), n_dld = integer(), n_mut = integer(),
    a_mut_and_dld = integer(), b_dld_only = integer(),
    c_mut_only = integer(), d_neither = integer(),
    fisher_p_two = numeric(), fisher_p_cooc = numeric(), fisher_p_mutexc = numeric(),
    fdr_two = numeric(), fdr_cooc = numeric(), fdr_mutexc = numeric()
  )
}

run_assoc_one <- function(cdef, prev_thr, results_dir_use,
                          cache_dir, mutation_set = "all_mutations",
                          min_mut = 5L, min_dld = 5L) {
  
  cohort <- cdef$cohort
  gdc_project <- cdef$gdc_project
  gdc_dir <- cdef$gdc_dir
  
  out_dir <- file.path(results_dir_use, cohort, paste0("prev_", prev_thr))
  
  dld_calls_path <- file.path(out_dir, "dld_calls_long.rds")
  dld_ps_path    <- file.path(out_dir, "dld_per_sample.rds")
  
  if (!file.exists(dld_calls_path)) stop("Missing: ", dld_calls_path)
  if (!file.exists(dld_ps_path))    stop("Missing: ", dld_ps_path)
  
  # Universe of patients = those in dld_per_sample (keeps consistent sample space)
  dld_per_sample <- readRDS(dld_ps_path)
  universe_patients <- dld_per_sample |>
    transmute(patient = substr(Sample, 1, 12)) |>
    distinct() |>
    pull(patient)
  
  N <- length(universe_patients)
  if (N == 0) {
    message(cohort, ": no patients available for testing.")
    return(empty_mutation_association())
  }
  
  # DLD calls long: (patient, Gene)
  dld_calls <- readRDS(dld_calls_path) |>
    transmute(patient = substr(Sample, 1, 12), Gene) |>
    distinct()
  
  # If no DLD calls in this setting, still write an empty CSV so downstream is clean.
  if (nrow(dld_calls) == 0) {
    message(cohort, ": no DLD calls available for testing.")
    out <- empty_mutation_association()
    readr::write_csv(out, file.path(out_dir, paste0(
      "mutation_x_DLD_gene_association_", mutation_set,
      "_minMut", min_mut, "_minDLD", min_dld, ".csv"
    )))
    return(out)
  }
  
  # MAF (patient, Gene), optionally restricted by mutation class.
  maf_data <- get_maf_cached_ensemble(gdc_project, cache_dir, gdc_dir)

  maf_min <- maf_data |>
    dplyr::transmute(
      patient = substr(Tumor_Sample_Barcode, 1, 12),
      Gene    = Hugo_Symbol,
      Variant_Classification = Variant_Classification
    ) |>
    dplyr::filter(
      dplyr::case_when(
        mutation_set == "all_mutations" ~ TRUE,
        mutation_set == "non_silent" ~ Variant_Classification %in% nonsilent_classes,
        mutation_set == "silent" ~ Variant_Classification %in% silent_classes,
        TRUE ~ FALSE
      )
    ) |>
    dplyr::select(patient, Gene) |>
    dplyr::distinct() |>
    dplyr::filter(patient %in% universe_patients)
  
  # Counts per gene
  dld_counts <- dld_calls |>
    count(Gene, name = "n_dld") |>
    filter(n_dld >= min_dld)
  
  mut_counts <- maf_min |>
    count(Gene, name = "n_mut") |>
    filter(n_mut >= min_mut)
  
  genes <- intersect(dld_counts$Gene, mut_counts$Gene)
  if (length(genes) == 0) {
    message(cohort, ": no genes meet the mutation and DLD count thresholds (", mutation_set, ").")
    return(empty_mutation_association())
  }
  
  dld_g  <- dld_calls |> filter(Gene %in% genes)
  maf_g <- maf_min  |> filter(Gene %in% genes)
  
  # a = mutated AND DLD for same gene
  a_counts <- dld_g |>
    inner_join(maf_g, by = c("patient", "Gene")) |>
    count(Gene, name = "a_mut_and_dld")
  
  stats <- tibble::tibble(Gene = genes) |>
    left_join(dld_counts, by = "Gene") |>
    left_join(mut_counts, by = "Gene") |>
    left_join(a_counts,  by = "Gene") |>
    mutate(
      a_mut_and_dld = dplyr::coalesce(a_mut_and_dld, 0L),
      b_dld_only    = n_dld  - a_mut_and_dld,
      c_mut_only   = n_mut - a_mut_and_dld,
      d_neither    = N - (a_mut_and_dld + b_dld_only + c_mut_only)
    )
  
  # Fisher tests per gene: two-sided + two one-sided (co-occurrence / mutual exclusivity)
  fisher_p_two    <- numeric(nrow(stats))
  fisher_p_cooc   <- numeric(nrow(stats))  # alternative = "greater"
  fisher_p_mutexc <- numeric(nrow(stats))  # alternative = "less"
  
  for (i in seq_len(nrow(stats))) {
    a <- stats$a_mut_and_dld[i]
    b <- stats$b_dld_only[i]
    c <- stats$c_mut_only[i]
    d <- stats$d_neither[i]
    
    tab <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    
    # two-sided: any association at all
    fisher_p_two[i]    <- fisher.test(tab, alternative = "two.sided")$p.value
    # co-occurrence: overlap bigger than expected
    fisher_p_cooc[i]   <- fisher.test(tab, alternative = "greater")  $p.value
    # mutual exclusivity: overlap smaller than expected
    fisher_p_mutexc[i] <- fisher.test(tab, alternative = "less")     $p.value
  }
  
  out <- stats |>
    rename(
      n_dld = n_dld,
      a_mut_and_dld = a_mut_and_dld,
      b_dld_only = b_dld_only
    ) |>
    mutate(
      fisher_p_two    = fisher_p_two,
      fisher_p_cooc   = fisher_p_cooc,
      fisher_p_mutexc = fisher_p_mutexc,
      fdr_two         = p.adjust(fisher_p_two,    method = "BH"),
      fdr_cooc        = p.adjust(fisher_p_cooc,   method = "BH"),
      fdr_mutexc      = p.adjust(fisher_p_mutexc, method = "BH")
    ) |>
    arrange(fdr_cooc, fdr_mutexc, fdr_two)
  
  out
}

# ----------------------------
# Match run_all.R results_dir logic exactly
# ----------------------------
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

# ----------------------------
# Optionally filter cancers here (e.g., only LGG)
# cancers_to_run <- cancers[names(cancers) %in% c("LGG")]
# ----------------------------
cancers_to_run <- cancers


# ----------------------------
# Output directory for Excel summaries
# ----------------------------
output_dir <- file.path(workflow_root, "analysis", "outputs", "mutation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ----------------------------
# Run across mutation sets × thresholds × cancers and write one Excel per setting
# ----------------------------
mutation_sets_to_run <- c("all_mutations", "non_silent", "silent")

for (mutation_set in mutation_sets_to_run) {
  for (thr in prevalence_thresholds) {
    message("\n[mutation_x_DLD] Building Excel for mutation_set=", mutation_set, " | prev_thr=", thr)
    
    wb <- openxlsx::createWorkbook()
    
    for (cancer_name in names(cancers_to_run)) {
      cdef <- cancers_to_run[[cancer_name]]
      message("[mutation_x_DLD]  Running: ", cancer_name, " | mutation_set=", mutation_set, " | prev_thr=", thr)
      
      res <- tryCatch(
        run_assoc_one(
          cdef = cdef,
          prev_thr = thr,
          results_dir_use = results_dir_use,
          cache_dir = cache_dir,
          mutation_set = mutation_set,
          min_mut = 5L,
          min_dld = 5L
        ),
        error = function(e) {
          stop("Mutation association failed for ", cancer_name,
               " | mutation_set=", mutation_set, " | prev_thr=", thr,
               ": ", conditionMessage(e), call. = FALSE)
        }
      )
      
      openxlsx::addWorksheet(wb, sheetName = cancer_name)
      openxlsx::writeData(wb, sheet = cancer_name, x = res)
    }
    
    thr_str  <- as.character(thr)
    
    gthr_tag <- if (isTRUE(use_gistic_thresholded_intersection)) "with_GISTICthr" else "no_GISTICthr"
    sig_tag  <- gistic_sig_folder_tag
    cnv_tag  <- if (isTRUE(use_cnv_median_normalization)) "CNVnorm" else "rawCNV"
    
    out_xlsx <- file.path(
      output_dir,
      paste0(
        "mutation_x_DLD_gene_association_", thr_str, "_",
        mutation_set, "_", gthr_tag, "_", sig_tag, "_", cnv_tag, ".xlsx"
      )
    )
    openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
    message("[mutation_x_DLD] Saved: ", out_xlsx)
  }
}

message("\nDone. Wrote one Excel per mutation set and threshold (sheets = cancers) into analysis/outputs/mutation.")
