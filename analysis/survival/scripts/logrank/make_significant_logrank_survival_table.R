# Purpose: Combine significant cohort-gene-state results from the log-rank screen.
# Inputs: Per-cohort log-rank results and config.R settings.
# Outputs: Significant-hit CSV/XLSX files under analysis/survival/outputs/logrank,
#   with a genes_in_peak_version subfolder for peak-restricted runs and an aldTPMgt<value>
#   subfolder when TPM filtering is enabled, followed by fdr_<threshold>.

source("config.R")
source(file.path("analysis", "survival", "R", "logrank_helpers.R"))

# ------------------------------------------------------------------------------
# Significant log-rank survival results table
# ------------------------------------------------------------------------------
# Builds one reproducible CSV containing every significant result from the
# log-rank/Kaplan-Meier survival screen under the current default settings.

cohorts_to_run <- names(cancers)
prev_thr <- prevalence_thresholds[1]
ald_tpm <- get_survival_ald_tpm_settings()
fdr_cutoff <- 0.05
min_group <- 10
min_ref <- 10

sample_filter_tag <- paste0("minPos", min_group, "_minRef", min_ref)
sig_tag <- if (isTRUE(use_gistic_significant)) {
  paste0("with_SIG_", gistic_sig_gene_source)
} else {
  "no_SIG"
}
gthr_tag <- if (isTRUE(use_gistic_thresholded_intersection)) "with_GISTICthr" else "no_GISTICthr"
cnv_tag <- if (isTRUE(use_cnv_median_normalization)) "cnvnorm" else "rawcnv"
prev_tag <- format(prev_thr, nsmall = 2)
fdr_tag <- format(fdr_cutoff, nsmall = 2)

logrank_output_root <- resolve_logrank_output_root()
output_dir <- file.path(logrank_output_root, paste0("fdr_", fdr_tag))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_csv <- file.path(
  output_dir,
  paste0(
    "survival_logrank_significant_results_",
    gthr_tag, "_", sig_tag, "_", cnv_tag,
    "_", sample_filter_tag,
    "_prev_", prev_tag,
    "_fdr_", fdr_tag,
    ald_tpm$suffix,
    ".csv"
  )
)
output_xlsx <- sub("\\.csv$", ".xlsx", output_csv)

comparison_defs <- list(
  ALD = list(
    file_part = "ALD_vs_otherDiploids",
    comparison = "ALD vs other diploid",
    group_label = "ALD",
    reference_label = "other diploid",
    n_group_col = "n_ald"
  ),
  DLD = list(
    file_part = "DLD_vs_otherDiploids",
    comparison = "DLD vs other diploid",
    group_label = "DLD",
    reference_label = "other diploid",
    n_group_col = "n_dld"
  ),
  AMP = list(
    file_part = "AMP_vs_otherDiploidsExcludingALD",
    comparison = "AMP vs other diploid excluding ALD",
    group_label = "AMP",
    reference_label = "other diploid excluding ALD",
    n_group_col = "n_amp"
  ),
  DEL = list(
    file_part = "DEL_vs_otherDiploidsExcludingDLD",
    comparison = "DEL vs other diploid excluding DLD",
    group_label = "DEL",
    reference_label = "other diploid excluding DLD",
    n_group_col = "n_del"
  )
)

read_significant_comparison <- function(cohort, state, def) {
  file_name <- paste0(
    cohort, "_", def$file_part, "_OS_logrank_",
    gthr_tag, "_", sig_tag, "_", cnv_tag,
    "_", sample_filter_tag,
    "_prev_", prev_tag,
    ald_tpm$suffix,
    ".csv"
  )
  path <- file.path(logrank_output_root, cohort, file_name)
  if (!file.exists(path)) {
    stop("Missing log-rank output: ", path)
  }
  
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c(
    "gene",
    def$n_group_col,
    "n_other_diploid",
    "n_events_group",
    "n_events_other_diploid",
    "km_median_group",
    "km_median_other_diploid",
    "logrank_chisq",
    "logrank_p_value",
    "logrank_status",
    "logrank_message",
    "FDR"
  )
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing columns in ", path, ": ", paste(missing_cols, collapse = ", "))
  }
  
  df <- df[
    df$logrank_status == "ok" &
      !is.na(df$FDR) &
      df$FDR < fdr_cutoff,
    ,
    drop = FALSE
  ]
  if (nrow(df) == 0) {
    return(NULL)
  }
  
  data.frame(
    cohort = cohort,
    state = state,
    comparison_type = def$comparison,
    gene = df$gene,
    n_group = df[[def$n_group_col]],
    n_reference = df$n_other_diploid,
    n_events_group = df$n_events_group,
    n_events_reference = df$n_events_other_diploid,
    km_median_group = df$km_median_group,
    km_median_reference = df$km_median_other_diploid,
    logrank_chisq = df$logrank_chisq,
    logrank_p_value = df$logrank_p_value,
    FDR = df$FDR,
    logrank_status = df$logrank_status,
    logrank_message = df$logrank_message,
    stringsAsFactors = FALSE
  )
}

all_significant <- do.call(
  rbind,
  unlist(
    lapply(cohorts_to_run, function(cohort) {
      lapply(names(comparison_defs), function(state) {
        read_significant_comparison(cohort, state, comparison_defs[[state]])
      })
    }),
    recursive = FALSE
  )
)

if (is.null(all_significant) || nrow(all_significant) == 0) {
  stop("No significant log-rank results found at FDR < ", fdr_cutoff)
}

all_significant <- all_significant[
  order(
    all_significant$cohort,
    all_significant$state,
    all_significant$FDR,
    all_significant$logrank_p_value
  ),
]

write.csv(all_significant, output_csv, row.names = FALSE)
message("[log-rank significant table] Saved CSV: ", output_csv)
message("[log-rank significant table] Rows: ", nrow(all_significant))

safe_sheet_name <- function(x) {
  x <- gsub("[\\\\/\\?\\*\\[\\]:]", "_", x)
  substr(x, 1, 31)
}

workbook_sheets <- c(
  list(list(name = "all_significant", data = all_significant)),
  lapply(sort(unique(all_significant$cohort)), function(cohort) {
    list(
      name = cohort,
      data = all_significant[all_significant$cohort == cohort, , drop = FALSE]
    )
  })
)

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Package 'openxlsx' is required to write the XLSX output. Install it with install.packages('openxlsx').")
}

wb <- openxlsx::createWorkbook()
header_style <- openxlsx::createStyle(
  fontName = "Arial",
  fontSize = 10,
  textDecoration = "bold",
  fgFill = "#D9EAF7"
)

for (sheet in workbook_sheets) {
  sheet_name <- safe_sheet_name(sheet$name)
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet_name, sheet$data, withFilter = TRUE)
  openxlsx::addStyle(
    wb,
    sheet = sheet_name,
    style = header_style,
    rows = 1,
    cols = seq_along(sheet$data),
    gridExpand = TRUE
  )
  openxlsx::freezePane(wb, sheet = sheet_name, firstRow = TRUE)
  openxlsx::setColWidths(wb, sheet = sheet_name, cols = seq_along(sheet$data), widths = "auto")
}

openxlsx::saveWorkbook(wb, output_xlsx, overwrite = TRUE)
message("[log-rank significant table] Saved XLSX: ", output_xlsx)
message("[log-rank significant table] Workbook sheets: ", length(workbook_sheets))

counts <- aggregate(
  gene ~ cohort + state,
  data = all_significant,
  FUN = length
)
names(counts)[names(counts) == "gene"] <- "significant_results"
print(counts[order(counts$cohort, counts$state), ], row.names = FALSE)
