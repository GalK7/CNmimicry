# Purpose:
#   Retrieve genes from significant GDAC GISTIC2 Level_4 peaks using the
#   wide_peak definition or the narrower genes_in_peak definition.
#
# Inputs:
#   - Cohort, GDAC run, peak q-value cutoff, gene source and resource/cache paths.
#   - Cached gene sets or downloaded GDAC peak files.
#
# Outputs:
#   - Amplification/deletion gene sets with source metadata.
#   - Parsed RDS caches under cache_dir/gistic_significant and downloaded
#     resources with a download manifest under gdac_resources_dir.
#
# Usage:
#   Defines get_gistic_significant_sets(), called by run_core() when significant
#   peak filtering is enabled. Existing caches are reused unless redownload
#   is explicitly requested.

get_gistic_significant_sets <- function(
    cohort,
    cache_dir,
    gdac_resources_dir,
    gdac_run_id = "analyses__2016_01_28",
    gistic_sig_q_cutoff = 0.25,
    gistic_sig_gene_source = "wide_peak",
    cohort_suffix = "-TP",
    force_redownload = FALSE
) {
  allowed_gene_sources <- c("wide_peak", "genes_in_peak")
  if (!gistic_sig_gene_source %in% allowed_gene_sources) {
    stop(
      "gistic_sig_gene_source must be one of: ",
      paste(allowed_gene_sources, collapse = ", ")
    )
  }
  
  update_gistic_manifest <- function(
      manifest_file,
      cohort,
      data_type,
      gdac_run_id,
      gistic_module,
      source_base_url,
      file_name,
      tar_path,
      notes
  ) {
    file_meta <- file.info(tar_path)
    
    new_row <- data.frame(
      cohort = cohort,
      data_type = data_type,
      gdac_run = gdac_run_id,
      gistic_module = gistic_module,
      source_base_url = source_base_url,
      file_name = file_name,
      file_size_bytes = unname(file_meta$size),
      file_created_time = format(file_meta$ctime, "%Y-%m-%d %H:%M:%S"),
      file_modified_time = format(file_meta$mtime, "%Y-%m-%d %H:%M:%S"),
      manifest_recorded_date = format(Sys.Date(), "%Y-%m-%d"),
      notes = notes,
      stringsAsFactors = FALSE
    )
    
    if (file.exists(manifest_file)) {
      manifest <- read.csv(manifest_file, stringsAsFactors = FALSE, check.names = FALSE)
      manifest <- manifest[manifest$file_name != file_name, , drop = FALSE]
      manifest <- rbind(manifest, new_row)
    } else {
      manifest <- new_row
    }
    
    manifest <- manifest[order(manifest$cohort, manifest$file_name), , drop = FALSE]
    utils::write.csv(manifest, manifest_file, row.names = FALSE)
  }
  
  # Cache
  cache_subdir <- file.path(cache_dir, "gistic_significant", gdac_run_id)
  dir.create(cache_subdir, showWarnings = FALSE, recursive = TRUE)
  
  q_tag <- gsub("\\.", "p", format(gistic_sig_q_cutoff, scientific = FALSE, trim = TRUE))
  source_tag <- gsub("[^A-Za-z0-9]+", "_", gistic_sig_gene_source)
  cache_amp <- file.path(cache_subdir, paste0(cohort, "_q", q_tag, "_", source_tag, "_sig_amp_genes.rds"))
  cache_del <- file.path(cache_subdir, paste0(cohort, "_q", q_tag, "_", source_tag, "_sig_del_genes.rds"))
  
  if (!force_redownload && file.exists(cache_amp) && file.exists(cache_del)) {
    return(list(
      sig_amp_genes = readRDS(cache_amp),
      sig_del_genes = readRDS(cache_del),
      meta = list(
        source = "cache",
        gdac_run_id = gdac_run_id,
        q_cutoff = gistic_sig_q_cutoff,
        gistic_sig_gene_source = gistic_sig_gene_source
      )
    ))
  }
  
  # GDAC cohort directory: usually COHORT-TP, but some cohorts are just COHORT (e.g., SKCM)
  candidate_cohorts <- unique(c(paste0(cohort, cohort_suffix), cohort))
  
  gdac_cohort <- NA_character_
  index_html  <- NA_character_
  
  for (cand in candidate_cohorts) {
    test_url <- paste0("https://gdac.broadinstitute.org/runs/", gdac_run_id, "/data/", cand, "/")
    tmp <- tryCatch(paste(readLines(test_url, warn = FALSE), collapse = "\n"),
                    error = function(e) NA_character_)
    if (!is.na(tmp)) {
      gdac_cohort <- cand
      index_html  <- tmp
      break
    }
  }
  
  if (is.na(index_html)) {
    stop("Failed to read GDAC index for ", cohort,
         ". Tried: ", paste(candidate_cohorts, collapse = ", "))
  }
  
  # Downloaded and extracted resource files
  raw_dir <- file.path(gdac_resources_dir, gdac_run_id, gdac_cohort)
  dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Find the tar.gz name from the GDAC directory index
  base_url <- paste0("https://gdac.broadinstitute.org/runs/", gdac_run_id, "/data/", gdac_cohort, "/")
  
  
  # Resolve date-specific subdirectories when present
  date_dirs <- unique(unlist(regmatches(index_html, gregexpr("[0-9]{8}/", index_html))))
  if (length(date_dirs) > 0) {
    date_dir <- sort(gsub("/", "", date_dirs), decreasing = TRUE)[1]
    base_url <- paste0(base_url, date_dir, "/")
    
    index_html <- tryCatch(paste(readLines(base_url, warn = FALSE), collapse = "\n"),
                           error = function(e) NA_character_)
    if (is.na(index_html)) {
      stop("Failed to read GDAC index inside date dir for ", gdac_cohort, ": ", base_url)
    }
  }
  # -------------------------------------------------------------------------------
  
  # Grab the tarball filename for CopyNumber_Gistic2.Level_4
  m <- regmatches(index_html, gregexpr("gdac\\.broadinstitute\\.org_[^\"<> ]+CopyNumber_Gistic2\\.Level_4[^\"<> ]+\\.tar\\.gz", index_html))
  tar_candidates <- unique(unlist(m))
  
  if (length(tar_candidates) == 0) {
    stop("No CopyNumber_Gistic2.Level_4 tar.gz found at: ", base_url)
  }
  
  # Prefer the most recent-looking filename
  tar_name <- sort(tar_candidates, decreasing = TRUE)[1]
  tar_path <- file.path(raw_dir, tar_name)
  
  # Download tarball if needed
  if (force_redownload || !file.exists(tar_path)) {
    message("Downloading GDAC GISTIC2 Level_4 for ", gdac_cohort, " from ", base_url)
    utils::download.file(url = paste0(base_url, tar_name), destfile = tar_path, mode = "wb", quiet = TRUE)
  }
  
  update_gistic_manifest(
    manifest_file = file.path(gdac_resources_dir, "gistic_download_manifest.csv"),
    cohort = cohort,
    data_type = "GISTIC2 Level 4",
    gdac_run_id = gdac_run_id,
    gistic_module = "CopyNumber_Gistic2.Level_4",
    source_base_url = base_url,
    file_name = tar_name,
    tar_path = tar_path,
    notes = paste(
      "SNP6 Copy number analysis (GISTIC2); package contains",
      "amp_genes.conf_99.txt, del_genes.conf_99.txt,",
      "table_amp.conf_99.txt, table_del.conf_99.txt,",
      "and focal_data_by_genes.txt when present."
    )
  )
  
  # Extract tarball (only if the target files are not already present)
  expected_amp <- list.files(raw_dir, pattern = "amp_genes\\.conf_99\\.txt$", recursive = TRUE, full.names = TRUE)
  expected_del <- list.files(raw_dir, pattern = "del_genes\\.conf_99\\.txt$", recursive = TRUE, full.names = TRUE)
  expected_table_amp <- list.files(raw_dir, pattern = "table_amp\\.conf_99\\.txt$", recursive = TRUE, full.names = TRUE)
  expected_table_del <- list.files(raw_dir, pattern = "table_del\\.conf_99\\.txt$", recursive = TRUE, full.names = TRUE)
  
  expected_missing <- if (gistic_sig_gene_source == "wide_peak") {
    length(expected_amp) == 0 || length(expected_del) == 0
  } else {
    length(expected_table_amp) == 0 || length(expected_table_del) == 0
  }
  
  if (force_redownload || expected_missing) {
    message("Extracting: ", basename(tar_path))
    utils::untar(tar_path, exdir = raw_dir, tar = "internal")
  }
  
  amp_file <- list.files(raw_dir, pattern = "amp_genes\\.conf_99\\.txt$", recursive = TRUE, full.names = TRUE)
  del_file <- list.files(raw_dir, pattern = "del_genes\\.conf_99\\.txt$", recursive = TRUE, full.names = TRUE)
  table_amp_file <- list.files(raw_dir, pattern = "table_amp\\.conf_99\\.txt$", recursive = TRUE, full.names = TRUE)
  table_del_file <- list.files(raw_dir, pattern = "table_del\\.conf_99\\.txt$", recursive = TRUE, full.names = TRUE)
  
  if (gistic_sig_gene_source == "wide_peak" && (length(amp_file) == 0 || length(del_file) == 0)) {
    stop("After extraction, could not find amp_genes.conf_99.txt / del_genes.conf_99.txt under: ", raw_dir)
  }
  if (gistic_sig_gene_source == "genes_in_peak" && (length(table_amp_file) == 0 || length(table_del_file) == 0)) {
    stop("After extraction, could not find table_amp.conf_99.txt / table_del.conf_99.txt under: ", raw_dir)
  }
  
  amp_file <- if (length(amp_file) > 0) amp_file[1] else NA_character_
  del_file <- if (length(del_file) > 0) del_file[1] else NA_character_
  table_amp_file <- if (length(table_amp_file) > 0) table_amp_file[1] else NA_character_
  table_del_file <- if (length(table_del_file) > 0) table_del_file[1] else NA_character_
  
  parse_wide_peak_genes <- function(path, q_cutoff) {
    gistic_file <- read.table(
      path,
      sep = "\t",
      header = TRUE,
      fill = TRUE,
      quote = "",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    
    q_vals <- suppressWarnings(as.numeric(gistic_file[gistic_file[, 1] == "q value", -1]))
    if (length(q_vals) == 0) return(character(0))
    
    sig_peaks <- which(q_vals <= q_cutoff) + 1
    if (length(sig_peaks) == 0) return(character(0))
    
    start_row <- which(gistic_file[, 1] == "genes in wide peak")
    if (length(start_row) == 0) return(character(0))
    
    gene_block <- gistic_file[start_row:nrow(gistic_file), sig_peaks, drop = FALSE]
    unique(na.omit(unlist(gene_block)))
  }
  
  parse_genes_in_peak <- function(path, q_cutoff) {
    peak_table <- read.table(
      path,
      sep = "\t",
      header = TRUE,
      fill = TRUE,
      quote = "",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    
    if (!"genes_in_peak" %in% colnames(peak_table)) {
      return(character(0))
    }
    
    q_col <- intersect(c("q value", "q_value", "q.values", "q.values.residual"), colnames(peak_table))
    if (length(q_col) > 0) {
      q_vals <- suppressWarnings(as.numeric(peak_table[[q_col[1]]]))
      peak_table <- peak_table[!is.na(q_vals) & q_vals <= q_cutoff, , drop = FALSE]
    }
    
    genes <- unlist(strsplit(peak_table$genes_in_peak, ",", fixed = TRUE))
    genes <- trimws(genes)
    unique(genes[!is.na(genes) & genes != ""])
  }
  
  if (gistic_sig_gene_source == "wide_peak") {
    sig_amp_genes <- parse_wide_peak_genes(amp_file, gistic_sig_q_cutoff)
    sig_del_genes <- parse_wide_peak_genes(del_file, gistic_sig_q_cutoff)
  } else {
    sig_amp_genes <- parse_genes_in_peak(table_amp_file, gistic_sig_q_cutoff)
    sig_del_genes <- parse_genes_in_peak(table_del_file, gistic_sig_q_cutoff)
  }
  
  saveRDS(sig_amp_genes, cache_amp)
  saveRDS(sig_del_genes, cache_del)
  
  list(
    sig_amp_genes = sig_amp_genes,
    sig_del_genes = sig_del_genes,
    meta = list(
      source = "gdac",
      gdac_run_id = gdac_run_id,
      gdac_cohort = gdac_cohort,
      q_cutoff = gistic_sig_q_cutoff,
      gistic_sig_gene_source = gistic_sig_gene_source,
      amp_file = amp_file,
      del_file = del_file,
      table_amp_file = table_amp_file,
      table_del_file = table_del_file,
      resolved_base_url = base_url
    )
  )
}
