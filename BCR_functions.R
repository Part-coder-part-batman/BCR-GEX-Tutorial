# =============================================================================
# BCR_functions.R
# Helper functions for BCR QC, clonal assignment, and mutation analysis
#
# These functions are sourced by BCR_GEX_Tutorial_Part1.Rmd -- do not run this file
# directly. All required packages are loaded by the pipeline script.
#
# Functions in this file (in order of use):
#   1. process_bcr_data()                    -- load AIRR TSV + Cell Ranger annotations
#   2. filter_bcr_data()                     -- QC filtering: productive, paired, deduplicated
#   3. plot_clonal_thresholds()              -- compute and visualize clonal threshold
#   4. define_clones_basic()                 -- hierarchical clonal assignment
#   5. annotate_clone_ids()                  -- human-readable clone ID encoding
#   6. qc_clone_assignment()                 -- verify clone assignment output
#   7. visualize_clones()                    -- rank-abundance, size, diversity plots
#   8. reconstruct_germlines_and_mutations() -- germline reconstruction + SHM
#   9. plot_mutation_frequencies()           -- mutation frequency plots
# =============================================================================


# =============================================================================
# 1. process_bcr_data
# =============================================================================

#' Load and prepare BCR data from Immcantation and Cell Ranger output
#'
#' Reads the AIRR-format TSV produced by MakeDb.py and the filtered contig
#' annotations CSV from Cell Ranger. Merges UMI counts into the BCR table,
#' then appends the sample ID to all cell and sequence IDs to prevent barcode
#' collisions when combining multiple samples downstream.
#'
#' @param airr_file Path to the AIRR TSV file ending in _airr_db-pass.tsv
#' @param annotations_file Path to filtered_contig_annotations.csv from Cell Ranger
#' @param sample_id Character string identifying this sample (e.g. "P1_LN")
#' @return A data frame of BCR contigs with updated cell_id and sequence_id
#'
#' @examples
#' bcr_data <- process_bcr_data(
#'   airr_file        = "P1_LN/results/P1_LN_airr_db-pass.tsv",
#'   annotations_file = "P1_LN/filtered_contig_annotations.csv",
#'   sample_id        = "P1_LN"
#' )

process_bcr_data <- function(airr_file, annotations_file, sample_id) {

  message("Loading AIRR BCR data from: ", airr_file)
  bcr_data <- airr::read_rearrangement(
    airr_file,
    aux_types = c(
      "v_germline_length" = "i",
      "d_germline_length" = "i",
      "j_germline_length" = "i"
    )
  )
  message("BCR data loaded. Total contigs: ", nrow(bcr_data))

  message("Reading Cell Ranger annotations from: ", annotations_file)
  annotations <- readr::read_csv(annotations_file, show_col_types = FALSE)

  if (!"contig_id" %in% colnames(annotations) || !"umis" %in% colnames(annotations)) {
    stop("Annotation file must contain 'contig_id' and 'umis' columns.")
  }
  message("Annotations loaded. Total contigs: ", nrow(annotations))

  # Merge UMI counts into BCR table
  umis <- annotations %>%
    dplyr::select(contig_id, umi_count = umis)

  bcr_data <- bcr_data %>%
    dplyr::left_join(umis, by = c("sequence_id" = "contig_id"))

  # Append sample ID to cell and sequence IDs
  # This prevents barcode collisions when multiple samples are combined later.
  # Example: "ACGT...TGCA-1_contig_1" becomes "ACGT...TGCA-1_P1_LN_contig_1"
  bcr_data$sample_id  <- sample_id
  base_cell_id        <- gsub("_contig_\\d+.*", "", bcr_data$sequence_id)
  bcr_data$cell_id    <- paste0(base_cell_id, "_", sample_id)
  bcr_data$sequence_id <- gsub("(_contig_\\d+)", paste0("_", sample_id, "\\1"),
                                bcr_data$sequence_id)

  message("Sample ID appended to cell and sequence IDs.")
  message("Example cell_id: ", bcr_data$cell_id[1])
  message("process_bcr_data complete.")

  return(bcr_data)
}


# =============================================================================
# 2. filter_bcr_data
# =============================================================================

#' Filter BCR contigs to retain productive, well-paired sequences
#'
#' Applies four sequential filters:
#'   1. Remove non-productive sequences
#'   2. Remove cells with more than one heavy chain (likely doublets)
#'   3. Remove light chains without a paired heavy chain
#'   4. For cells with multiple light chains, retain only the one with the
#'      highest UMI count
#'
#' Prints a QC summary of chain pairing at the end.
#'
#' @param bcr_data Data frame from process_bcr_data()
#' @return Filtered BCR data frame ready for clonal assignment
#'
#' @examples
#' bcr_data <- filter_bcr_data(bcr_data)

filter_bcr_data <- function(bcr_data) {

  # ---- Step 1: Remove non-productive sequences ------------------------------
  message("Step 1: Removing non-productive sequences...")
  n_before <- nrow(bcr_data)
  bcr_data <- dplyr::filter(bcr_data, productive)
  message("  Removed: ", n_before - nrow(bcr_data), " | Remaining: ", nrow(bcr_data))

  # ---- Step 2: Remove cells with multiple heavy chains ----------------------
  message("Step 2: Removing cells with multiple heavy chains...")
  multi_heavy_cells <- names(which(
    table(dplyr::filter(bcr_data, locus == "IGH")$cell_id) > 1
  ))
  message("  Cells with >1 heavy chain: ", length(multi_heavy_cells))
  n_before <- nrow(bcr_data)
  bcr_data <- dplyr::filter(bcr_data, !cell_id %in% multi_heavy_cells)
  message("  Removed: ", n_before - nrow(bcr_data), " | Remaining: ", nrow(bcr_data))

  # ---- Step 3: Remove light chains without a paired heavy chain -------------
  message("Step 3: Removing unpaired light chains...")
  bcr_data <- bcr_data %>%
    dplyr::mutate(base_cell_id = gsub("_contig_\\d+.*", "", cell_id))

  cell_summary <- bcr_data %>%
    dplyr::group_by(base_cell_id) %>%
    dplyr::summarize(
      has_heavy = any(locus == "IGH"),
      has_light = any(locus %in% c("IGK", "IGL")),
      .groups = "drop"
    )

  unpaired_light_cells <- cell_summary %>%
    dplyr::filter(!has_heavy & has_light) %>%
    dplyr::pull(base_cell_id)

  message("  Cells with light chain but no heavy chain: ", length(unpaired_light_cells))
  bcr_data <- dplyr::filter(bcr_data, !base_cell_id %in% unpaired_light_cells)
  bcr_data <- dplyr::select(bcr_data, -base_cell_id)
  message("  Remaining: ", nrow(bcr_data))

  # ---- Step 4: Retain highest-UMI light chain per cell ----------------------
  message("Step 4: Collapsing multiple light chains -- keeping highest UMI per cell...")
  light <- bcr_data %>%
    dplyr::filter(locus %in% c("IGK", "IGL")) %>%
    dplyr::group_by(cell_id) %>%
    dplyr::slice_max(umi_count, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

  heavy <- dplyr::filter(bcr_data, locus == "IGH")
  bcr_data <- dplyr::bind_rows(heavy, light)
  message("  Remaining after light chain collapse: ", nrow(bcr_data))

  # QC check: no cell should have more than one light chain
  light_chain_qc <- bcr_data %>%
    dplyr::filter(locus %in% c("IGK", "IGL")) %>%
    dplyr::count(cell_id) %>%
    dplyr::filter(n > 1)

  if (nrow(light_chain_qc) == 0) {
    message("QC passed: all cells have at most one light chain.")
  } else {
    warning("QC failed: ", nrow(light_chain_qc), " cells still have >1 light chain.")
    print(light_chain_qc)
  }

  # ---- Chain pairing summary ------------------------------------------------
  heavy_cells  <- unique(dplyr::filter(bcr_data, locus == "IGH")$cell_id)
  light_cells  <- unique(dplyr::filter(bcr_data, locus %in% c("IGK", "IGL"))$cell_id)
  paired_cells <- intersect(heavy_cells, light_cells)

  pairing_summary <- data.frame(
    Status = c("Heavy + Light (paired)", "Heavy only", "Light only"),
    Cells  = c(length(paired_cells),
               length(setdiff(heavy_cells, light_cells)),
               length(setdiff(light_cells, heavy_cells)))
  )
  message("Chain pairing summary:")
  print(pairing_summary)

  return(bcr_data)
}


# =============================================================================
# 3. plot_clonal_thresholds
# =============================================================================

#' Compute and visualize the clonal distance threshold
#'
#' Calculates pairwise Hamming distances between IGH sequences and fits a
#' Gaussian-mixture model (via SHazaM) to identify the threshold separating
#' clonally related from unrelated sequences.
#'
#' Returns two plots:
#'   - A histogram of nearest-neighbor distances with the automatic threshold
#'   - The SHazaM model fit used to determine that threshold
#'
#' @param bcr_data Filtered BCR data frame (from filter_bcr_data())
#' @param user_defined_threshold Optional numeric value to overlay as a
#'   reference line. Set to NULL to omit. Default: NULL
#' @param binwidth Histogram bin width. Default: 0.02
#' @return A named list with elements: user_plot, shazam_plot, auto_threshold
#'
#' @examples
#' thr_res   <- plot_clonal_thresholds(bcr_data, user_defined_threshold = NULL)
#' threshold <- thr_res$auto_threshold
#' print(thr_res$shazam_plot)

plot_clonal_thresholds <- function(bcr_data, user_defined_threshold = NULL, binwidth = 0.02) {

  message("Calculating nearest-neighbor distances for IGH sequences...")
  dist_nearest <- shazam::distToNearest(dplyr::filter(bcr_data, locus == "IGH"))

  if (all(is.na(dist_nearest$dist_nearest))) {
    stop("All dist_nearest values are NA -- check that IGH sequences are present.")
  }

  # ---- Plot 1: Histogram with optional user threshold -----------------------
  p_user <- ggplot2::ggplot(
    subset(dist_nearest, !is.na(dist_nearest)),
    ggplot2::aes(x = dist_nearest)
  ) +
    ggplot2::geom_histogram(fill = "lightblue", color = "white", binwidth = binwidth) +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.1)) +
    ggplot2::labs(
      x     = "Hamming distance to nearest neighbor",
      y     = "Count",
      title = "Nearest-Neighbor Distance Distribution"
    ) +
    ggplot2::theme_bw()

  if (!is.null(user_defined_threshold)) {
    p_user <- p_user +
      ggplot2::geom_vline(
        xintercept = user_defined_threshold,
        color = "red", linetype = "dashed", linewidth = 1
      )
  }

  # ---- Step 2: Automatic threshold via GMM ----------------------------------
  message("Fitting GMM to determine automatic threshold...")
  threshold_output <- shazam::findThreshold(
    dist_nearest$dist_nearest,
    method = "gmm",
    model  = "gamma-norm",
    cutoff = "user",
    spc    = 0.995
  )
  auto_threshold <- threshold_output@threshold
  message("Automatic threshold: ", signif(auto_threshold, 4))

  # ---- Plot 2: SHazaM model fit with threshold line -------------------------
  p_shazam <- plot(threshold_output, binwidth = binwidth, silent = TRUE) +
    ggplot2::geom_vline(
      xintercept = auto_threshold,
      color = "darkblue", linetype = "dashed", linewidth = 1
    ) +
    ggplot2::labs(title = "SHazaM GMM Threshold Fit") +
    ggplot2::theme_bw()

  return(list(
    user_plot      = p_user,
    shazam_plot    = p_shazam,
    auto_threshold = auto_threshold
  ))
}


# =============================================================================
# 4. define_clones_basic
# =============================================================================

#' Assign BCR sequences to clones using hierarchical clustering
#'
#' Wraps scoper::hierarchicalClones() with settings appropriate for
#' multi-sample 10x data: clustering is driven by heavy chains only
#' (only_heavy = TRUE), with light chain splitting enabled to resolve
#' ambiguous heavy-chain clones that carry different light chains
#' (split_light = TRUE).
#'
#' @param bcr_data Filtered BCR data frame (from filter_bcr_data())
#' @param threshold Numeric distance threshold (from plot_clonal_thresholds())
#' @param sample_id Character string identifying this sample
#' @return BCR data frame with clone_id column added
#'
#' @examples
#' bcr_data <- define_clones_basic(bcr_data, threshold = threshold, sample_id = "P1_LN")

define_clones_basic <- function(bcr_data, threshold, sample_id) {

  message("Running hierarchical clonal assignment for sample: ", sample_id)

  bcr_data$subject_id <- sample_id
  bcr_data <- bcr_data %>%
    dplyr::mutate(cell_id = gsub("_contig_\\d+.*", "", cell_id))

  results <- scoper::hierarchicalClones(
    bcr_data,
    cell_id         = "cell_id",
    threshold       = threshold,
    only_heavy      = TRUE,
    split_light     = TRUE,
    summarize_clones = FALSE,
    fields          = "subject_id"
  )

  n_clones   <- length(unique(results$clone_id))
  n_excluded <- length(setdiff(bcr_data$cell_id, results$cell_id))
  message("Clones assigned: ", n_clones)
  message("Cells excluded from clonal assignment: ", n_excluded)

  return(results)
}


# =============================================================================
# 5. annotate_clone_ids
# =============================================================================

#' Replace numeric clone IDs with human-readable encoded identifiers
#'
#' The default clone_id from hierarchicalClones() is an integer. This function
#' replaces it with an informative string encoding four pieces of information:
#'
#'   SAMPLEID_RANDOMCODE_CLONESIZE_ISOTYPES
#'
#' Example: "P1_LN_aB3x_14_G1G2M"
#'
#' The random 4-character code ensures uniqueness across samples when data
#' are later merged, while the encoded clone size and isotypes provide
#' at-a-glance information without requiring a separate lookup table.
#'
#' @param results BCR data frame with numeric clone_id (from define_clones_basic())
#' @param sample_id Character string identifying this sample
#' @return BCR data frame with informative clone_id and clone_count columns
#'
#' @examples
#' bcr_data <- annotate_clone_ids(bcr_data, sample_id = "P1_LN")

annotate_clone_ids <- function(results, sample_id) {

  message("Annotating clone IDs for sample: ", sample_id)

  # Simplify c_call to gene-level (strip allele, collapse lambda subtypes)
  simplify_c_call <- function(c_call) {
    sapply(strsplit(c_call, ","), function(x) {
      gene <- strsplit(x[1], "\\*")[[1]][1]
      if (grepl("^IGLC", gene)) return("IGLC")
      return(gene)
    })
  }
  results$c_call <- simplify_c_call(results$c_call)

  # Generate a unique 4-character random code per clone
  unique_clones  <- unique(results$clone_id)
  # Exclude E and e to prevent codes like "1E82" being read as scientific notation
  safe_chars <- c(LETTERS[LETTERS != "E"], letters[letters != "e"], as.character(0:9))
  random_codes <- setNames(
    replicate(length(unique_clones), {
      paste(sample(safe_chars, 4, replace = TRUE), collapse = "")
    }),
    unique_clones
  )

  # Compute per-clone stats from IGH rows only
  heavy_data <- results[results$locus == "IGH", ]
  clone_stats <- lapply(unique(heavy_data$clone_id), function(cid) {
    cd          <- heavy_data[heavy_data$clone_id == cid, ]
    clone_count <- length(unique(cd$cell_id))
    isotypes    <- paste(
      sort(gsub("^IGH", "", unique(cd$c_call[!is.na(cd$c_call)]))),
      collapse = ""
    )
    data.frame(clone_id = cid, clone_count = clone_count,
               isotypes = isotypes, stringsAsFactors = FALSE)
  })
  clone_stats <- do.call(rbind, clone_stats)

  # Handle any clones without heavy chain rows
  missing_heavy <- setdiff(unique_clones, clone_stats$clone_id)
  if (length(missing_heavy) > 0) {
    clone_stats <- rbind(clone_stats, data.frame(
      clone_id = missing_heavy, clone_count = 0L,
      isotypes = "NoHeavy", stringsAsFactors = FALSE
    ))
  }

  # Build new informative clone ID
  clone_stats$random_code  <- random_codes[clone_stats$clone_id]
  clone_stats$new_clone_id <- paste(
    sample_id, clone_stats$random_code,
    clone_stats$clone_count, clone_stats$isotypes,
    sep = "_"
  )

  # Merge back and replace clone_id
  results <- merge(
    results,
    clone_stats[, c("clone_id", "new_clone_id", "clone_count")],
    by = "clone_id", all.x = TRUE
  )
  results$clone_id     <- results$new_clone_id
  results$new_clone_id <- NULL

  message("Clone ID annotation complete. Total clones: ",
          length(unique(results$clone_id)))

  # Show a few example clone IDs so the user can verify the encoding
  # Format: SAMPLEID_RANDOMCODE_CLONESIZE_ISOTYPES
  example_clones <- results %>%
    dplyr::filter(locus == "IGH") %>%
    dplyr::group_by(clone_id) %>%
    dplyr::summarise(
      size    = dplyr::n_distinct(cell_id),
      isotype = paste(sort(unique(c_call[!is.na(c_call)])), collapse = "/"),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(size)) %>%
    head(5)

  message("Example clone IDs (format: SAMPLEID_CODE_SIZE_ISOTYPES):")
  print(example_clones)

  return(results)
}


# =============================================================================
# 6. qc_clone_assignment
# =============================================================================

#' QC check on clone assignment output
#'
#' Verifies that:
#'   - The number of cells assigned to clones is as expected
#'   - All clone IDs follow the SAMPLEID_RANDOM_SIZE_ISOTYPES format
#'   - The clone size encoded in each ID matches the actual cell count
#'
#' @param bcr_data The filtered BCR data before clonal assignment
#' @param results The BCR data after clonal assignment and annotation
#' @return Invisibly returns a summary table of clone size consistency
#'
#' @examples
#' qc_clone_assignment(bcr_data_pre, bcr_data_post)

qc_clone_assignment <- function(bcr_data, results) {

  message("Running clone assignment QC...")

  # Cell counts
  n_start <- length(unique(bcr_data$cell_id))
  n_final <- length(unique(results$cell_id))
  message("  Cells before clonal assignment: ", n_start)
  message("  Cells assigned to clones:       ", n_final)
  message("  Cells lost:                     ", n_start - n_final)

  # Clone ID format check
  # Format is SAMPLEID_RANDOMCODE_SIZE_ISOTYPES. The sample ID can itself
  # contain underscores (e.g. "P1_LN"), so counting underscores is unreliable.
  # Match the invariant tail instead: a 4-character code, then the numeric
  # size, then the isotype field.
  format_ok <- all(grepl("_[A-Za-z0-9]{4}_[0-9]+_[^_]+$", results$clone_id))
  if (format_ok) {
    message("  Clone ID format: OK (SAMPLEID_RANDOM_SIZE_ISOTYPES)")
  } else {
    warning("  Some clone IDs do not follow the expected format.")
  }

  # Encoded size vs actual size
  # Extract the numeric size field directly. It is the digit run immediately
  # preceded by an underscore and immediately followed by the final isotype
  # field, so this is robust to underscores inside the sample ID.
  results <- results %>%
    dplyr::mutate(
      clone_size_encoded = as.numeric(
        stringr::str_extract(clone_id, "(?<=_)[0-9]+(?=_[^_]+$)")
      )
    )

  size_check <- results %>%
    dplyr::group_by(clone_id) %>%
    dplyr::summarise(
      actual_size  = dplyr::n_distinct(cell_id),
      encoded_size = dplyr::first(clone_size_encoded),
      match        = actual_size == encoded_size,
      .groups      = "drop"
    )

  # Count an unparseable size (NA) as a mismatch rather than dropping it,
  # so a malformed clone ID surfaces here instead of passing silently.
  mismatches <- dplyr::filter(size_check, !match | is.na(match))
  if (nrow(mismatches) == 0) {
    message("  Encoded clone sizes match actual sizes: OK")
  } else {
    warning("  ", nrow(mismatches), " clone(s) have mismatched encoded vs actual size.")
  }

  message("  Total unique clones: ", dplyr::n_distinct(results$clone_id))

  invisible(size_check)
}


# =============================================================================
# 7. visualize_clones
# =============================================================================

#' Plot clone rank-abundance, size distribution, and clonal diversity
#'
#' All three plots use heavy chain sequences only, grouped by sample_id.
#'
#' @param bcr_data Clone-assigned BCR data (must have clone_id and sample_id)
#' @param nboot Number of bootstrap replicates for abundance/diversity. Default: 100
#' @return Named list of ggplot objects: rank_abundance_plot, clone_size_plot,
#'   diversity_plot
#'
#' @examples
#' clone_plots <- visualize_clones(bcr_data)
#' print(clone_plots$rank_abundance_plot)

visualize_clones <- function(bcr_data, nboot = 100) {

  if (!all(c("clone_id", "sample_id") %in% colnames(bcr_data))) {
    stop("Input must contain 'clone_id' and 'sample_id' columns.")
  }

  heavy_data <- dplyr::filter(bcr_data, locus == "IGH")
  if (nrow(heavy_data) == 0) stop("No IGH sequences found in input.")

  # Rank-abundance
  message("Calculating rank-abundance curves...")
  abund <- alakazam::estimateAbundance(heavy_data, group = "sample_id", nboot = nboot)
  rank_abundance_plot <- plot(abund, silent = TRUE) +
    ggplot2::facet_wrap(~sample_id) +
    ggplot2::labs(title = "Clone Rank-Abundance", x = "Clone Rank", y = "Abundance") +
    ggplot2::theme_bw()

  # Clone size distribution
  message("Calculating clone size distribution...")
  clone_sizes <- alakazam::countClones(heavy_data, groups = "sample_id")
  clone_size_plot <- ggplot2::ggplot(clone_sizes, ggplot2::aes(x = seq_count)) +
    ggplot2::geom_bar(fill = "skyblue", color = "black") +
    ggplot2::facet_wrap(~sample_id) +
    ggplot2::labs(title = "Clone Size Distribution",
                  x = "Sequences per clone", y = "Count") +
    ggplot2::theme_bw()

  # Clonal diversity
  message("Calculating clonal diversity...")
  diversity <- alakazam::alphaDiversity(heavy_data, group = "sample_id", nboot = nboot)
  diversity_plot <- plot(diversity, silent = TRUE) +
    ggplot2::facet_wrap(~sample_id) +
    ggplot2::labs(title = "Clonal Diversity") +
    ggplot2::theme_bw()

  message("visualize_clones complete.")
  return(list(
    rank_abundance_plot = rank_abundance_plot,
    clone_size_plot     = clone_size_plot,
    diversity_plot      = diversity_plot
  ))
}


# =============================================================================
# 8. reconstruct_germlines_and_mutations
# =============================================================================

#' Reconstruct germline sequences and calculate somatic hypermutation (SHM)
#'
#' Uses dowser::createGermlines() to infer the unmutated germline for each
#' sequence, then calculates both mutation frequency and raw mutation count
#' across the V gene region (IMGT_V definition) using shazam::observedMutations().
#'
#' Light chain clone IDs are resolved by inheriting the heavy chain clone ID
#' for the same cell, which is required for correct germline reconstruction.
#'
#' @param bcr_data Clone-assigned BCR data (from annotate_clone_ids())
#' @param reference_dir Path to IMGT germline VDJ reference directory
#' @param nproc Number of CPU cores. Default: 1
#' @return Named list:
#'   - results_with_mut: BCR data frame with mu_freq and mu_count columns added
#'   - mut_freq_by_clone: per-clone median mutation frequency summary
#'   - mut_histogram_plot: ggplot histogram of median mutation frequencies
#'
#' @examples
#' mut_outputs <- reconstruct_germlines_and_mutations(
#'   bcr_data,
#'   reference_dir = "C:/Users/YourName/Documents/Immcantation/imgt/human/vdj",
#'   nproc         = 1
#' )
#' bcr_data <- mut_outputs$results_with_mut

reconstruct_germlines_and_mutations <- function(bcr_data, reference_dir, nproc = 1) {

  message("Loading IMGT germline references from: ", reference_dir)
  references <- dowser::readIMGT(dir = reference_dir)

  # Resolve light chain clone IDs from paired heavy chains
  bcr_data <- bcr_data %>%
    dplyr::mutate(base_cell_id = gsub("_contig_\\d+.*", "", sequence_id))

  heavy <- dplyr::filter(bcr_data, locus == "IGH")
  light <- dplyr::filter(bcr_data, locus %in% c("IGK", "IGL"))

  heavy_clone_map <- heavy %>%
    dplyr::group_by(base_cell_id) %>%
    dplyr::summarise(clone_id_heavy = dplyr::first(clone_id), .groups = "drop")

  light <- light %>%
    dplyr::left_join(heavy_clone_map, by = "base_cell_id") %>%
    dplyr::mutate(clone_id = dplyr::coalesce(clone_id, clone_id_heavy)) %>%
    dplyr::select(-clone_id_heavy)

  # Remove light chains that mapped to multiple cells (ambiguous)
  ambiguous <- light %>%
    dplyr::group_by(sequence_id) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::pull(sequence_id)

  light <- dplyr::filter(light, !sequence_id %in% ambiguous)

  bcr_clean <- dplyr::bind_rows(heavy, light) %>%
    dplyr::select(-base_cell_id)

  # Simplify c_call
  simplify_c_call <- function(c_call) {
    sapply(strsplit(c_call, ","), function(x) {
      gene <- strsplit(x[1], "\\*")[[1]][1]
      if (grepl("^IGLC", gene)) return("IGLC")
      return(gene)
    })
  }
  bcr_clean$c_call <- simplify_c_call(bcr_clean$c_call)

  # Reconstruct germlines
  message("Reconstructing germline sequences...")
  bcr_germlines <- dowser::createGermlines(
    bcr_clean,
    references = references,
    fields     = "subject_id",
    nproc      = nproc
  )

  # Calculate mutation frequency
  message("Calculating mutation frequencies...")
  data_mut_freq <- shazam::observedMutations(
    bcr_germlines,
    sequenceColumn  = "sequence_alignment",
    germlineColumn  = "germline_alignment_d_mask",
    regionDefinition = IMGT_V,
    frequency       = TRUE,
    combine         = TRUE,
    nproc           = nproc
  )

  # Calculate raw mutation count
  message("Calculating mutation counts...")
  data_mut_count <- shazam::observedMutations(
    bcr_germlines,
    sequenceColumn  = "sequence_alignment",
    germlineColumn  = "germline_alignment_d_mask",
    regionDefinition = IMGT_V,
    frequency       = FALSE,
    combine         = TRUE,
    nproc           = nproc
  )

  # Merge mutation metrics back
  mut_metrics <- dplyr::left_join(
    dplyr::select(data_mut_freq, sequence_id, mu_freq),
    dplyr::select(data_mut_count, sequence_id, mu_count),
    by = "sequence_id"
  )
  results_with_mut <- dplyr::left_join(bcr_germlines, mut_metrics, by = "sequence_id")

  # Per-clone median mutation frequency
  mut_freq_by_clone <- results_with_mut %>%
    dplyr::group_by(clone_id, locus) %>%
    dplyr::summarise(
      median_mut_freq = median(mu_freq, na.rm = TRUE),
      .groups = "drop"
    )

  mut_histogram_plot <- ggplot2::ggplot(
    mut_freq_by_clone,
    ggplot2::aes(x = median_mut_freq, fill = locus)
  ) +
    ggplot2::geom_histogram(binwidth = 0.005, color = "black", alpha = 0.7,
                            position = "dodge") +
    ggplot2::labs(
      x     = "Median mutation frequency",
      y     = "Count",
      fill  = "Locus",
      title = "Mutation Frequency by Clone"
    ) +
    ggplot2::theme_bw()

  message("reconstruct_germlines_and_mutations complete.")
  return(list(
    results_with_mut   = results_with_mut,
    mut_freq_by_clone  = mut_freq_by_clone,
    mut_histogram_plot = mut_histogram_plot
  ))
}


# =============================================================================
# 9. plot_mutation_frequencies
# =============================================================================

#' Plot somatic hypermutation frequency distributions
#'
#' Produces three plots:
#'   1. Histogram of per-clone median mutation frequency, split by locus
#'   2. Boxplot of per-sequence mutation frequency by subject
#'   3. Boxplot of per-sequence mutation frequency by isotype
#'
#' @param bcr_data BCR data frame with mu_freq, mu_count, clone_id, locus,
#'   subject_id, and c_call columns (from reconstruct_germlines_and_mutations())
#' @param binwidth Bin width for histogram. Default: 0.005
#' @return Named list of ggplot objects: histogram, by_subject, by_isotype
#'
#' @examples
#' mf_plots <- plot_mutation_frequencies(bcr_data)
#' print(mf_plots$by_isotype)

plot_mutation_frequencies <- function(bcr_data, binwidth = 0.005) {

  # Per-clone median histogram
  mut_freq_clone <- bcr_data %>%
    dplyr::group_by(clone_id, locus) %>%
    dplyr::summarize(
      median_mut_freq = median(mu_freq, na.rm = TRUE),
      .groups = "drop"
    )

  histogram <- ggplot2::ggplot(
    mut_freq_clone,
    ggplot2::aes(x = median_mut_freq, fill = locus)
  ) +
    ggplot2::geom_histogram(binwidth = binwidth, color = "black",
                            alpha = 0.7, position = "dodge") +
    ggplot2::labs(x = "Median mutation frequency", y = "Count", fill = "Locus",
                  title = "Mutation Frequency Distribution by Clone") +
    ggplot2::theme_bw()

  # By subject
  by_subject <- ggplot2::ggplot(
    bcr_data,
    ggplot2::aes(y = mu_freq, x = subject_id, fill = locus)
  ) +
    ggplot2::geom_boxplot(outlier.size = 0.5) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.3, color = "darkgray") +
    ggplot2::labs(x = "Subject", y = "Mutation frequency", fill = "Locus",
                  title = "Mutation Frequency by Subject") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  # By isotype
  by_isotype <- ggplot2::ggplot(
    bcr_data,
    ggplot2::aes(y = mu_freq, x = c_call, fill = locus)
  ) +
    ggplot2::geom_boxplot(outlier.size = 0.5) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.3, color = "darkgray") +
    ggplot2::labs(x = "Isotype", y = "Mutation frequency", fill = "Locus",
                  title = "Mutation Frequency by Isotype") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  message("plot_mutation_frequencies complete.")
  return(list(
    histogram  = histogram,
    by_subject = by_subject,
    by_isotype = by_isotype
  ))
}
