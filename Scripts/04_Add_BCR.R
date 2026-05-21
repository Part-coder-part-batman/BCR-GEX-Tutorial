# =============================================================================
# 04_Add_BCR.R
# Add BCR metadata to GEX Seurat objects
# =============================================================================
#
# HOW TO USE:
#   1. Fill in the sample definitions and output directory below
#   2. Run the script (Ctrl+Shift+Enter in RStudio to source the whole file)
#
# INPUT FILES (one per sample):
#   - <SAMPLE_ID>_S3_postDoublet.rds   produced by 03_Doublet_Removal.R
#       Provide the full path as rds_in in the sample definitions below.
#   - <SAMPLE_ID>_bcr_data.tsv         produced by 01_BCR_Pipeline.R
#       Provide the full path as bcr_tsv in the sample definitions below.
#
# OUTPUT:
#   - <SAMPLE_ID>_S4_postBCR.rds
#       Saved to the full path you provide as rds_out in the sample definitions.
#
# WHAT THIS SCRIPT DOES:
#   Loads the BCR TSV produced in step 01 and joins it onto the Seurat object
#   by matching cell barcodes. Each cell receives BCR metadata columns:
#   isotype (c_call), clone ID, clone size, V/J genes, and mutation frequency
#   for both heavy and light chains. Cells with no BCR match are kept but
#   flagged as BCR = FALSE.
#
#   After joining, produces two summary plots per run:
#   - BCR+ vs BCR- fraction per sample (percentage stacked bar)
#   - Isotype distribution per sample (percentage bar chart, BCR+ cells only)
#
# NOTE -- BARCODE MATCHING:
#   Cell barcodes in the Seurat object were suffixed with the sample name in
#   step 02 (e.g. ACGT...TGCA-1_MySample1). The BCR TSV cell_id column must
#   use the same suffix. If you used 01_BCR_Pipeline.R this was done
#   automatically. If match rates are unexpectedly low (<30%), check that
#   the suffixes are consistent between the two pipelines.
#
# NOTE -- BCR- CELLS:
#   Cells with BCR = FALSE are genuine B cells (CD79A+) that did not yield a
#   detectable BCR sequence. They are kept in the object and will contribute
#   to GEX analysis. All BCR metadata columns are NA for these cells.
#
# SEURAT VERSION NOTE:
#   Compatible with Seurat v4 and v5.
# =============================================================================


# ---- 0. PACKAGES ------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(data.table)
  library(tibble)
  library(ggplot2)
  library(scales)
})


# =============================================================================
# !! ONLY MODIFY THIS SECTION !!
# =============================================================================

# --- Sample definitions -------------------------------------------------------
# Add one entry per sample.
# rds_in  : full path to the _S3_postDoublet.rds file produced by 03_Doublet_Removal.R
# rds_out : full path where the _S4_postBCR.rds file will be saved
# bcr_tsv : full path to the _bcr_data.tsv file produced by 01_BCR_Pipeline.R

samples <- list(
  list(
    name    = "MySample1",
    rds_in  = "C:/Users/YourName/Documents/MyProject/RDS_Objects/MySample1_S3_postDoublet.rds",
    rds_out = "C:/Users/YourName/Documents/MyProject/RDS_Objects/MySample1_S4_postBCR.rds",
    bcr_tsv = "C:/Users/YourName/Documents/MyProject/MySample1/Output/MySample1_bcr_data.tsv"
  ),
  list(
    name    = "MySample2",
    rds_in  = "C:/Users/YourName/Documents/MyProject/RDS_Objects/MySample2_S3_postDoublet.rds",
    rds_out = "C:/Users/YourName/Documents/MyProject/RDS_Objects/MySample2_S4_postBCR.rds",
    bcr_tsv = "C:/Users/YourName/Documents/MyProject/MySample2/Output/MySample2_bcr_data.tsv"
  )
)

# =============================================================================
# END OF MODIFIABLE SECTION -- do not change anything below
# =============================================================================


if (!dir.exists(dirname(samples[[1]]$rds_out))) {
  dir.create(dirname(samples[[1]]$rds_out), recursive = TRUE)
}

summary_rows <- list()


# =============================================================================
# LOOP OVER SAMPLES
# =============================================================================

for (s in samples) {

  message("\n", paste(rep("=", 60), collapse = ""))
  message("Processing sample: ", s$name)
  message(paste(rep("=", 60), collapse = ""))

  # ---- 1. LOAD SEURAT OBJECT -------------------------------------------------
  if (!file.exists(s$rds_in)) {
    stop("Cannot find RDS file:\n  ", s$rds_in,
         "\nCheck rds_in for sample '", s$name, "' and confirm 03_Doublet_Removal.R has been run.")
  }
  obj <- readRDS(s$rds_in)
  message("  B cells loaded (post-doublet removal): ", ncol(obj))


  # ---- 2. LOAD BCR TSV -------------------------------------------------------
  if (!file.exists(s$bcr_tsv)) {
    stop("Cannot find BCR TSV:\n  ", s$bcr_tsv,
         "\nCheck bcr_tsv path and confirm 01_BCR_Pipeline.R has been run.")
  }
  bcr <- data.table::fread(s$bcr_tsv, sep = "\t", data.table = FALSE, nThread = 1)
  message("  BCR TSV loaded: ", nrow(bcr), " rows")


  # ---- 3. SANITY CHECK: DUPLICATES PER CELL ----------------------------------
  heavy_rows <- dplyr::filter(bcr, locus == "IGH")
  dup_heavy  <- heavy_rows %>% dplyr::count(cell_id) %>% dplyr::filter(n > 1)

  if (nrow(dup_heavy) > 0) {
    warning("  ", nrow(dup_heavy), " cells have >1 heavy chain -- keeping highest UMI.",
            " This should have been resolved in 01_BCR_Pipeline.R.")
    heavy_rows <- heavy_rows %>%
      dplyr::group_by(cell_id) %>%
      dplyr::slice_max(umi_count, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
  } else {
    message("  Sanity check passed: no duplicate heavy chains")
  }

  light_rows <- dplyr::filter(bcr, locus %in% c("IGK", "IGL"))
  dup_light  <- light_rows %>% dplyr::count(cell_id) %>% dplyr::filter(n > 1)

  if (nrow(dup_light) > 0) {
    warning("  ", nrow(dup_light), " cells have >1 light chain -- keeping highest UMI.")
    light_rows <- light_rows %>%
      dplyr::group_by(cell_id) %>%
      dplyr::slice_max(umi_count, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
  } else {
    message("  Sanity check passed: no duplicate light chains")
  }


  # ---- 4. BUILD PER-CELL METADATA TABLE --------------------------------------
  heavy_meta <- heavy_rows %>%
    dplyr::transmute(
      cell_id,
      c_call      = c_call,
      clone_id    = clone_id,
      clone_count = clone_count,
      v_call      = v_call,
      j_call      = j_call,
      mu_freq_H   = mu_freq,
      umi_count_H = umi_count
    )

  light_meta <- light_rows %>%
    dplyr::transmute(
      cell_id,
      kappa_lambda = ifelse(locus == "IGK", "kappa", "lambda"),
      mu_freq_L    = mu_freq,
      umi_count_L  = umi_count
    )


  # ---- 5. CHECK CELL ID OVERLAP ----------------------------------------------
  b_cells     <- colnames(obj)
  bcr_cells   <- unique(heavy_meta$cell_id)
  n_matched   <- length(intersect(b_cells, bcr_cells))
  n_bcr_only  <- length(setdiff(bcr_cells, b_cells))
  n_gex_only  <- length(setdiff(b_cells, bcr_cells))
  pct_matched <- round(100 * n_matched / length(b_cells), 1)

  message("  BCR matching (B cells only):")
  message("    B cells in Seurat:          ", length(b_cells))
  message("    Cells with BCR detected:    ", length(bcr_cells))
  message("    Matched:                    ", n_matched,
          "  (", pct_matched, "% of B cells)")
  message("    BCR only (filtered in GEX): ", n_bcr_only)
  message("    B cells with no BCR:        ", n_gex_only)

  if (pct_matched < 30) {
    warning("Less than 30% of B cells matched to BCR data for ", s$name,
            ". Check that barcode suffixes are consistent between 01_BCR_Pipeline.R and 02_GEX_QC.R.")
  }


  # ---- 6. JOIN BCR METADATA ONTO SEURAT CELLS --------------------------------
  meta_add <- tibble::tibble(cell_id = b_cells) %>%
    dplyr::mutate(BCR = cell_id %in% bcr_cells) %>%
    dplyr::left_join(heavy_meta, by = "cell_id") %>%
    dplyr::left_join(light_meta, by = "cell_id")

  meta_df <- as.data.frame(meta_add)
  rownames(meta_df) <- meta_df$cell_id
  meta_df$cell_id   <- NULL

  obj <- AddMetaData(obj, metadata = meta_df)
  message("  BCR+ cells: ", sum(obj$BCR, na.rm = TRUE))


  # ---- 7. INSPECT: ISOTYPE AND TOP CLONES ------------------------------------
  bcr_pos <- subset(obj, subset = BCR == TRUE)

  message("  Isotype distribution (BCR+ cells with resolved isotype):")
  print(table(bcr_pos$c_call, useNA = "no"))

  message("  Top 10 largest clones:")
  top_clones <- bcr_pos@meta.data %>%
    dplyr::filter(!is.na(clone_id)) %>%
    dplyr::distinct(clone_id, clone_count) %>%
    dplyr::arrange(dplyr::desc(clone_count)) %>%
    head(10)
  print(top_clones)


  # ---- 8. SAVE ---------------------------------------------------------------
  saveRDS(obj, file = s$rds_out)
  message("  Saved: ", s$rds_out)

  summary_rows[[s$name]] <- data.frame(
    Sample      = s$name,
    B_cells     = length(b_cells),
    BCR_cells   = length(bcr_cells),
    Matched     = n_matched,
    Pct_matched = pct_matched,
    BCR_plus    = sum(obj$BCR, na.rm = TRUE)
  )
}


# =============================================================================
# FINAL SUMMARY TABLE
# =============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("SUMMARY ACROSS ALL SAMPLES")
message(paste(rep("=", 60), collapse = ""))
summary_df <- do.call(rbind, summary_rows)
rownames(summary_df) <- NULL
print(summary_df)
message("\nDone. All RDS files saved.")


# =============================================================================
# SUMMARY PLOTS
# =============================================================================

# Build a lookup: sample name -> rds_out path
rds_out_lookup <- setNames(
  sapply(samples, `[[`, "rds_out"),
  sapply(samples, `[[`, "name")
)

# ---- BCR+ vs BCR- per sample (percentage stacked bar) -----------------------
bcr_counts <- do.call(rbind, lapply(summary_df$Sample, function(sname) {
  obj <- readRDS(rds_out_lookup[[sname]])
  data.frame(
    Sample = sname,
    Status = ifelse(obj$BCR, "BCR+", "BCR-"),
    stringsAsFactors = FALSE
  )
}))

print(
  ggplot(bcr_counts, aes(x = Sample, fill = Status)) +
    geom_bar(position = "fill") +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(values = c("BCR+" = "#2E86AB", "BCR-" = "#D9D9D9")) +
    labs(title = "B cells with and without detected BCR",
         x = NULL, y = "Percentage of B cells", fill = NULL) +
    theme_classic(base_size = 13)
)

# ---- Isotype distribution per sample (% of BCR+ cells) ---------------------
isotype_counts <- do.call(rbind, lapply(summary_df$Sample, function(sname) {
  obj     <- readRDS(rds_out_lookup[[sname]])
  bcr_pos <- subset(obj, subset = BCR == TRUE)
  df <- as.data.frame(table(bcr_pos$c_call), stringsAsFactors = FALSE)
  colnames(df) <- c("Isotype", "Count")
  df$Pct    <- 100 * df$Count / sum(df$Count)
  df$Sample <- sname
  df
}))

isotype_order <- c("IGHM", "IGHD", "IGHA1", "IGHA2",
                   "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHE")
isotype_counts$Isotype <- factor(isotype_counts$Isotype, levels = isotype_order)

print(
  ggplot(isotype_counts, aes(x = Isotype, y = Pct, fill = Sample)) +
    geom_col(position = "dodge") +
    labs(title = "Isotype distribution per sample (BCR+ cells)",
         x = NULL, y = "% of BCR+ cells", fill = "Sample") +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)
