# =============================================================================
# 08_BCR_Visualization.R
# Join cell type annotations onto BCR data and generate clonality visualizations
# =============================================================================
#
# HOW TO USE:
#   1. Fill in the file paths and settings in the USER SETTINGS section below
#   2. Run the script (Ctrl+Shift+Enter in RStudio to source the whole file)
#
# This script is split into three parts that you can run sequentially or
# independently:
#
#   PART A -- Join cell type and compute duplicate counts
#     Loads the BCR TSV from step 01 and the annotated Seurat object from
#     step 07. Joins cell_type onto each BCR contig by barcode, and computes
#     VDJ duplicate counts (how many cells in the same clone share an identical
#     VDJ DNA sequence). Saves the annotated long-format TSV and a wide-format
#     Excel summary table.
#
#   PART B -- Donut plots
#     Generates clonality donut plots at three levels of granularity:
#     combined (whole sample), per isotype, and per cell type cluster.
#     Requires Part A to have been run first (bcr_annotated must be in workspace).
#
#   PART C -- Phylogenetic trees
#     Builds germline-rooted phylogenetic trees for the top expanded clones
#     using IQ-TREE 2 via the dowser package. Requires IQ-TREE 2 to be
#     installed (see note below). Saves trees as PDFs.
#     Requires Part A to have been run first (bcr_annotated must be in workspace).
#
# INPUT FILES:
#   - integrated_S7_annotated.rds    from step 07
#   - <SAMPLE_ID>_bcr_data.tsv       from step 01, one per sample
#   - BCR_viz_functions.R            companion functions file (same folder as
#                                    this script)
#
# OUTPUT:
#   - <SAMPLE_ID>_S8_bcr_annotated.tsv    long-format BCR with cell_type joined
#   - <SAMPLE_ID>_S8_final_table.xlsx     wide-format one-row-per-cell summary
#   - Plots/Part8/Trees/<SAMPLE_ID>_<CLONE_ID>_tree.pdf   tree PDFs
#
# NOTE -- IQ-TREE 2 (required for Part C only):
#   Download the binary for your OS from:
#     https://github.com/Cibiv/IQ-TREE/releases
#   Unzip and note the path to the executable, then set IQTREE_EXEC below.
#   If you add the IQ-TREE bin folder to your system PATH you can set
#   IQTREE_EXEC <- "iqtree2" and it will be found automatically.
#   Parts A and B do not require IQ-TREE and can be run without it.
# =============================================================================


# =============================================================================
# !! ONLY EDIT THIS SECTION !!
# =============================================================================

# Path to BCR_viz_functions.R -- update to wherever you saved it
BCR_VIZ_FUNCTIONS_FILE <- "C:/Users/YourName/Documents/BCR-GEX-Tutorial/BCR_viz_functions.R"

# Annotated Seurat object from step 07
RDS_IN <- "C:/Users/YourName/Documents/MyProject/RDS_Objects/integrated_S7_annotated.rds"

# BCR TSV files from step 01 -- one entry per sample
# Format: list(<sample_id> = <path_to_tsv>, ...)
BCR_TSV <- list(
  MySample1 = "C:/Users/YourName/Documents/MyProject/MySample1/Output/MySample1_bcr_data.tsv",
  MySample2 = "C:/Users/YourName/Documents/MyProject/MySample2/Output/MySample2_bcr_data.tsv"
)

# Output directory for annotated BCR files and Excel tables
OUT_DIR <- "C:/Users/YourName/Documents/MyProject/BCR_Data/Annotated"

# Output directory for plots
PLOT_DIR <- "C:/Users/YourName/Documents/MyProject/Plots/Part8"

# --- Part C: IQ-TREE settings -------------------------------------------------
# Path to IQ-TREE 2 executable. Use "iqtree2" if it is on your system PATH,
# or supply the full path:
#   Windows : IQTREE_EXEC <- "C:/tools/iqtree2/bin/iqtree2.exe"
#   macOS   : IQTREE_EXEC <- "/usr/local/bin/iqtree2"
IQTREE_EXEC <- "iqtree2"

# Number of trees to build per sample (top N clones by size)
TOP_N_TREES <- 5

# Minimum clone size (paired cells) required to attempt a tree
MIN_CLONE_SIZE <- 3

# Cell type color palette -- update labels and colors to match your cell_type
# annotations from step 07. Names must match the cell_type values exactly.
# Add or remove entries to match the clusters in your dataset.
cluster_colors <- c(
  "CellType_1" = "#F8766D",
  "CellType_2" = "#DB8E00",
  "CellType_3" = "#AEA200",
  "CellType_4" = "#64B200",
  "CellType_5" = "#00BD5C",
  "CellType_6" = "#00C1A7",
  "CellType_7" = "#00A6FF",
  "CellType_8" = "#B385FF",
  "CellType_9" = "#EF67EB",
  "CellType_10" = "#FF63B6"
  # Add or remove entries to match the clusters in your dataset
)

# =============================================================================
# END OF EDITABLE SECTION -- do not change anything below
# =============================================================================


# ---- Load packages ----------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(writexl)
  library(tibble)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(dowser)
  library(ggtree)
  library(ape)
  library(RColorBrewer)
  library(stringr)
})


# ---- Load BCR_viz_functions.R -----------------------------------------------
if (!file.exists(BCR_VIZ_FUNCTIONS_FILE)) {
  stop(
    "Cannot find BCR_viz_functions.R at:\n  ", BCR_VIZ_FUNCTIONS_FILE,
    "\nPlease update BCR_VIZ_FUNCTIONS_FILE at the top of this script."
  )
}
source(BCR_VIZ_FUNCTIONS_FILE)


# ---- Check all input files exist --------------------------------------------
to_check <- c(
  "Seurat RDS" = RDS_IN,
  setNames(unlist(BCR_TSV), paste0("BCR TSV [", names(BCR_TSV), "]"))
)
missing_files <- to_check[!file.exists(to_check)]
if (length(missing_files) > 0) {
  stop(
    "The following paths could not be found:\n  ",
    paste(names(missing_files), "->", missing_files, collapse = "\n  "),
    "\nPlease check the paths at the top of this script."
  )
}

dir.create(OUT_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# PART A -- JOIN CELL TYPE AND COMPUTE DUPLICATE COUNTS
# =============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("PART A -- Join cell type and compute duplicate counts")
message(paste(rep("=", 60), collapse = ""))


# ---- A1. Load Seurat object and extract cell type map -----------------------
message("\n>> A1: Loading annotated Seurat object...")

obj <- readRDS(RDS_IN)
message("Loaded: ", ncol(obj), " cells | ",
        length(unique(obj$cell_type)), " cell types")

cell_type_map <- obj@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::select(cell_id, cell_type, sample_id)

message("Cell type distribution:")
print(table(cell_type_map$cell_type))


# ---- A2. Join cell type onto BCR data per sample ----------------------------
message("\n>> A2: Joining cell type onto BCR data...")

bcr_annotated <- list()

for (sid in names(BCR_TSV)) {

  message("\n--- Processing: ", sid, " ---")

  # fread is faster than read_tsv for large files on OneDrive-backed paths
  bcr <- data.table::fread(BCR_TSV[[sid]], sep = "\t", data.table = FALSE)
  message("  BCR contigs loaded: ", nrow(bcr))

  # Join cell_type by barcode
  ct_this_sample <- dplyr::filter(cell_type_map, sample_id == sid) %>%
    dplyr::select(cell_id, cell_type)

  bcr <- dplyr::left_join(bcr, ct_this_sample, by = "cell_id")

  # Report match rate on heavy chains (one per cell)
  n_matched   <- sum(!is.na(bcr$cell_type[bcr$locus == "IGH"]))
  n_igh_total <- sum(bcr$locus == "IGH")
  pct_matched <- round(100 * n_matched / n_igh_total, 1)
  message("  Cell type matched: ", n_matched, " / ", n_igh_total,
          " heavy chain cells (", pct_matched, "%)")

  if (pct_matched < 50) {
    warning("  Less than 50% of heavy chain cells matched to a cell type. ",
            "Check that barcode suffixes are consistent between step 01 and step 02.")
  }

  # Compute VDJ duplicate counts within each clone.
  # For each chain, trim the sequence to V-start -> J-end (same coordinates
  # used by the tree builder) and count how many cells in the same clone share
  # an identical VDJ DNA sequence. This is stored in Duplicate_DNA_H / _L.
  bcr <- bcr %>%
    dplyr::mutate(
      VDJ_DNA = dplyr::if_else(
        !is.na(v_sequence_start) & !is.na(j_sequence_end) &
          v_sequence_start > 0   & j_sequence_end > 0,
        substr(sequence, v_sequence_start, j_sequence_end),
        sequence
      )
    ) %>%
    dplyr::group_by(clone_id, locus, VDJ_DNA) %>%
    dplyr::mutate(
      Duplicate_DNA_H = dplyr::if_else(locus == "IGH",             dplyr::n(), NA_integer_),
      Duplicate_DNA_L = dplyr::if_else(locus %in% c("IGK", "IGL"), dplyr::n(), NA_integer_)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-VDJ_DNA)

  bcr_annotated[[sid]] <- bcr
  message("  Done: ", nrow(bcr), " total contigs | ",
          length(unique(bcr$cell_id[bcr$locus == "IGH"])), " cells with heavy chain")
}


# ---- A3. Save long-format annotated TSV -------------------------------------
message("\n>> A3: Saving annotated long-format TSV files...")

for (sid in names(bcr_annotated)) {
  out_path <- file.path(OUT_DIR, paste0(sid, "_S8_bcr_annotated.tsv"))
  data.table::fwrite(bcr_annotated[[sid]], file = out_path, sep = "\t")
  message("Saved: ", out_path)
}


# ---- A4. Build and save wide-format Excel table -----------------------------
# One row per cell with heavy and light chain information side by side.
# Includes: cell barcode, clone ID, clone size, sample, cell type, isotype,
# light chain type, VDJ duplicate counts (heavy and light), SHM counts and
# frequencies, V/D/J gene calls, and full nucleotide sequences.
message("\n>> A4: Building and saving wide-format Excel tables...")

for (sid in names(bcr_annotated)) {

  bcr <- bcr_annotated[[sid]]

  heavy <- bcr %>%
    dplyr::filter(locus == "IGH") %>%
    dplyr::transmute(
      cell_id,
      clone_id,
      clone_count,
      sample_id,
      cell_type,
      isotype         = c_call,
      Duplicate_DNA_H,
      mu_count_H      = mu_count,
      mu_freq_H       = mu_freq,
      v_call_H        = v_call,
      d_call_H        = d_call,
      j_call_H        = j_call,
      sequence_H      = sequence
    )

  light <- bcr %>%
    dplyr::filter(locus %in% c("IGK", "IGL")) %>%
    dplyr::transmute(
      cell_id,
      light_chain = dplyr::case_when(
        locus == "IGK" ~ "Kappa",
        locus == "IGL" ~ "Lambda",
        TRUE           ~ locus
      ),
      Duplicate_DNA_L,
      mu_count_L = mu_count,
      mu_freq_L  = mu_freq,
      v_call_L   = v_call,
      j_call_L   = j_call,
      sequence_L = sequence
    )

  wide <- dplyr::left_join(heavy, light, by = "cell_id") %>%
    dplyr::arrange(dplyr::desc(clone_count))

  out_path <- file.path(OUT_DIR, paste0(sid, "_S8_final_table.xlsx"))
  writexl::write_xlsx(wide, out_path)
  message(sid, ": wide-format table saved (", nrow(wide), " cells) -- ", out_path)
}

message("\nPart A complete.")


# =============================================================================
# PART B -- DONUT PLOTS
# =============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("PART B -- Donut plots")
message(paste(rep("=", 60), collapse = ""))

for (sid in names(bcr_annotated)) {

  message("\n--- ", sid, " ---")

  # Combined donut: all clones across all isotypes
  message("Combined donut...")
  p <- plot_combined_donut(bcr_annotated[[sid]], sample_name = sid)
  print(p$all_clones)

  # Per-isotype donuts
  message("Per-isotype donuts...")
  iso_plots <- plot_isotype_donuts(bcr_annotated[[sid]], sample_name = sid)
  for (iso in names(iso_plots)) {
    print(iso_plots[[iso]]$all_clones)
  }

  # Per-cluster donuts (excludes NA cell_type cells automatically)
  message("Per-cluster donuts...")
  cl_plots <- plot_cluster_donuts(
    bcr_annotated[[sid]],
    sample_name   = sid,
    cell_type_col = "cell_type"
  )
  for (cl in names(cl_plots)) {
    print(cl_plots[[cl]]$all_clones)
  }
}

message("\nPart B complete.")


# =============================================================================
# PART C -- PHYLOGENETIC TREES
# =============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("PART C -- Phylogenetic trees")
message(paste(rep("=", 60), collapse = ""))

tree_dir <- file.path(PLOT_DIR, "Trees")
dir.create(tree_dir, recursive = TRUE, showWarnings = FALSE)

tree_plots_all <- list()

for (sid in names(bcr_annotated)) {

  message("\n--- Building trees for: ", sid, " ---")

  tree_plots_all[[sid]] <- build_bcr_trees(
    bcr_data       = bcr_annotated[[sid]],
    cluster_colors = cluster_colors,
    top_n          = TOP_N_TREES,
    cell_type_col  = "cell_type",
    exec           = IQTREE_EXEC,
    min_clone_size = MIN_CLONE_SIZE
  )

  # Print to viewer
  for (cid in names(tree_plots_all[[sid]])) {
    print(tree_plots_all[[sid]][[cid]])
  }

  # Save to PDF
  for (cid in names(tree_plots_all[[sid]])) {
    out_file <- file.path(tree_dir, paste0(cid, "_tree.pdf"))
    ggplot2::ggsave(out_file, tree_plots_all[[sid]][[cid]], width = 10, height = 8)
    message("Saved: ", out_file)
  }
}

message("\nPart C complete.")
message("\n", paste(rep("=", 60), collapse = ""))
message("Done. All outputs saved to: ", OUT_DIR)
message(paste(rep("=", 60), collapse = ""))
