# =============================================================================
# 02_GEX_QC.R
# GEX quality control and filtering
# =============================================================================
#
# HOW TO USE:
#   1. Fill in the file paths and thresholds in the USER SETTINGS section below
#   2. Run the script (Ctrl+Shift+Enter in RStudio to source the whole file)
#
# INPUT FILES (one per sample):
#   - sample_filtered_feature_bc_matrix.h5
#       Found at: per_sample_outs > <sample> > count
#       Alternative: filtered_feature_bc_matrix/ folder (set h5 = NULL and
#       provide folder path instead -- see sample definitions below)
#
# OUTPUT:
#   - <SAMPLE_ID>_S2_postQC.rds   saved in <RDS_DIR>
#
# Run this script once. All samples defined below are processed in a loop.
#
# SEURAT VERSION NOTE:
#   This script uses Seurat v5. Before NormalizeData() is run, Seurat v5 may
#   print: "Default search for 'data' layer in 'RNA' assay yielded no results"
#   This is harmless and expected at this stage.
# =============================================================================


# ---- 0. PACKAGES ------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
})


# =============================================================================
# !! ONLY MODIFY THIS SECTION !!
# =============================================================================

# --- Sample definitions -------------------------------------------------------
# Add one entry per sample.
# Set h5 to the full path of the .h5 file, OR set h5 = NULL and provide
# a folder path to filtered_feature_bc_matrix/ instead.

samples <- list(
  list(
    name   = "MySample1",
    h5     = "C:/Users/YourName/Documents/MySample1/count/sample_filtered_feature_bc_matrix.h5",
    folder = NULL
  ),
  list(
    name   = "MySample2",
    h5     = "C:/Users/YourName/Documents/MySample2/count/sample_filtered_feature_bc_matrix.h5",
    folder = NULL
  )
)

# --- Output directory ---------------------------------------------------------
# All RDS files are saved here (flat folder, no subfolders)
RDS_DIR <- "C:/Users/YourName/Documents/MyProject/RDS_Objects"

# --- Organism -----------------------------------------------------------------
mito_pattern <- "^MT-"   # "^MT-" for human  |  "^mt-" for mouse

# --- Gene-level filter --------------------------------------------------------
# Keep only genes detected in at least this many cells
min_cells_gene <- 5

# --- QC thresholds ------------------------------------------------------------
# Inspect the violin plots produced by this script before deciding on values.
# These thresholds are applied to ALL samples.
#
# Typical starting points for human lymphoid / tumor tissue:
#   min_features : 200-500    (remove empty droplets)
#   max_features : 5000-8000  (remove doublets)
#   min_UMI      : 500-1000
#   max_UMI      : 20000-50000
#   filter_mt    : 10-20%     (remove damaged / dying cells)

min_features <- 250
max_features <- 6500
min_UMI      <- 500
max_UMI      <- 35000
filter_mt    <- 10

# =============================================================================
# END OF MODIFIABLE SECTION — do not change anything below
# =============================================================================


if (!dir.exists(RDS_DIR)) dir.create(RDS_DIR, recursive = TRUE)

filter_summary <- list()


# =============================================================================
# LOOP OVER SAMPLES
# =============================================================================

for (s in samples) {

  sample_name <- s$name
  message("\n", paste(rep("=", 60), collapse = ""))
  message("Processing sample: ", sample_name)
  message(paste(rep("=", 60), collapse = ""))

  # ---- 1. READ COUNTS -------------------------------------------------------
  message("Reading Cell Ranger counts...")

  if (!is.null(s$h5)) {
    if (!file.exists(s$h5)) stop("H5 file not found:\n  ", s$h5)
    counts <- Read10X_h5(s$h5)
  } else if (!is.null(s$folder)) {
    if (!dir.exists(s$folder)) stop("Matrix folder not found:\n  ", s$folder)
    counts <- Read10X(data.dir = s$folder)
  } else {
    stop("Sample '", sample_name, "': set either h5 or folder path.")
  }

  # Handle multi-assay outputs (e.g. Feature Barcoding)
  if (is.list(counts)) counts <- counts[["Gene Expression"]]


  # ---- 2. CREATE SEURAT OBJECT ----------------------------------------------
  obj <- CreateSeuratObject(
    counts       = counts,
    project      = sample_name,
    min.cells    = min_cells_gene,
    min.features = 0
  )

  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mito_pattern)

  # Suffix barcodes with sample name to prevent collisions during integration
  obj <- RenameCells(obj, new.names = paste0(colnames(obj), "_", sample_name))

  message("Seurat object created: ", ncol(obj), " cells | ", nrow(obj), " genes")


  # ---- 3. PRE-FILTER QC PLOTS -----------------------------------------------
  message("Plotting QC metrics BEFORE filtering...")

  print(
    VlnPlot(obj,
            features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
            pt.size  = 0.1,
            ncol     = 3) &
      ggplot2::labs(title = paste0(sample_name, " — QC metrics (pre-filter)"))
  )

  print(
    FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
      ggplot2::ggtitle(paste0(sample_name, " — UMI vs mt%"))
  )

  print(
    FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
      ggplot2::ggtitle(paste0(sample_name, " — UMI vs genes detected"))
  )

  pre_stats <- data.frame(
    Cells         = ncol(obj),
    Median_genes  = median(obj$nFeature_RNA),
    Median_UMI    = median(obj$nCount_RNA),
    Median_mt_pct = round(median(obj$percent.mt), 2),
    Max_mt_pct    = round(max(obj$percent.mt), 2),
    Cells_mt_gt10 = sum(obj$percent.mt > filter_mt),
    Cells_gene_hi = sum(obj$nFeature_RNA > max_features),
    Cells_gene_lo = sum(obj$nFeature_RNA < min_features),
    Cells_UMI_hi  = sum(obj$nCount_RNA > max_UMI),
    Cells_UMI_lo  = sum(obj$nCount_RNA < min_UMI)
  )
  message("--- Pre-filter summary ---")
  print(pre_stats)


  # ---- 4. APPLY QC FILTERS --------------------------------------------------
  cells_before <- ncol(obj)

  obj <- subset(
    obj,
    subset = nFeature_RNA >= min_features &
             nFeature_RNA <= max_features &
             nCount_RNA   >= min_UMI      &
             nCount_RNA   <= max_UMI      &
             percent.mt   <= filter_mt
  )

  cells_after <- ncol(obj)
  message(sprintf(
    "Cells kept: %d / %d  (%.1f%% retained, %d removed)",
    cells_after, cells_before,
    100 * cells_after / cells_before,
    cells_before - cells_after
  ))

  filter_summary[[sample_name]] <- data.frame(
    Sample        = sample_name,
    Cells_before  = cells_before,
    Cells_after   = cells_after,
    Cells_removed = cells_before - cells_after,
    Pct_retained  = round(100 * cells_after / cells_before, 1)
  )


  # ---- 5. POST-FILTER QC PLOTS ----------------------------------------------
  message("Plotting QC metrics AFTER filtering...")

  print(
    VlnPlot(obj,
            features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
            pt.size  = 0.1,
            ncol     = 3) &
      ggplot2::labs(title = paste0(sample_name, " — QC metrics (post-filter)"))
  )


  # ---- 6. SAVE --------------------------------------------------------------
  out_path <- file.path(RDS_DIR, paste0(sample_name, "_S2_postQC.rds"))
  saveRDS(obj, file = out_path)
  message("Saved: ", out_path)
}


# =============================================================================
# FINAL SUMMARY TABLE
# =============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("SUMMARY ACROSS ALL SAMPLES")
message(paste(rep("=", 60), collapse = ""))
print(do.call(rbind, filter_summary))
message("\nDone. All RDS files saved to: ", RDS_DIR)
