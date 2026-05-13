# =============================================================================
# 03_Doublet_Removal.R
# B cell enrichment and marker-based doublet removal
# =============================================================================
#
# HOW TO USE:
#   1. Fill in the file paths in the USER SETTINGS section below
#   2. Run the script (Ctrl+Shift+Enter in RStudio to source the whole file)
#
# INPUT FILES (one per sample):
#   - <SAMPLE_ID>_S2_postQC.rds   produced by 02_GEX_QC.R
#       Provide the full path as rds_in in the sample definitions below.
#
# OUTPUT:
#   - <SAMPLE_ID>_S3_postDoublet.rds
#       Saved to the full path you provide as rds_out in the sample definitions.
#
# WHAT THIS SCRIPT DOES:
#   Step 1 — Keeps only B cells by filtering for CD79A > 0 in the counts layer
#   Step 2 — Removes likely doublets: CD79A+ cells that also express a non-B
#             cell marker (T cell, myeloid, or NK cell markers)
#
# SEURAT VERSION NOTE:
#   Requires Seurat v5.
#
# MARKER PANEL NOTE:
#   The default non-B markers (CD3E, CD3D, CD4, CD8A, LYZ, CST3, KLRD1,
#   KLRB1, CD247) work well for lymphoid and tumor tissue. Depending on your
#   tissue context you may want to add or substitute markers -- for example,
#   epithelial markers for mucosal tissue.
# =============================================================================


# ---- 0. PACKAGES ------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})


# =============================================================================
# !! ONLY MODIFY THIS SECTION !!
# =============================================================================

# --- Sample definitions -------------------------------------------------------
# Add one entry per sample.
# rds_in  : full path to the _S2_postQC.rds file produced by 02_GEX_QC.R
# rds_out : full path where the _S3_postDoublet.rds file will be saved

samples <- list(
  list(
    name    = "MySample1",
    rds_in  = "C:/Users/YourName/Documents/MyProject/RDS_Objects/MySample1_S2_postQC.rds",
    rds_out = "C:/Users/YourName/Documents/MyProject/RDS_Objects/MySample1_S3_postDoublet.rds"
  ),
  list(
    name    = "MySample2",
    rds_in  = "C:/Users/YourName/Documents/MyProject/RDS_Objects/MySample2_S2_postQC.rds",
    rds_out = "C:/Users/YourName/Documents/MyProject/RDS_Objects/MySample2_S3_postDoublet.rds"
  )
)

# =============================================================================
# END OF MODIFIABLE SECTION — do not change anything below
# =============================================================================


# ---- Helper: marker-based doublet detection ---------------------------------
# Returns barcodes of cells that are CD79A+ AND express any non-B cell marker
get_doublet_cells <- function(seurat_obj) {

  non_b_markers <- c("CD3E", "CD3D", "CD4", "CD8A",   # T cells
                     "LYZ",  "CST3",                   # Monocytes / myeloid
                     "KLRD1", "KLRB1", "CD247")        # NK cells

  cd79a <- tryCatch(
    FetchData(seurat_obj, vars = "CD79A", layer = "counts")[[1]],
    error = function(e) return(NULL)
  )
  if (is.null(cd79a)) return(NULL)

  doublet_cells <- lapply(non_b_markers, function(marker) {
    tryCatch({
      marker_expr <- FetchData(seurat_obj, vars = marker, layer = "counts")[[1]]
      rownames(seurat_obj@meta.data)[cd79a > 0 & marker_expr > 0]
    }, error = function(e) return(NULL))
  })

  unique(unlist(doublet_cells))
}


# =============================================================================
# LOAD SAMPLES
# =============================================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("Seurat version: ", as.character(packageVersion("Seurat")))
message("Loading ", length(samples), " sample(s)...")
message(paste(rep("=", 60), collapse = ""))

sample_list <- lapply(samples, function(s) {
  if (!file.exists(s$rds_in)) {
    stop("Cannot find input file:\n  ", s$rds_in,
         "\nCheck rds_in for sample '", s$name, "' and confirm 02_GEX_QC.R has been run.")
  }
  message("  Loaded: ", s$name)
  readRDS(s$rds_in)
})
names(sample_list) <- sapply(samples, `[[`, "name")


# =============================================================================
# LOOP OVER SAMPLES
# =============================================================================

summary_rows <- list()

for (s in samples) {

  message("\n", paste(rep("=", 60), collapse = ""))
  message("Processing sample: ", s$name)
  message(paste(rep("=", 60), collapse = ""))

  obj     <- sample_list[[s$name]]
  n_total <- ncol(obj)
  message("  Cells in (post-QC):          ", n_total)

  # ---- 1. KEEP ONLY B CELLS (CD79A > 0) -------------------------------------
  b_cells <- tryCatch({
    cd79a <- FetchData(obj, vars = "CD79A", layer = "counts")[[1]]
    rownames(obj@meta.data)[cd79a > 0]
  }, error = function(e) {
    message("  WARNING: CD79A not found in this object — skipping B cell filter.")
    return(colnames(obj))
  })

  obj      <- subset(obj, cells = b_cells)
  n_bcells <- ncol(obj)
  message("  B cells (CD79A+):            ", n_bcells,
          "  (", round(100 * n_bcells / n_total, 1), "% of input)")

  # ---- 2. IDENTIFY AND REMOVE DOUBLETS --------------------------------------
  doublets   <- get_doublet_cells(obj)
  n_doublets <- length(doublets)
  message("  Doublets detected:           ", n_doublets)

  if (n_doublets > 0) {
    obj <- subset(obj, cells = setdiff(colnames(obj), doublets))
  }

  n_final <- ncol(obj)
  message("  Cells out (post-doublet):    ", n_final,
          "  (", round(100 * n_final / n_bcells, 1), "% of B cells retained)")

  # ---- 3. SAVE --------------------------------------------------------------
  saveRDS(obj, file = s$rds_out)
  message("  Saved: ", s$rds_out)

  summary_rows[[s$name]] <- data.frame(
    Sample           = s$name,
    Cells_in         = n_total,
    B_cells_CD79Apos = n_bcells,
    Pct_B_cells      = round(100 * n_bcells / n_total, 1),
    Doublets_removed = n_doublets,
    Cells_out        = n_final,
    Pct_retained     = round(100 * n_final / n_bcells, 1)
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
