# =============================================================================
# 05_Solo_GEX.R
# Per-sample normalization, clustering, UMAP, and BCR overlay
# =============================================================================
#
# HOW TO USE:
#   1. Fill in the sample definitions and settings below
#   2. Run the script (Ctrl+Shift+Enter in RStudio to source the whole file)
#
# INPUT FILES (one per sample):
#   - <SAMPLE_ID>_S4_postBCR.rds   produced by 04_Add_BCR.R
#       Provide the full path as rds_in in the sample definitions below.
#
# OUTPUT:
#   - <SAMPLE_ID>_S5_soloProcessed.rds
#       Saved to the full path you provide as rds_out in the sample definitions.
#
# WHAT THIS SCRIPT DOES:
#   For each sample: normalizes, finds HVGs (excluding IG/TCR variable genes),
#   runs PCA and UMAP, clusters at the specified resolution, scores cell cycle,
#   and generates a full panel of QC and BCR overlay plots.
#
#   This step is pre-integration. Each sample is processed independently so
#   you can inspect the data before merging samples in step 06.
#
# NOTE — IG/TCR GENE EXCLUSION:
#   V(D)J variable gene transcripts are excluded from HVG selection because
#   their variability reflects clonal diversity (which V gene a cell expresses)
#   rather than cell state. They are NOT removed from the object -- they remain
#   available for feature plots and differential expression. They are only
#   excluded from the features used for PCA and UMAP.
#
# NOTE — CELL CYCLE:
#   Cell cycle phase is scored and stored as metadata but is NOT regressed out.
#   In B cell data, proliferating cells are a biologically meaningful population
#   (typically dark zone GC B cells) and should remain visible in the UMAP.
#
# SEURAT VERSION NOTE:
#   Requires Seurat v5. Before NormalizeData() runs, Seurat v5 may print:
#   "Default search for 'data' layer in 'RNA' assay yielded no results"
#   This is harmless and expected -- the data layer is created by NormalizeData().
# =============================================================================


# ---- 0. PACKAGES ------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(viridis)
  library(scales)
  library(patchwork)
})


# =============================================================================
# !! ONLY MODIFY THIS SECTION !!
# =============================================================================

# --- Sample definitions -------------------------------------------------------
# Add one entry per sample.
# rds_in  : full path to the _S4_postBCR.rds file produced by 04_Add_BCR.R
# rds_out : full path where the _S5_soloProcessed.rds file will be saved

samples <- list(
  list(
    name    = "MySample1",
    rds_in  = "C:/Users/YourName/Documents/MyProject/RDS_Objects/MySample1_S4_postBCR.rds",
    rds_out = "C:/Users/YourName/Documents/MyProject/RDS_Objects/MySample1_S5_soloProcessed.rds"
  ),
  list(
    name    = "MySample2",
    rds_in  = "C:/Users/YourName/Documents/MyProject/RDS_Objects/MySample2_S4_postBCR.rds",
    rds_out = "C:/Users/YourName/Documents/MyProject/RDS_Objects/MySample2_S5_soloProcessed.rds"
  )
)

# --- PCA settings -------------------------------------------------------------
# Variance threshold for automatic PC selection (0.8 = 80% of variance explained)
# The script selects the fewest PCs that together explain this fraction of
# variance, with a floor of 10 and a ceiling of 25.
var_expl <- 0.8

# --- Clustering resolution ----------------------------------------------------
# Higher values produce more clusters. 0.5 is a reasonable starting point
# for a B cell dataset. Adjust after inspecting the UMAP if needed.
cluster_resolution <- 0.5

# =============================================================================
# END OF MODIFIABLE SECTION — do not change anything below
# =============================================================================


# ---- IG/TCR variable gene exclusion patterns --------------------------------
# These match V, D, and J gene transcripts for IG heavy, kappa, lambda, and
# TCR alpha/beta chains. They are removed from the HVG list before PCA.
# The ^IGHD pattern catches IGHD diversity genes (e.g. IGHD3-10) but NOT
# the IgD isotype transcript (IGHD with no number suffix), which is kept.
variable_gene_patterns <- c(
  "^IGHV", "^IGKV", "^IGKJ", "^IGLV", "^IGLJ",
  "^IGHD", "^IGHJ",
  "^TRAV", "^TRBV", "^TRAJ", "^TRBJ", "^TRBD"
)
variable_gene_regex <- paste(variable_gene_patterns, collapse = "|")


# ---- Isotype color scheme ---------------------------------------------------
isotype_order <- c("IGHD", "IGHM", "IGHA1", "IGHA2",
                   "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHE")

isotype_colors <- c(
  "IGHD"  = "#e41a1c",
  "IGHM"  = "#377eb8",
  "IGHA1" = "#a65628",
  "IGHA2" = "#f781bf",
  "IGHG1" = "#4daf4a",
  "IGHG2" = "#984ea3",
  "IGHG3" = "#ff7f00",
  "IGHG4" = "#ffff33",
  "IGHE"  = "#66c2a5"
)


# ---- B cell marker panel ----------------------------------------------------
b_cell_markers <- c(
  # B cell identity
  "CD79A", "CD19", "MS4A1",
  # Naive
  "IGHD", "IGHM", "IL4R", "TCL1A",
  # Memory
  "CD27", "CD24", "TNFRSF13B",
  # Germinal center
  "BCL6", "AICDA", "MEF2B",
  # Plasmablast / plasma cell
  "JCHAIN", "MZB1", "XBP1",
  # Proliferating
  "MKI67", "UBE2C", "PCNA",
  # Early activation
  "CD83", "MYC", "CCND2", "TRAF4", "MIR155HG",
  # Interferon-stimulated
  "IFI6", "IFI44", "IFIT1", "IFIT2", "IFIT3",
  # Atypical memory / ABCs
  "FCRL4", "FCRL5", "ZBTB32", "ITGAX", "ZEB2",
  # Exhaustion / regulatory
  "HAVCR1", "HAVCR2", "TIGIT", "LAG3", "IL10"
)


# ---- Cell cycle gene sets ---------------------------------------------------
s.genes   <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes


# =============================================================================
# LOOP OVER SAMPLES
# =============================================================================

for (s in samples) {

  message("\n", paste(rep("=", 60), collapse = ""))
  message("Processing sample: ", s$name)
  message(paste(rep("=", 60), collapse = ""))


  # ---- 1. LOAD ---------------------------------------------------------------
  if (!file.exists(s$rds_in)) {
    stop("Cannot find input RDS:\n  ", s$rds_in,
         "\nCheck rds_in for sample '", s$name,
         "' and confirm 04_Add_BCR.R has been run.")
  }
  obj <- readRDS(s$rds_in)
  DefaultAssay(obj) <- "RNA"
  message("  Loaded: ", ncol(obj), " cells | ", nrow(obj), " genes")


  # ---- 2. NORMALIZE ----------------------------------------------------------
  message("  Normalizing...")
  obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                       scale.factor = 10000, verbose = FALSE)


  # ---- 3. FIND HVGs AND EXCLUDE IG/TCR VARIABLE GENES -----------------------
  message("  Finding variable features...")
  obj <- FindVariableFeatures(obj, selection.method = "vst",
                              nfeatures = 2000, verbose = FALSE)

  genes_to_exclude <- grep(variable_gene_regex, rownames(obj), value = TRUE)
  n_excluded       <- length(intersect(genes_to_exclude, VariableFeatures(obj)))
  message("  Excluding ", n_excluded, " IG/TCR variable genes from HVG list")

  hvg_filtered <- setdiff(VariableFeatures(obj), genes_to_exclude)
  if (length(hvg_filtered) < 50) {
    warning("  Fewer than 50 HVGs remain after filtering. Using unfiltered HVG list.")
    hvg_filtered <- VariableFeatures(obj)
  }
  VariableFeatures(obj) <- hvg_filtered
  message("  HVGs after filtering: ", length(hvg_filtered))


  # ---- 4. SCALE AND PCA ------------------------------------------------------
  message("  Scaling and running PCA...")
  obj <- ScaleData(obj, features = hvg_filtered, verbose = FALSE)
  obj <- RunPCA(obj, features = hvg_filtered, npcs = 25, verbose = FALSE)

  # Auto-select PCs: fewest that explain >= var_expl of variance
  eigs    <- obj@reductions$pca@stdev^2
  cum_var <- cumsum(eigs) / sum(eigs)
  n_pc    <- which(cum_var >= var_expl)[1]
  if (is.na(n_pc)) n_pc <- 25
  n_pc <- max(n_pc, 10)
  n_pc <- min(n_pc, 25)
  message("  PCs selected: ", n_pc,
          " (explains ", round(100 * cum_var[n_pc], 1), "% of variance)")

  print(
    ElbowPlot(obj, ndims = 25) +
      ggtitle(paste0(s$name, " — Elbow plot")) +
      geom_vline(xintercept = n_pc, linetype = "dashed", color = "red") +
      labs(caption = paste0("Red dashed line = PC ", n_pc,
                            " (", round(100 * cum_var[n_pc], 1),
                            "% cumulative variance)"))
  )


  # ---- 5. CLUSTER AND UMAP ---------------------------------------------------
  message("  Clustering and running UMAP...")
  obj <- FindNeighbors(obj, dims = 1:n_pc, verbose = FALSE)
  obj <- FindClusters(obj, resolution = cluster_resolution, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:n_pc, verbose = FALSE)

  n_clusters <- length(unique(obj$seurat_clusters))
  message("  Clusters found: ", n_clusters)


  # ---- 6. CELL CYCLE SCORING -------------------------------------------------
  message("  Scoring cell cycle...")
  obj <- CellCycleScoring(obj, s.features = s.genes,
                          g2m.features = g2m.genes, set.ident = FALSE)
  obj$Phase <- factor(obj$Phase, levels = c("G1", "S", "G2M"))


  # ---- 7. PLOTS --------------------------------------------------------------
  message("  Generating plots...")

  # --- 7a. Cluster UMAP ---
  num_clusters <- table(obj$seurat_clusters)
  print(
    DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1.2) +
      scale_color_hue(labels = paste0(names(num_clusters),
                                      " (", num_clusters, ")")) +
      ggtitle(paste0(s$name, " — Clusters (res = ", cluster_resolution, ")")) +
      theme(legend.text = element_text(size = 9))
  )

  # --- 7b. QC metrics on UMAP ---
  print(
    FeaturePlot(obj, reduction = "umap",
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                order = TRUE, pt.size = 1.2) +
      patchwork::plot_annotation(
        title = paste0(s$name, " — QC metrics on UMAP"))
  )

  # --- 7c. DotPlot: B cell marker panel ---
  markers_present <- intersect(b_cell_markers, rownames(obj))
  print(
    DotPlot(obj, features = markers_present, cols = c("gray90", "red"),
            dot.scale = 6) +
      RotatedAxis() +
      ggtitle(paste0(s$name, " — B cell marker panel")) +
      theme(axis.text.x = element_text(size = 8))
  )

  # --- 7d. FeaturePlots: key B cell markers ---
  markers_for_feature <- intersect(
    c("CD79A", "IGHD", "IGHM", "CD27", "BCL6", "AICDA",
      "MKI67", "PRDM1", "XBP1", "FCRL4", "FCRL5"),
    rownames(obj)
  )
  if (length(markers_for_feature) > 0) {
    print(
      FeaturePlot(obj, reduction = "umap", features = markers_for_feature,
                  order = TRUE, pt.size = 1.2) &
        scale_color_gradientn(colors = viridis::inferno(256))
    )
  }

  # --- 7e. Cell cycle UMAP ---
  Idents(obj) <- "Phase"
  print(
    DimPlot(obj, reduction = "umap", label = FALSE, pt.size = 1.2,
            cols = c("G1" = "#AAAAAA", "S" = "#E69F00", "G2M" = "#CC79A7")) +
      ggtitle(paste0(s$name, " — Cell cycle phase"))
  )
  Idents(obj) <- "seurat_clusters"

  # --- 7f. Cell cycle composition per cluster ---
  cc_df <- obj@meta.data %>%
    dplyr::select(Phase, seurat_clusters) %>%
    dplyr::group_by(seurat_clusters, Phase) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop_last") %>%
    dplyr::mutate(pct = 100 * n / sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Phase = factor(Phase, levels = c("G1", "S", "G2M")))

  print(
    ggplot(cc_df, aes(x = seurat_clusters, y = pct, fill = Phase)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("G1" = "#AAAAAA", "S" = "#E69F00",
                                   "G2M" = "#CC79A7")) +
      ylab("% of cells") + xlab("Cluster") +
      theme_classic(base_size = 12) +
      ggtitle(paste0(s$name, " — Cell cycle composition per cluster"))
  )

  # --- 7g. BCR overlay plots (only if BCR metadata present) ---
  if ("BCR" %in% colnames(obj@meta.data)) {

    # Isotype UMAP (BCR+ cells with resolved isotype only)
    obj_bcr <- subset(obj, subset = !is.na(c_call) & c_call != "")
    if (ncol(obj_bcr) > 0) {
      obj_bcr$c_call <- factor(obj_bcr$c_call, levels = isotype_order)
      present_iso    <- intersect(isotype_order, levels(droplevels(obj_bcr$c_call)))
      print(
        DimPlot(obj_bcr, group.by = "c_call", label = FALSE, pt.size = 1.2,
                cols = isotype_colors[present_iso],
                order = rev(present_iso)) +
          ggtitle(paste0(s$name, " — Isotype (BCR+ cells)"))
      )
    }

    # SHM frequency UMAP (heavy chain)
    if ("mu_freq_H" %in% colnames(obj@meta.data)) {
      print(
        FeaturePlot(obj, features = "mu_freq_H", reduction = "umap",
                    order = TRUE, pt.size = 1.2) +
          scale_color_gradientn(colors = c("#FFFFCC", "#41B6C4", "#0C2C84"),
                                na.value = "grey90") +
          ggtitle(paste0(s$name, " — SHM frequency, heavy chain"))
      )
    }

    # Clone size UMAP
    if ("clone_count" %in% colnames(obj@meta.data)) {
      print(
        FeaturePlot(obj, features = "clone_count", reduction = "umap",
                    order = TRUE, pt.size = 1.2) +
          scale_color_gradientn(colors = c("#FFFFCC", "#FD8D3C", "#800026"),
                                na.value = "grey90") +
          ggtitle(paste0(s$name, " — Clone size"))
      )
    }

    # Isotype composition per cluster (stacked bar)
    if ("c_call" %in% colnames(obj@meta.data)) {
      iso_df <- obj@meta.data %>%
        dplyr::filter(!is.na(c_call), c_call != "", !is.na(seurat_clusters)) %>%
        dplyr::mutate(c_call          = factor(c_call, levels = isotype_order),
                      seurat_clusters = as.factor(seurat_clusters))
      present_iso2 <- intersect(isotype_order, levels(droplevels(iso_df$c_call)))
      print(
        ggplot(iso_df, aes(x = seurat_clusters, y = 1, fill = c_call)) +
          geom_col(position = "fill") +
          scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
          scale_fill_manual(values = isotype_colors[present_iso2],
                            breaks = present_iso2) +
          ylab("Percentage of BCR+ cells") + xlab("Cluster") +
          theme_classic(base_size = 12) +
          ggtitle(paste0(s$name, " — Isotype composition per cluster"))
      )
    }

  } else {
    message("  No BCR metadata found — skipping BCR overlay plots.")
    message("  Run 04_Add_BCR.R first if you want BCR plots.")
  }

  # Reset idents to clusters
  Idents(obj) <- "seurat_clusters"


  # ---- 8. SAVE ---------------------------------------------------------------
  saveRDS(obj, file = s$rds_out)
  message("  Saved: ", s$rds_out)
}

message("\n", paste(rep("=", 60), collapse = ""))
message("Done. All samples processed.")
message(paste(rep("=", 60), collapse = ""))
