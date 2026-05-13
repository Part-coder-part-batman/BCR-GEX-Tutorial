# =============================================================================
# 06_Integrate_Samples.R
# Merge per-sample Seurat objects and integrate across samples using RPCA
# =============================================================================
#
# HOW TO USE:
#   1. Fill in the sample definitions and settings below
#   2. Run the script (Ctrl+Shift+Enter in RStudio to source the whole file)
#
# INPUT FILES (one per sample):
#   - <SAMPLE_ID>_S5_soloProcessed.rds   produced by 05_Solo_GEX.R
#       Provide the full path as rds_in in the sample definitions below.
#
# OUTPUT:
#   - integrated_S6.rds
#       Saved to the full path you provide as rds_out below.
#
# WHAT THIS SCRIPT DOES:
#   Loads all per-sample objects, merges them, normalizes and finds HVGs
#   (excluding IG/TCR variable genes), runs an unintegrated PCA and UMAP
#   as a pre-integration batch check, then runs RPCA integration via
#   Seurat v5's IntegrateLayers(). After integration, runs clustering and
#   UMAP on the corrected embedding and generates a full panel of plots
#   for evaluating integration quality and inspecting the joint cluster
#   structure.
#
# NOTE — INTEGRATION METHOD:
#   This script uses RPCA by default. It is the right choice for same-platform
#   same-experiment data. For more divergent datasets (different conditions,
#   different tissues, many samples) consider CCA or Harmony. To switch,
#   change the IntegrateLayers() call:
#     RPCA    : method = RPCAIntegration,    new.reduction = "integrated.rpca"
#     CCA     : method = CCAIntegration,     new.reduction = "integrated.cca"
#     Harmony : method = HarmonyIntegration, new.reduction = "harmony"
#   Then update the reduction name in FindNeighbors() and RunUMAP() to match.
#
# NOTE — IG/TCR GENE EXCLUSION:
#   V(D)J variable gene transcripts are excluded from HVG selection for the
#   same reason as in step 05: their variability reflects clonal diversity,
#   not cell state. They are NOT removed from the object.
#
# NOTE — SEURAT VERSION:
#   Requires Seurat v5. Uses the IntegrateLayers() framework with per-sample
#   RNA layers. Not compatible with Seurat v4.
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
# rds_in : full path to the _S5_soloProcessed.rds file produced by 05_Solo_GEX.R
# The name field is used as the sample_id stored in object metadata and as the
# barcode prefix -- it must be unique across samples.

samples <- list(
  list(
    name   = "MySample1",
    rds_in = "C:/Users/YourName/Documents/MyProject/RDS_Objects/MySample1_S5_soloProcessed.rds"
  ),
  list(
    name   = "MySample2",
    rds_in = "C:/Users/YourName/Documents/MyProject/RDS_Objects/MySample2_S5_soloProcessed.rds"
  )
)

# --- Output path --------------------------------------------------------------
rds_out <- "C:/Users/YourName/Documents/MyProject/RDS_Objects/integrated_S6.rds"

# --- PCA settings -------------------------------------------------------------
# Variance threshold for automatic PC selection (0.8 = 80% of variance explained)
# The script selects the fewest PCs that together explain this fraction of
# variance, with a floor of 10 and a ceiling of 25.
var_expl <- 0.8

# --- Clustering resolution ----------------------------------------------------
# Higher values produce more clusters. 0.5 is a reasonable starting point.
# Adjust after inspecting the UMAP if needed.
cluster_resolution <- 0.5

# =============================================================================
# END OF MODIFIABLE SECTION — do not change anything below
# =============================================================================


# ---- IG/TCR variable gene exclusion patterns --------------------------------
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


# =============================================================================
# 1. LOAD
# =============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("Loading ", length(samples), " samples...")
message(paste(rep("=", 60), collapse = ""))

seurat_list  <- list()

for (s in samples) {
  if (!file.exists(s$rds_in)) {
    stop("Cannot find input RDS:\n  ", s$rds_in,
         "\nCheck rds_in for sample '", s$name,
         "' and confirm 05_Solo_GEX.R has been run.")
  }
  obj            <- readRDS(s$rds_in)
  obj$sample_id  <- s$name
  seurat_list[[s$name]] <- obj
  message("  Loaded: ", s$name, " — ", ncol(obj), " cells | ", nrow(obj), " genes")
}


# =============================================================================
# 2. MERGE
# =============================================================================
message("\nMerging samples...")

merged <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1]
)

message("Merged object: ", ncol(merged), " cells | ", nrow(merged), " genes")


# =============================================================================
# 3. PREPARE LAYERS FOR SEURAT V5 INTEGRATION
# =============================================================================
DefaultAssay(merged) <- "RNA"
merged[["RNA"]] <- JoinLayers(merged[["RNA"]])
merged[["RNA"]] <- split(merged[["RNA"]], f = merged$sample_id)
message("Layers after split: ", paste(Layers(merged), collapse = ", "))


# =============================================================================
# 4. NORMALIZE + HVGs
# =============================================================================
message("\nNormalizing and finding variable features...")

merged <- NormalizeData(merged, normalization.method = "LogNormalize",
                        scale.factor = 10000, verbose = FALSE)
merged <- FindVariableFeatures(merged, selection.method = "vst",
                               nfeatures = 2000, verbose = FALSE)

genes_to_exclude <- grep(variable_gene_regex, rownames(merged), value = TRUE)
n_excluded       <- length(intersect(genes_to_exclude, VariableFeatures(merged)))
hvg_filtered     <- setdiff(VariableFeatures(merged), genes_to_exclude)

if (length(hvg_filtered) < 50) {
  warning("Fewer than 50 HVGs after IG/TCR filtering — using unfiltered HVG list.")
  hvg_filtered <- VariableFeatures(merged)
}

VariableFeatures(merged) <- hvg_filtered
message("  IG/TCR variable genes excluded from HVG list: ", n_excluded)
message("  HVGs after filtering: ", length(hvg_filtered))


# =============================================================================
# 5. SCALE + PCA (UNINTEGRATED)
# =============================================================================
message("\nScaling and running PCA...")

merged <- ScaleData(merged, features = hvg_filtered, verbose = FALSE)
merged <- RunPCA(merged, features = hvg_filtered, npcs = 25, verbose = FALSE)

# Auto-select PCs
eigs    <- merged@reductions$pca@stdev^2
cum_var <- cumsum(eigs) / sum(eigs)
n_pc    <- which(cum_var >= var_expl)[1]
if (is.na(n_pc)) n_pc <- 25
n_pc    <- max(n_pc, 10)
n_pc    <- min(n_pc, 25)
message("  PCs selected: ", n_pc,
        " (explains ", round(100 * cum_var[n_pc], 1), "% of variance)")

print(
  ElbowPlot(merged, ndims = 25) +
    ggtitle("Elbow plot — merged unintegrated") +
    geom_vline(xintercept = n_pc, linetype = "dashed", color = "red") +
    labs(caption = paste0("Red dashed line = PC ", n_pc,
                          " (", round(100 * cum_var[n_pc], 1),
                          "% cumulative variance)"))
)


# =============================================================================
# 6. PRE-INTEGRATION UMAP (batch check)
# =============================================================================
message("\nRunning pre-integration UMAP...")

merged <- RunUMAP(merged, dims = 1:n_pc, reduction = "pca",
                  reduction.name = "umap.unintegrated", verbose = FALSE)

print(
  DimPlot(merged, reduction = "umap.unintegrated", group.by = "sample_id",
          pt.size = 0.8) +
    ggtitle("UMAP — Unintegrated (colored by sample)") +
    theme(plot.title = element_text(hjust = 0.5))
)


# =============================================================================
# 7. RPCA INTEGRATION
# =============================================================================
message("\nRunning RPCA integration...")

options(future.globals.maxSize = 2000 * 1024^2)

merged <- IntegrateLayers(
  object         = merged,
  method         = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction  = "integrated.rpca",
  verbose        = FALSE
)

merged[["RNA"]] <- JoinLayers(merged[["RNA"]])
message("  Integration complete. Layers rejoined.")


# =============================================================================
# 8. CLUSTER + UMAP (INTEGRATED)
# =============================================================================
message("\nClustering and running UMAP on integrated embedding...")

merged <- FindNeighbors(merged, dims = 1:n_pc,
                        reduction = "integrated.rpca", verbose = FALSE)
merged <- FindClusters(merged, resolution = cluster_resolution, verbose = FALSE)
merged <- RunUMAP(merged, dims = 1:n_pc,
                  reduction = "integrated.rpca", verbose = FALSE)

n_clusters <- length(unique(merged$seurat_clusters))
message("  Clusters found: ", n_clusters)


# =============================================================================
# 9. PLOTS
# =============================================================================
message("\nGenerating plots...")

# --- 9a. Integrated UMAP colored by sample ---
print(
  DimPlot(merged, reduction = "umap", group.by = "sample_id", pt.size = 0.8) +
    ggtitle("UMAP — Integrated (colored by sample)") +
    theme(plot.title = element_text(hjust = 0.5))
)

# --- 9b. Integrated UMAP colored by cluster ---
num_clusters <- table(merged$seurat_clusters)
print(
  DimPlot(merged, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.8) +
    scale_color_hue(labels = paste0(names(num_clusters), " (", num_clusters, ")")) +
    ggtitle(paste0("UMAP — Integrated clusters (res = ", cluster_resolution, ")")) +
    theme(legend.text = element_text(size = 9))
)

# --- 9c. Cluster composition by sample ---
comp_df <- merged@meta.data %>%
  dplyr::select(seurat_clusters, sample_id) %>%
  dplyr::group_by(seurat_clusters, sample_id) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop_last") %>%
  dplyr::mutate(pct = 100 * n / sum(n)) %>%
  dplyr::ungroup()

print(
  ggplot(comp_df, aes(x = seurat_clusters, y = pct, fill = sample_id)) +
    geom_bar(stat = "identity") +
    ylab("% of cells") + xlab("Cluster") +
    theme_classic(base_size = 12) +
    ggtitle("Cluster composition by sample")
)

# --- 9d. DotPlot: B cell marker panel ---
markers_present <- intersect(b_cell_markers, rownames(merged))
print(
  DotPlot(merged, features = markers_present,
          cols = c("gray90", "red"), dot.scale = 6) +
    RotatedAxis() +
    ggtitle("B cell marker panel — integrated clusters") +
    theme(axis.text.x = element_text(size = 8))
)

# --- 9e. FeaturePlots: key B cell markers ---
markers_for_feature <- intersect(
  c("CD79A", "IGHD", "IGHM", "CD27", "BCL6", "AICDA",
    "MKI67", "PRDM1", "XBP1", "FCRL4", "FCRL5"),
  rownames(merged)
)
if (length(markers_for_feature) > 0) {
  print(
    FeaturePlot(merged, reduction = "umap", features = markers_for_feature,
                order = TRUE, pt.size = 0.8) &
      scale_color_gradientn(colors = viridis::inferno(256))
  )
}

# --- 9f. BCR overlay plots (only if BCR metadata present) ---
# c_call is added by 04_Add_BCR.R — if it is absent, BCR overlay plots are skipped.
if ("c_call" %in% colnames(merged@meta.data)) {

  # Isotype UMAP
  obj_bcr <- subset(merged, subset = !is.na(c_call) & c_call != "")
  if (ncol(obj_bcr) > 0) {
    obj_bcr$c_call <- factor(obj_bcr$c_call, levels = isotype_order)
    present_iso    <- intersect(isotype_order, levels(droplevels(obj_bcr$c_call)))
    print(
      DimPlot(obj_bcr, group.by = "c_call", label = FALSE, pt.size = 0.8,
              cols = isotype_colors[present_iso],
              order = rev(present_iso)) +
        ggtitle("Integrated UMAP — Isotype (BCR+ cells)")
    )
  }

  # SHM frequency
  if ("mu_freq_H" %in% colnames(merged@meta.data)) {
    print(
      FeaturePlot(merged, features = "mu_freq_H", reduction = "umap",
                  order = TRUE, pt.size = 0.8) +
        scale_color_gradientn(colors = c("#FFFFCC", "#41B6C4", "#0C2C84"),
                              na.value = "grey90") +
        ggtitle("Integrated UMAP — SHM frequency, heavy chain")
    )
  }

  # Clone size
  if ("clone_count" %in% colnames(merged@meta.data)) {
    print(
      FeaturePlot(merged, features = "clone_count", reduction = "umap",
                  order = TRUE, pt.size = 0.8) +
        scale_color_gradientn(colors = c("#FFFFCC", "#FD8D3C", "#800026"),
                              na.value = "grey90") +
        ggtitle("Integrated UMAP — Clone size")
    )
  }

  # Isotype per cluster barplot
  iso_df <- merged@meta.data %>%
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
      ggtitle("Isotype composition per integrated cluster")
  )

} else {
  message("  No BCR metadata found — skipping BCR overlay plots.")
  message("  Run 04_Add_BCR.R first if you want BCR plots.")
}


# =============================================================================
# 10. SAVE
# =============================================================================
saveRDS(merged, file = rds_out)
message("\n", paste(rep("=", 60), collapse = ""))
message("Done.")
message("Saved: ", rds_out)
message("Cells: ", ncol(merged), " | Clusters: ", n_clusters)
message(paste(rep("=", 60), collapse = ""))
