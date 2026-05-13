# =============================================================================
# 07_Annotate_Clusters.R
# Annotate integrated clusters with cell type labels
# =============================================================================
#
# HOW TO USE:
#   This script is split into three parts that you run sequentially:
#
#   PART A -- Annotation plots
#     Run this first. It generates all the plots you need to decide what
#     each cluster is: marker DotPlot, FeaturePlots, isotype, SHM, clonality,
#     cell cycle, and a DEG heatmap. Review these before filling in Part C.
#
#   PART B -- DEG finder
#     Run this to get differentially expressed genes for any cluster(s).
#     Change deg_mode, cluster_A, and cluster_B as needed and re-source.
#     Modes: "all_clusters", "one_vs_rest", "one_vs_one"
#     Note: cluster_B accepts a single value or a vector e.g. c(1, 2)
#
#   PART C -- Assign annotations and save
#     Fill in cluster_annotations with one label per cluster, then run.
#     Adds a cell_type column to the object and saves the annotated RDS.
#
# INPUT:
#   - integrated_S6.rds   produced by 06_Integrate_Samples.R
#
# OUTPUT:
#   - integrated_S7_annotated.rds   Seurat object with cell_type metadata
#   - DEG TSV files saved to deg_save_dir (Part B)
#
# NOTE — SMALL CLUSTERS:
#   Clusters with fewer than ~100 cells should be interpreted with caution.
#   Statistical power for both clustering and DEGs is limited at this size.
#
# NOTE — PRESTO:
#   FindAllMarkers() is slow by default. For a significant speed improvement:
#     install.packages("devtools")
#     devtools::install_github("immunogenomics/presto")
#   Seurat will use it automatically after installation.
#
# SEURAT VERSION NOTE:
#   Requires Seurat v5.
# =============================================================================


# ---- 0. PACKAGES ------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(viridis)
  library(scales)
  library(patchwork)
  library(pheatmap)
})


# =============================================================================
# !! ONLY MODIFY THIS SECTION !!
# =============================================================================

# --- Input / output paths -----------------------------------------------------
rds_in  <- "C:/Users/YourName/Documents/MyProject/RDS_Objects/integrated_S6.rds"
rds_out <- "C:/Users/YourName/Documents/MyProject/RDS_Objects/integrated_S7_annotated.rds"

# --- Part B: DEG settings -----------------------------------------------------
# Change deg_mode and cluster numbers, then re-source to run different comparisons.
# "all_clusters" : FindAllMarkers across all clusters (ignore cluster_A and cluster_B)
# "one_vs_rest"  : cluster_A vs all other cells
# "one_vs_one"   : cluster_A vs cluster_B. cluster_B can be a single cluster
#                  (e.g. cluster_B <- 1) or multiple clusters
#                  (e.g. cluster_B <- c(1, 2))
deg_mode     <- "all_clusters"
cluster_A    <- 0
cluster_B    <- 1
deg_min_pct  <- 0.25
deg_logfc    <- 0.25
deg_only_pos <- TRUE
deg_top_n    <- 20
deg_save_dir <- "C:/Users/YourName/Documents/MyProject/RDS_Objects"

# --- Part C: Cluster annotations ----------------------------------------------
# Fill in one label per cluster after reviewing the Part A plots and Part B DEGs.
# Use "Unknown" for any cluster you cannot confidently identify.
# If two clusters are the same cell type, give them the same label.

cluster_annotations <- c(
  "0"  = "MySample_CellType_0",
  "1"  = "MySample_CellType_1",
  "2"  = "MySample_CellType_2",
  "3"  = "MySample_CellType_3",
  "4"  = "MySample_CellType_4",
  "5"  = "MySample_CellType_5"
  # Add or remove entries to match the number of clusters in your object
)

# =============================================================================
# END OF MODIFIABLE SECTION — do not change anything below
# =============================================================================


# ---- Shared constants -------------------------------------------------------
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

s.genes   <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes


# =============================================================================
# LOAD
# =============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("Loading integrated object...")
if (!file.exists(rds_in)) {
  stop("Cannot find input RDS:\n  ", rds_in,
       "\nCheck rds_in and confirm 06_Integrate_Samples.R has been run.")
}

obj <- readRDS(rds_in)
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "seurat_clusters"

num_clusters <- table(obj$seurat_clusters)
message("Loaded: ", ncol(obj), " cells | ", length(num_clusters), " clusters")
message(paste(rep("=", 60), collapse = ""))


# ---- Cell cycle scoring (if not already present) ----------------------------
if (!"Phase" %in% colnames(obj@meta.data)) {
  message("Scoring cell cycle phases...")
  obj <- CellCycleScoring(obj, s.features = s.genes,
                          g2m.features = g2m.genes, set.ident = FALSE)
  message("Cell cycle scoring complete.")
} else {
  message("Phase column already present -- skipping scoring.")
}
obj$Phase <- factor(obj$Phase, levels = c("G1", "S", "G2M"))
Idents(obj) <- "seurat_clusters"


# =============================================================================
# PART A — ANNOTATION PLOTS
# =============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("PART A — Annotation plots")
message(paste(rep("=", 60), collapse = ""))


# --- A1. UMAP by cluster and by sample ---------------------------------------
message("Plot A1: UMAP by cluster and sample...")

print(
  DimPlot(obj, reduction = "umap", label = TRUE, label.size = 5, repel = TRUE) +
    scale_color_hue(labels = paste0(names(num_clusters), "  (n=", num_clusters, ")")) +
    ggtitle("Integrated UMAP — Clusters") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 9))
)

sample_col <- if ("sample_id" %in% colnames(obj@meta.data)) "sample_id" else "orig.ident"
print(
  DimPlot(obj, reduction = "umap", group.by = sample_col,
          label = FALSE, pt.size = 0.8) +
    ggtitle("Integrated UMAP — Sample") +
    theme(plot.title = element_text(hjust = 0.5))
)


# --- A2. Sample composition per cluster --------------------------------------
message("Plot A2: Sample composition per cluster...")

sample_df <- obj@meta.data %>%
  dplyr::rename(Sample = !!sample_col) %>%
  dplyr::select(Sample, seurat_clusters) %>%
  dplyr::group_by(seurat_clusters, Sample) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop_last") %>%
  dplyr::mutate(pct = 100 * n / sum(n)) %>%
  dplyr::ungroup()

print(
  ggplot(sample_df, aes(x = seurat_clusters, y = pct, fill = Sample)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(labels = scales::label_percent(scale = 1, accuracy = 1)) +
    ylab("% of cluster") + xlab("Cluster") +
    theme_classic(base_size = 12) +
    ggtitle("Sample composition per cluster")
)


# --- A3. B cell marker DotPlot -----------------------------------------------
message("Plot A3: B cell marker DotPlot...")

markers_present <- intersect(b_cell_markers, rownames(obj))
print(
  DotPlot(obj, features = markers_present,
          cols = c("gray90", "red"), dot.scale = 6) +
    RotatedAxis() +
    ggtitle("B cell marker panel — integrated clusters") +
    theme(axis.text.x = element_text(size = 8))
)


# --- A4. FeaturePlots --------------------------------------------------------
message("Plot A4: FeaturePlots of key markers...")

markers_for_feature <- intersect(
  c("CD79A", "IGHD", "IGHM", "CD27", "BCL6", "AICDA",
    "MKI67", "PRDM1", "XBP1", "FCRL4", "FCRL5"),
  rownames(obj)
)
if (length(markers_for_feature) > 0) {
  print(
    FeaturePlot(obj, reduction = "umap", features = markers_for_feature,
                order = TRUE, pt.size = 0.8) &
      scale_color_gradientn(colors = viridis::inferno(256))
  )
}


# --- A5. Isotype per cluster -------------------------------------------------
message("Plot A5: Isotype per cluster...")

if ("c_call" %in% colnames(obj@meta.data)) {

  obj_bcr <- subset(obj, subset = !is.na(c_call) & c_call != "")
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
      ggtitle("Isotype composition per cluster")
  )

} else {
  message("  No c_call column found -- skipping isotype plots.")
}


# --- A6. SHM frequency and clone size ----------------------------------------
message("Plot A6: SHM frequency and clone size...")

if ("mu_freq_H" %in% colnames(obj@meta.data)) {
  print(
    FeaturePlot(obj, features = "mu_freq_H", reduction = "umap",
                order = TRUE, pt.size = 0.8) +
      scale_color_gradientn(colors = c("#FFFFCC", "#41B6C4", "#0C2C84"),
                            na.value = "grey90") +
      ggtitle("Integrated UMAP — SHM frequency, heavy chain")
  )

  shm_df <- obj@meta.data %>%
    dplyr::filter(!is.na(mu_freq_H)) %>%
    dplyr::mutate(seurat_clusters = as.factor(seurat_clusters))

  print(
    ggplot(shm_df, aes(x = seurat_clusters, y = mu_freq_H)) +
      geom_boxplot(outlier.size = 0.5, fill = "steelblue", alpha = 0.6) +
      ylab("SHM frequency (heavy chain)") + xlab("Cluster") +
      theme_classic(base_size = 12) +
      ggtitle("SHM frequency per cluster")
  )
}

if ("clone_count" %in% colnames(obj@meta.data)) {
  print(
    FeaturePlot(obj, features = "clone_count", reduction = "umap",
                order = TRUE, pt.size = 0.8) +
      scale_color_gradientn(colors = c("#FFFFCC", "#FD8D3C", "#800026"),
                            na.value = "grey90") +
      ggtitle("Integrated UMAP — Clone size")
  )
}


# --- A7. Clonality per cluster -----------------------------------------------
message("Plot A7: Clonality per cluster...")

if ("clone_id" %in% colnames(obj@meta.data)) {

  clone_df <- obj@meta.data %>%
    dplyr::filter(!is.na(clone_id), !is.na(seurat_clusters)) %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::summarise(
      n_cells         = dplyr::n(),
      n_clones        = dplyr::n_distinct(clone_id),
      clonality       = 1 - (n_clones / n_cells),
      mean_clone_size = n_cells / n_clones,
      max_clone_size  = max(table(clone_id)),
      pct_expanded    = 100 * mean(table(clone_id)[clone_id] > 1),
      .groups = "drop"
    ) %>%
    dplyr::arrange(as.numeric(as.character(seurat_clusters)))

  message("  Clonality summary:")
  print(as.data.frame(clone_df))

  print(
    ggplot(clone_df, aes(x = as.factor(seurat_clusters), y = clonality)) +
      geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
      ylab("Clonality (1 - clones/cells)") + xlab("Cluster") +
      theme_classic(base_size = 12) +
      ggtitle("Clonality per cluster\n(higher = more clonal expansion)")
  )

  print(
    ggplot(clone_df, aes(x = as.factor(seurat_clusters), y = mean_clone_size)) +
      geom_bar(stat = "identity", fill = "darkorange", alpha = 0.8) +
      ylab("Mean clone size") + xlab("Cluster") +
      theme_classic(base_size = 12) +
      ggtitle("Mean clone size per cluster")
  )

} else {
  message("  No clone_id column found -- skipping clonality plots.")
}


# --- A8. Cell cycle per cluster ----------------------------------------------
message("Plot A8: Cell cycle per cluster...")

Idents(obj) <- "Phase"
print(
  DimPlot(obj, reduction = "umap", label = FALSE, pt.size = 0.8,
          cols = c("G1" = "#AAAAAA", "S" = "#E69F00", "G2M" = "#CC79A7")) +
    ggtitle("Integrated UMAP — Cell cycle phase")
)
Idents(obj) <- "seurat_clusters"

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
    ggtitle("Cell cycle composition per cluster")
)


# --- A9. DEG heatmap: top 5 markers per cluster ------------------------------
message("Plot A9: DEG heatmap (FindAllMarkers -- may take a few minutes)...")

all_markers <- FindAllMarkers(
  obj,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  verbose         = FALSE
)

top5 <- all_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 5) %>%
  dplyr::ungroup()

message("  FindAllMarkers complete. Total DEGs: ", nrow(all_markers))

top20 <- all_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 20) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(as.numeric(as.character(cluster)), desc(avg_log2FC))

message("  Top 20 markers per cluster:")
print(as.data.frame(top20[, c("cluster", "gene", "avg_log2FC",
                               "pct.1", "pct.2", "p_val_adj")]),
      row.names = FALSE)

if (nrow(top5) > 0) {
  genes_to_plot <- unique(as.character(top5$gene))
  genes_to_plot <- intersect(genes_to_plot, rownames(obj))

  avg_expr <- AverageExpression(
    obj,
    features = genes_to_plot,
    group.by = "seurat_clusters",
    assays   = "RNA",
    layer    = "data",
    verbose  = FALSE
  )$RNA

  # Preserve row and column names explicitly before z-scoring.
  # t(scale(t())) can drop dimnames in some Seurat/R version combinations.
  # AverageExpression also prepends 'g' to column names when group.by starts
  # with a number -- we store and restore names to handle both issues cleanly.
  col_names  <- colnames(avg_expr)
  row_names  <- rownames(avg_expr)

  avg_scaled <- t(scale(t(avg_expr)))
  avg_scaled[!is.finite(avg_scaled)] <- 0

  rownames(avg_scaled) <- row_names
  colnames(avg_scaled) <- col_names

  gene_order <- unique(
    as.character(top5$gene)[as.character(top5$gene) %in% rownames(avg_scaled)]
  )
  cluster_order <- colnames(avg_scaled)[
    order(as.numeric(gsub("^g", "", colnames(avg_scaled))))
  ]

  avg_scaled <- avg_scaled[gene_order, cluster_order, drop = FALSE]

  pheatmap::pheatmap(
    avg_scaled,
    cluster_rows  = FALSE,
    cluster_cols  = FALSE,
    color         = colorRampPalette(c("dodgerblue", "white", "firebrick"))(100),
    fontsize_row  = 7,
    fontsize_col  = 9,
    angle_col     = 0,
    main          = "Top 5 DEGs per cluster (z-scored average expression)",
    border_color  = NA
  )
} else {
  message("  No DEGs found above thresholds.")
}

message("\nPart A complete. Review all plots, then run Part B and/or Part C.")


# =============================================================================
# PART B — DEG FINDER
# =============================================================================
# Change deg_mode / cluster_A / cluster_B in the USER SETTINGS above,
# then re-source from here to run different comparisons.
# =============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("PART B — DEG finder (mode: ", deg_mode, ")")
message(paste(rep("=", 60), collapse = ""))

if (!dir.exists(deg_save_dir)) dir.create(deg_save_dir, recursive = TRUE)
Idents(obj) <- "seurat_clusters"

if (deg_mode == "all_clusters") {

  message("Running FindAllMarkers across all clusters...")
  deg_result <- FindAllMarkers(
    obj,
    min.pct         = deg_min_pct,
    logfc.threshold = deg_logfc,
    only.pos        = deg_only_pos,
    verbose         = FALSE
  )
  deg_result <- deg_result %>% dplyr::arrange(cluster, desc(avg_log2FC))

  top_per_cluster <- deg_result %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(order_by = avg_log2FC, n = deg_top_n) %>%
    dplyr::ungroup()

  message("Top ", deg_top_n, " markers per cluster:")
  print(as.data.frame(
    top_per_cluster[, c("cluster", "gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
  ))

  out_file <- file.path(deg_save_dir, "DEGs_all_clusters.tsv")
  write.table(deg_result, file = out_file, sep = "\t",
              quote = FALSE, row.names = FALSE)
  message("Full DEG table saved: ", out_file)

} else if (deg_mode == "one_vs_rest") {

  message("Running FindMarkers: cluster ", cluster_A, " vs all others...")
  deg_result <- FindMarkers(
    obj,
    ident.1         = cluster_A,
    min.pct         = deg_min_pct,
    logfc.threshold = deg_logfc,
    only.pos        = deg_only_pos,
    verbose         = FALSE
  )
  deg_result$gene    <- rownames(deg_result)
  deg_result$cluster <- cluster_A
  deg_result         <- deg_result %>% dplyr::arrange(desc(avg_log2FC))

  message("Top ", deg_top_n, " markers for cluster ", cluster_A, ":")
  print(head(
    deg_result[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")],
    deg_top_n
  ))

  out_file <- file.path(deg_save_dir,
                        paste0("DEGs_cluster", cluster_A, "_vs_rest.tsv"))
  write.table(deg_result, file = out_file, sep = "\t",
              quote = FALSE, row.names = FALSE)
  message("Full DEG table saved: ", out_file)

} else if (deg_mode == "one_vs_one") {

  message("Running FindMarkers: cluster ", cluster_A,
          " vs cluster(s) ", paste(cluster_B, collapse = ", "), "...")
  deg_result <- FindMarkers(
    obj,
    ident.1         = cluster_A,
    ident.2         = cluster_B,
    min.pct         = deg_min_pct,
    logfc.threshold = deg_logfc,
    only.pos        = deg_only_pos,
    verbose         = FALSE
  )
  deg_result$gene       <- rownames(deg_result)
  deg_result$cluster    <- cluster_A
  deg_result$comparison <- paste0("cluster", cluster_A,
                                  "_vs_cluster", paste(cluster_B, collapse = "_"))
  deg_result            <- deg_result %>% dplyr::arrange(desc(avg_log2FC))

  message("Top ", deg_top_n, " markers for cluster ", cluster_A,
          " vs cluster(s) ", paste(cluster_B, collapse = ", "), ":")
  print(head(
    deg_result[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")],
    deg_top_n
  ))

  out_file <- file.path(deg_save_dir,
                        paste0("DEGs_cluster", cluster_A,
                               "_vs_cluster", paste(cluster_B, collapse = "_"), ".tsv"))
  write.table(deg_result, file = out_file, sep = "\t",
              quote = FALSE, row.names = FALSE)
  message("Full DEG table saved: ", out_file)

} else {
  stop("Unknown deg_mode: '", deg_mode,
       "'. Choose from: all_clusters, one_vs_rest, one_vs_one")
}


# =============================================================================
# PART C — ASSIGN ANNOTATIONS AND SAVE
# =============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("PART C — Assigning annotations and saving")
message(paste(rep("=", 60), collapse = ""))

# Validate coverage
clusters_in_obj  <- as.character(sort(unique(obj$seurat_clusters)))
clusters_defined <- names(cluster_annotations)
missing_clusters <- setdiff(clusters_in_obj, clusters_defined)
extra_clusters   <- setdiff(clusters_defined, clusters_in_obj)

if (length(missing_clusters) > 0) {
  warning(
    "The following clusters have no annotation defined: ",
    paste(missing_clusters, collapse = ", "),
    "\nThese cells will be labeled 'Unlabeled'. Add them to cluster_annotations."
  )
}
if (length(extra_clusters) > 0) {
  message("Note: annotation includes cluster(s) not in object: ",
          paste(extra_clusters, collapse = ", "), " -- ignored.")
}

# Apply annotations
# unname() is required -- Seurat matches metadata by cell barcode, not by the
# cluster names that would otherwise be carried as vector names
full_annotation_map               <- cluster_annotations
full_annotation_map[missing_clusters] <- "Unlabeled"
obj$cell_type <- unname(full_annotation_map[as.character(obj$seurat_clusters)])

# Summary table
annotation_summary <- obj@meta.data %>%
  dplyr::group_by(seurat_clusters, cell_type) %>%
  dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") %>%
  dplyr::arrange(as.numeric(as.character(seurat_clusters)))

message("Annotation summary:")
print(as.data.frame(annotation_summary))

# Diagnostic UMAP: cell type label with cluster number and n=
annotation_summary$label <- paste0(
  annotation_summary$cell_type,
  "\n(C", annotation_summary$seurat_clusters,
  ", n=", annotation_summary$n_cells, ")"
)
label_map           <- setNames(annotation_summary$label,
                                as.character(annotation_summary$seurat_clusters))
obj$cell_type_label <- unname(label_map[as.character(obj$seurat_clusters)])

Idents(obj) <- "cell_type_label"
print(
  DimPlot(obj, reduction = "umap", label = TRUE, label.size = 3.5,
          repel = TRUE) +
    ggtitle("Integrated UMAP — Annotated (diagnostic)") +
    theme(plot.title  = element_text(hjust = 0.5),
          legend.text = element_text(size = 8)) +
    NoLegend()
)

Idents(obj) <- "cell_type"
print(
  DimPlot(obj, reduction = "umap", label = TRUE, label.size = 4,
          repel = TRUE) +
    ggtitle("Integrated UMAP — Cell types") +
    theme(plot.title = element_text(hjust = 0.5)) +
    guides(color = guide_legend(override.aes = list(size = 4)))
)

Idents(obj) <- "seurat_clusters"

# Save
saveRDS(obj, file = rds_out)

message("\n", paste(rep("=", 60), collapse = ""))
message("Done.")
message("Saved: ", rds_out)
message("Cell types: ", paste(sort(unique(obj$cell_type)), collapse = ", "))
message("Cells: ", ncol(obj))
message(paste(rep("=", 60), collapse = ""))
