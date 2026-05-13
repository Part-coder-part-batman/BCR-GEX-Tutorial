# =============================================================================
# 09_Shared_Clone_Analysis.R
# Cross-sample clonal analysis with scRepertoire
# =============================================================================
#
# HOW TO USE:
#   1. Fill in the file paths and settings in the USER SETTINGS section below
#   2. Run the script (Ctrl+Shift+Enter in RStudio to source the whole file)
#
# This script is split into four parts that you can run sequentially or
# independently (Parts B-D require Part A to have been run first):
#
#   PART A -- Prepare AIRR data and run combineBCR()
#     Loads the annotated BCR TSVs from step 08, maps Immcantation column
#     names to the fields scRepertoire expects, and runs combineBCR() to
#     assign cross-sample CTstrict clonotypes. Builds a master BCR table
#     with CTstrict joined onto every cell. Saves combined_BCR_S9.rds and
#     All_S9_bcr_with_CTstrict.tsv.
#
#   PART B -- Repertoire characterization and shared clone quantification
#     Clonal homeostasis plot (clone size distribution per sample) and a
#     stacked bar chart of cells in shared vs. tissue-exclusive clones.
#     Prints shared clone counts to the console and saves a summary table.
#
#   PART C -- Attach clonal data to the Seurat object
#     Runs combineExpression() to add CTstrict, Frequency, and cloneType
#     to every cell in the integrated Seurat object. Saves the updated
#     Seurat object and generates clonal network and clonal occupancy plots.
#
#   PART D -- Shared clone dot plot
#     For each patient, shows which cell type clusters each shared LN/PT
#     clone occupies. Dots are colored by dominant isotype, shaped by tissue
#     of origin (LN / PT / both), and sized by cell count.
#
# INPUT FILES:
#   - <SAMPLE_ID>_S8_bcr_annotated.tsv    from step 08, one per sample
#   - integrated_S7_annotated.rds          from step 07
#
# OUTPUT:
#   - combined_BCR_S9.rds
#   - All_S9_bcr_with_CTstrict.tsv
#   - integrated_S9_scRep.rds
#   - Plots/Part9/clonal_homeostasis.pdf
#   - Plots/Part9/shared_cells_bar_<PATIENT_ID>.pdf
#   - Plots/Part9/clonal_network_all.pdf
#   - Plots/Part9/clonal_network_<CLUSTER>.pdf
#   - Plots/Part9/clonal_occupy_counts.pdf
#   - Plots/Part9/clonal_occupy_proportion.pdf
#   - Plots/Part9/shared_clones_dotplot_<PATIENT_ID>.pdf
#
# NOTE — clonalOccupy / occupiedscRepertoire:
#   clonalOccupy() was renamed to occupiedscRepertoire() in scRepertoire v2.
#   This script calls occupiedscRepertoire(). If you are on an older version
#   replace it with clonalOccupy() and the same arguments will work.
# =============================================================================


# =============================================================================
# !! ONLY EDIT THIS SECTION !!
# =============================================================================

# Annotated BCR TSV files from step 08 — one entry per sample
# Format: list(<sample_id> = <path_to_tsv>, ...)
BCR_TSV <- list(
  P1_LN = "C:/Users/YourName/Documents/MyProject/BCR_Data/Annotated/P1_LN_S8_bcr_annotated.tsv",
  P1_PT = "C:/Users/YourName/Documents/MyProject/BCR_Data/Annotated/P1_PT_S8_bcr_annotated.tsv"
)

# Annotated Seurat object from step 07
RDS_IN <- "C:/Users/YourName/Documents/MyProject/RDS_Objects/integrated_S7_annotated.rds"

# Output directory for BCR tables
OUT_DIR <- "C:/Users/YourName/Documents/MyProject/BCR_Data/Annotated"

# Output directory for RDS objects
RDS_DIR <- "C:/Users/YourName/Documents/MyProject/RDS_Objects"

# Output directory for plots
PLOT_DIR <- "C:/Users/YourName/Documents/MyProject/Plots/Part9"

# Patient label used in plot titles and output filenames (e.g. "P1")
PATIENT_ID <- "P1"

# Cluster of interest for the filtered clonal network plot (Part C)
# Replace with whichever cell type you want to examine
CLUSTER_OF_INTEREST <- "PC"

# Cell type cluster order for the dot plot x-axis.
# Replace with the cell type labels from your own annotation, in the order
# you want them to appear left-to-right on the plot.
CLUSTER_ORDER <- c(
  "CellType_1",
  "CellType_2",
  "CellType_3",
  "CellType_4",
  "CellType_5",
  "CellType_6",
  "CellType_7",
  "CellType_8",
  "CellType_9",
  "CellType_10",
  "Unsure"
  # Add or remove entries to match the clusters in your dataset
)

# Isotype color palette — consistent with the donut palette in step 08
ISOTYPE_COLORS <- c(
  IGHD  = "#E41A1C",
  IGHM  = "#377EB8",
  IGHG1 = "#4DAF4A",
  IGHG2 = "#984EA3",
  IGHG3 = "#FF7F00",
  IGHG4 = "#FFFF33",
  IGHA1 = "#A65628",
  IGHA2 = "#F781BF"
)

# =============================================================================
# END OF EDITABLE SECTION — do not change anything below
# =============================================================================


# ---- Load packages ----------------------------------------------------------
suppressPackageStartupMessages({
  library(scRepertoire)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(data.table)
  library(stringr)
  library(ggraph)     # required by clonalNetwork — load before plotting
  library(patchwork)
})


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
dir.create(RDS_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)


# ---- Helper: map Immcantation column names to scRepertoire expectations -----
# combineBCR() expects specific column names that differ from the AIRR/
# Immcantation convention. This function adds the expected aliases without
# removing the originals.
prepare_for_scRepertoire <- function(df) {

  if (!"barcode"  %in% names(df) && "cell_id"    %in% names(df)) df$barcode  <- df$cell_id
  if (!"sample"   %in% names(df) && "sample_id"  %in% names(df)) df$sample   <- df$sample_id
  if (!"chain"    %in% names(df) && "locus"      %in% names(df)) df$chain    <- df$locus
  if (!"v_gene"   %in% names(df) && "v_call"     %in% names(df)) df$v_gene   <- df$v_call
  if (!"d_gene"   %in% names(df) && "d_call"     %in% names(df)) df$d_gene   <- df$d_call
  if (!"j_gene"   %in% names(df) && "j_call"     %in% names(df)) df$j_gene   <- df$j_call
  if (!"c_gene"   %in% names(df) && "c_call"     %in% names(df)) df$c_gene   <- df$c_call

  # CDR3 sequences — scRepertoire uses cdr3_nt for the CTstrict distance calc
  if (!"cdr3"    %in% names(df) && "junction_aa" %in% names(df)) df$cdr3    <- df$junction_aa
  if (!"cdr3_nt" %in% names(df) && "junction"    %in% names(df)) df$cdr3_nt <- df$junction
  if (!"cdr3_aa" %in% names(df) && "junction_aa" %in% names(df)) df$cdr3_aa <- df$junction_aa

  # UMI count as reads proxy (avoids NAs in scRepertoire internals)
  if (!"reads"   %in% names(df) && "umi_count"   %in% names(df)) df$reads   <- df$umi_count
  if (!"reads"   %in% names(df))                                  df$reads   <- 1L

  return(df)
}


# ---- Helper: pick dominant isotype for a clone x cluster --------------------
# Ties broken by isotype priority (IgG first, then IgA, IgD, IgM).
ISOTYPE_PRIORITY <- c("IGHG1","IGHG2","IGHG3","IGHG4","IGHA1","IGHA2","IGHD","IGHM")

pick_dominant_isotype <- function(x) {
  x <- x[!is.na(x) & x != "IGHE"]
  if (length(x) == 0) return(NA_character_)
  tab  <- table(x)
  mx   <- max(tab)
  tied <- names(tab)[tab == mx]
  if (length(tied) == 1) return(tied)
  tied <- intersect(ISOTYPE_PRIORITY, tied)
  if (length(tied) > 0) return(tied[1])
  sort(names(tab)[tab == mx])[1]
}


# =============================================================================
# PART A — PREPARE AIRR DATA AND RUN combineBCR()
# =============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("PART A — Prepare AIRR data and run combineBCR()")
message(paste(rep("=", 60), collapse = ""))


# ---- A1. Load and map BCR TSVs ----------------------------------------------
message("\n>> A1: Loading and mapping BCR TSVs...")

airr_list <- lapply(BCR_TSV, function(path) {
  data.table::fread(path, sep = "\t", data.table = FALSE)
})
airr_list <- lapply(airr_list, prepare_for_scRepertoire)

# Confirm required columns are present after mapping
required_cols <- c("barcode", "chain", "v_gene", "j_gene", "cdr3_nt")
for (sid in names(airr_list)) {
  missing <- setdiff(required_cols, colnames(airr_list[[sid]]))
  if (length(missing) > 0)
    stop(sid, ": missing required columns after mapping: ",
         paste(missing, collapse = ", "))
}
message("Column mapping OK for all samples.")


# ---- A2. Run combineBCR() ---------------------------------------------------
message("\n>> A2: Running combineBCR() across ", length(airr_list), " samples...")

# Note on barcode format: combineBCR() prepends the sample name to each
# barcode, producing e.g. "P1_LN_ACTGCTCA...-1_P1_LN". This double
# suffix/prefix is expected — our cell_ids already carry the sample suffix
# from step 01, and combineBCR() adds its own prefix on top. The ct_lookup
# step in A3 strips the leading prefix to recover the original cell_id format.

combined_BCR <- combineBCR(
  airr_list,
  samples     = names(airr_list),
  removeNA    = FALSE,
  removeMulti = TRUE
)

message("combineBCR() complete.")
message("Cells per sample:")
print(sapply(combined_BCR, nrow))
message("Example CTstrict values (", names(combined_BCR)[1], "):")
print(head(combined_BCR[[1]]$CTstrict, 6))

saveRDS(combined_BCR, file.path(RDS_DIR, "combined_BCR_S9.rds"))
message("Saved: combined_BCR_S9.rds")


# ---- A3. Build CTstrict lookup ----------------------------------------------
message("\n>> A3: Building CTstrict lookup table...")

# combineBCR() prepended "<sample>_" to each barcode. Strip it so the
# cell_id matches the Part 8 format (ACGT...-1_P1_LN).
ct_lookup <- do.call(rbind, lapply(names(combined_BCR), function(sid) {
  df <- combined_BCR[[sid]]
  data.frame(
    sample_id = df$sample,
    cell_id   = sub(paste0("^", sid, "_"), "", df$barcode),
    CTstrict  = df$CTstrict,
    CTaa      = df$CTaa,
    CTnt      = df$CTnt,
    CTgene    = df$CTgene,
    stringsAsFactors = FALSE
  )
}))

message("CTstrict lookup: ", nrow(ct_lookup), " rows")
message("Example cell_id values after stripping:")
print(head(ct_lookup$cell_id, 6))


# ---- A4. Build master BCR table and join CTstrict ---------------------------
message("\n>> A4: Building master BCR table and joining CTstrict...")

bcr_all <- do.call(rbind, lapply(names(BCR_TSV), function(sid) {
  data.table::fread(BCR_TSV[[sid]], sep = "\t", data.table = FALSE)
}))

message("Master BCR table: ", nrow(bcr_all), " rows | ",
        length(unique(bcr_all$cell_id[bcr_all$locus == "IGH"])),
        " cells with heavy chain")

bcr_all <- dplyr::left_join(bcr_all, ct_lookup, by = c("sample_id", "cell_id"))

# Check match rate on heavy chains
n_igh     <- sum(bcr_all$locus == "IGH")
n_matched <- sum(!is.na(bcr_all$CTstrict[bcr_all$locus == "IGH"]))
message("CTstrict matched: ", n_matched, " / ", n_igh,
        " heavy chain cells (", round(100 * n_matched / n_igh, 1), "%)")

if (n_matched / n_igh < 0.5) {
  warning("Less than 50% of heavy chain cells matched a CTstrict. ",
          "Check barcode format consistency between step 01 and step 08.")
}

# Tissue label derived from sample_id
bcr_all <- bcr_all %>%
  dplyr::mutate(
    tissue = dplyr::case_when(
      grepl("_LN$|^LN", sample_id) ~ "LN",
      grepl("_PT$|^PT", sample_id) ~ "PT",
      TRUE ~ "UNK"
    )
  )

data.table::fwrite(bcr_all,
                   file.path(OUT_DIR, "All_S9_bcr_with_CTstrict.tsv"),
                   sep = "\t")
message("Saved: All_S9_bcr_with_CTstrict.tsv")

message("\nPart A complete.")


# =============================================================================
# PART B — REPERTOIRE CHARACTERIZATION AND SHARED CLONE QUANTIFICATION
# =============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("PART B — Repertoire characterization and shared clone quantification")
message(paste(rep("=", 60), collapse = ""))


# ---- B1. Clonal homeostasis -------------------------------------------------
message("\n>> B1: Clonal homeostasis plot...")

p_homeostasis <- clonalHomeostasis(
  combined_BCR,
  cloneCall = "strict"
)
print(p_homeostasis)

ggsave(file.path(PLOT_DIR, "clonal_homeostasis.pdf"),
       p_homeostasis, width = 6, height = 5)
message("Saved: clonal_homeostasis.pdf")


# ---- B2. Count shared clones and cells --------------------------------------
message("\n>> B2: Counting shared clones and cells...")

bcr_igh <- bcr_all %>%
  dplyr::filter(locus == "IGH", !is.na(CTstrict), CTstrict != "")

clones_ln  <- unique(bcr_igh$CTstrict[bcr_igh$tissue == "LN"])
clones_pt  <- unique(bcr_igh$CTstrict[bcr_igh$tissue == "PT"])
shared_cts <- intersect(clones_ln, clones_pt)

message("CTstrict clones in LN only:   ", length(setdiff(clones_ln, clones_pt)))
message("CTstrict clones in PT only:   ", length(setdiff(clones_pt, clones_ln)))
message("CTstrict clones shared LN+PT: ", length(shared_cts))

cells_ln_shared <- sum(bcr_igh$CTstrict[bcr_igh$tissue == "LN"] %in% shared_cts)
cells_pt_shared <- sum(bcr_igh$CTstrict[bcr_igh$tissue == "PT"] %in% shared_cts)

# Raw counts only — no percentages (denominator choice is non-trivial when
# comparing a pure CTstrict approach against a hybrid Immcantation method)
message("LN cells in shared clones:    ", cells_ln_shared)
message("PT cells in shared clones:    ", cells_pt_shared)
message("Total shared cells:           ", cells_ln_shared + cells_pt_shared)

# Summary table
shared_clone_tbl <- bcr_igh %>%
  dplyr::filter(CTstrict %in% shared_cts) %>%
  dplyr::count(CTstrict, tissue, name = "n_cells") %>%
  tidyr::pivot_wider(
    names_from  = tissue,
    values_from = n_cells,
    values_fill = 0L
  ) %>%
  dplyr::mutate(total_cells = LN + PT) %>%
  dplyr::arrange(dplyr::desc(total_cells))

message("\nShared clone table (top 10):")
print(head(shared_clone_tbl, 10))


# ---- B3. Shared vs. exclusive cell count bar chart -------------------------
message("\n>> B3: Shared vs. exclusive cell count bar chart...")

shared_bar_df <- bcr_igh %>%
  dplyr::mutate(
    status = dplyr::if_else(CTstrict %in% shared_cts, "Shared", "Exclusive")
  ) %>%
  dplyr::count(tissue, status, name = "n_cells") %>%
  dplyr::mutate(
    tissue = factor(tissue, levels = c("LN", "PT")),
    status = factor(status, levels = c("Shared", "Exclusive"))
  )

p_shared_bar <- ggplot(shared_bar_df,
                       aes(x = tissue, y = n_cells, fill = status)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = n_cells),
            position = position_stack(vjust = 0.5),
            color = "white", size = 4, fontface = "bold") +
  scale_fill_manual(
    values = c(Shared = "#2166AC", Exclusive = "#B2B2B2"),
    name   = "Clone status"
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.x = element_blank()) +
  labs(
    x     = "Tissue",
    y     = "Number of cells (IGH)",
    title = paste0("Cells in shared vs. exclusive clones — ", PATIENT_ID)
  )
print(p_shared_bar)

ggsave(file.path(PLOT_DIR, paste0("shared_cells_bar_", PATIENT_ID, ".pdf")),
       p_shared_bar, width = 5, height = 4)
message("Saved: shared_cells_bar_", PATIENT_ID, ".pdf")

message("\nPart B complete.")


# =============================================================================
# PART C — ATTACH CLONAL DATA TO SEURAT OBJECT
# =============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("PART C — Attach clonal data to Seurat object")
message(paste(rep("=", 60), collapse = ""))


# ---- C1. Load Seurat object -------------------------------------------------
message("\n>> C1: Loading Seurat object...")

seu <- readRDS(RDS_IN)
message("Loaded: ", ncol(seu), " cells")


# ---- C2. Harmonize barcodes and run combineExpression() ---------------------
message("\n>> C2: Harmonizing barcodes and running combineExpression()...")

# Strip the leading "<sample>_" prefix that combineBCR() added so barcodes
# match the Seurat cell name format (ACGT...-1_P1_LN)
combined_BCR_harmonized <- lapply(combined_BCR, function(df) {
  df$barcode <- sub(paste0("^", df$sample, "_"), "", df$barcode)
  df
})

# Overlap check before proceeding
airr_barcodes <- unique(unlist(lapply(combined_BCR_harmonized, function(x) x$barcode)))
n_overlap <- sum(colnames(seu) %in% airr_barcodes)
message("Barcode overlap (Seurat vs AIRR): ", n_overlap, " / ", ncol(seu))

if (n_overlap == 0) {
  stop(
    "Zero overlap after barcode harmonization.\n",
    "Check that cell_id format in step 08 TSVs matches Seurat colnames.\n",
    "Example Seurat names:  ", paste(head(colnames(seu), 3), collapse = ", "), "\n",
    "Example AIRR barcodes: ", paste(head(airr_barcodes, 3), collapse = ", ")
  )
}

# Match Seurat cell names to sample IDs for combineExpression grouping.
# Cell names look like "ACGT...-1_P1_LN"; we match against the known sample IDs.
sample_ids  <- names(combined_BCR_harmonized)
seu$sample  <- NA_character_
for (sid in sample_ids) {
  seu$sample[grepl(paste0("_", sid, "$"), colnames(seu))] <- sid
}
message("Seurat sample labels:")
print(table(seu$sample, useNA = "ifany"))

# proportion = TRUE is required so Frequency is stored as a proportion.
# The cloneType bins (Rare/Small/Medium/Large/Hyperexpanded) use proportion
# thresholds (0.0001 / 0.001 / 0.01 / 0.1 / 1). With proportion = FALSE,
# raw counts are stored and every clone exceeds the 0.1 threshold, collapsing
# all cells into "Hyperexpanded".
# group.by is intentionally omitted: in scRepertoire v2, combining group.by
# with proportion = TRUE causes frequencies to be calculated incorrectly,
# also resulting in all cells landing in "Hyperexpanded".
seu_scRep <- combineExpression(
  combined_BCR_harmonized,
  sc         = seu,
  cloneCall  = "strict",
  proportion = TRUE,
  filterNA   = FALSE
)

message("combineExpression() complete.")
message("New metadata columns added:")
print(setdiff(colnames(seu_scRep@meta.data), colnames(seu@meta.data)))

if (!"CTstrict" %in% colnames(seu_scRep@meta.data)) {
  warning("CTstrict column not found in Seurat metadata. ",
          "Check that cloneCall = 'strict' matches the combined_BCR columns.")
} else {
  message("Cells with CTstrict assigned: ",
          sum(!is.na(seu_scRep$CTstrict)), " / ", ncol(seu_scRep))
}

saveRDS(seu_scRep, file.path(RDS_DIR, "integrated_S9_scRep.rds"))
message("Saved: integrated_S9_scRep.rds")


# ---- C3. Clonal network on UMAP ---------------------------------------------
message("\n>> C3: Clonal network plots...")

# Full network — all clusters
p_net_all <- clonalNetwork(
  seu_scRep,
  reduction       = "umap",
  identity        = "cell_type",
  filter.clones   = NULL,
  filter.identity = NULL,
  cloneCall       = "strict"
)
print(p_net_all)

ggsave(file.path(PLOT_DIR, "clonal_network_all.pdf"),
       p_net_all, width = 9, height = 7)
message("Saved: clonal_network_all.pdf")

# Filtered network — single cluster of interest
p_net_filtered <- clonalNetwork(
  seu_scRep,
  reduction       = "umap",
  identity        = "cell_type",
  filter.identity = CLUSTER_OF_INTEREST,
  filter.clones   = NULL,
  cloneCall       = "strict"
)
print(p_net_filtered)

ggsave(file.path(PLOT_DIR,
                 paste0("clonal_network_", CLUSTER_OF_INTEREST, ".pdf")),
       p_net_filtered, width = 9, height = 7)
message("Saved: clonal_network_", CLUSTER_OF_INTEREST, ".pdf")


# ---- C4. Clonal occupancy by cluster ----------------------------------------
message("\n>> C4: Clonal occupancy plots...")

# Note: clonalOccupy() was renamed to occupiedscRepertoire() in scRepertoire
# v2.0. If you are on an older version, replace occupiedscRepertoire() with
# clonalOccupy() — the arguments are the same.

p_occupy <- occupiedscRepertoire(
  seu_scRep,
  x.axis     = "cell_type",
  proportion = FALSE,
  label      = TRUE
)
print(p_occupy + theme(axis.text.x = element_text(angle = 45, hjust = 1)))

ggsave(file.path(PLOT_DIR, "clonal_occupy_counts.pdf"),
       p_occupy + theme(axis.text.x = element_text(angle = 45, hjust = 1)),
       width = 9, height = 5)

p_occupy_prop <- occupiedscRepertoire(
  seu_scRep,
  x.axis     = "cell_type",
  proportion = TRUE,
  label      = FALSE
)
print(p_occupy_prop + theme(axis.text.x = element_text(angle = 45, hjust = 1)))

ggsave(file.path(PLOT_DIR, "clonal_occupy_proportion.pdf"),
       p_occupy_prop + theme(axis.text.x = element_text(angle = 45, hjust = 1)),
       width = 9, height = 5)

message("Saved: clonal_occupy_counts.pdf and clonal_occupy_proportion.pdf")

message("\nPart C complete.")


# =============================================================================
# PART D — SHARED CLONE DOT PLOT
# =============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("PART D — Shared clone dot plot")
message(paste(rep("=", 60), collapse = ""))


# ---- D1. Build plot data ----------------------------------------------------
message("\n>> D1: Building plot data...")

# Heavy chain rows with valid CTstrict, cell_type, and isotype
bcr_plot <- bcr_all %>%
  dplyr::filter(
    locus     == "IGH",
    !is.na(CTstrict),  CTstrict  != "",
    !is.na(cell_type), cell_type != "",
    !is.na(c_call),    c_call    != "",
    c_call    != "IGHE"
  ) %>%
  dplyr::filter(!is.na(tissue))

# Identify clones present in both LN and PT
shared_cts_plot <- bcr_plot %>%
  dplyr::group_by(CTstrict) %>%
  dplyr::summarise(
    has_LN = any(tissue == "LN"),
    has_PT = any(tissue == "PT"),
    .groups = "drop"
  ) %>%
  dplyr::filter(has_LN & has_PT) %>%
  dplyr::pull(CTstrict)

message("Shared CTstrict clones with cell type annotation: ",
        length(shared_cts_plot))

# Collapse to one row per (CTstrict x cell_type)
plot_df <- bcr_plot %>%
  dplyr::filter(CTstrict %in% shared_cts_plot) %>%
  dplyr::group_by(CTstrict, cell_type) %>%
  dplyr::summarise(
    Count       = dplyr::n(),
    has_LN_here = any(tissue == "LN"),
    has_PT_here = any(tissue == "PT"),
    ccall_dom   = pick_dominant_isotype(c_call),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Shape = dplyr::case_when(
      has_LN_here & has_PT_here ~ "Both",
      has_PT_here               ~ "PT",
      has_LN_here               ~ "LN"
    ),
    ccall_dom  = factor(ccall_dom, levels = ISOTYPE_PRIORITY),
    cell_type  = factor(cell_type, levels = CLUSTER_ORDER),
    Count_plot = pmin(pmax(Count, 1L), max(Count))
  ) %>%
  dplyr::filter(!is.na(Shape), !is.na(ccall_dom))

# Order clones by total cell count — largest at the top of the y axis
clone_order_tbl <- plot_df %>%
  dplyr::group_by(CTstrict) %>%
  dplyr::summarise(total_cells = sum(Count), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(total_cells))

plot_df <- plot_df %>%
  dplyr::mutate(CloneID = factor(CTstrict, levels = clone_order_tbl$CTstrict))

message("Plot data: ", dplyr::n_distinct(plot_df$CTstrict),
        " clones x ", dplyr::n_distinct(plot_df$cell_type), " clusters")


# ---- D2. Build and save the plot --------------------------------------------
message("\n>> D2: Building shared clone dot plot...")

y_lines    <- seq_along(levels(plot_df$CloneID))
max_count  <- max(plot_df$Count)
mid_count  <- round(max_count / 2)

p_shared <- ggplot(plot_df,
                   aes(x     = cell_type,
                       y     = CloneID,
                       color = ccall_dom,
                       shape = Shape,
                       size  = Count_plot)) +
  geom_hline(yintercept = y_lines,
             linetype = "dotted", color = "grey70", linewidth = 0.4) +
  geom_point(stroke = 0.8) +
  scale_x_discrete(limits = CLUSTER_ORDER, drop = FALSE) +
  scale_color_manual(values = ISOTYPE_COLORS, name = "Isotype") +
  scale_shape_manual(
    values = c(PT = 15, LN = 16, Both = 17),
    breaks = c("PT", "LN", "Both"),
    labels = c(PT = "PT", LN = "LN", Both = "LN+PT"),
    name   = "Origin"
  ) +
  scale_size_continuous(
    limits = c(1, max_count),
    range  = c(3, 8),
    breaks = c(1, mid_count, max_count),
    labels = c("1", as.character(mid_count), as.character(max_count)),
    name   = "Cells / clone / cluster"
  ) +
  theme_minimal() +
  theme(
    legend.position    = "bottom",
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  labs(
    x     = "Cell type cluster",
    y     = NULL,
    title = paste0("Shared LN/PT clones — ", PATIENT_ID)
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 4)),
    shape = guide_legend(override.aes = list(size = 4)),
    size  = guide_legend(override.aes = list(shape = 16))
  )

print(p_shared)

ggsave(file.path(PLOT_DIR, paste0("shared_clones_dotplot_", PATIENT_ID, ".pdf")),
       p_shared, width = 10, height = 8)
message("Saved: shared_clones_dotplot_", PATIENT_ID, ".pdf")

message("\nPart D complete.")
message("\n", paste(rep("=", 60), collapse = ""))
message("Done. All outputs saved to:")
message("  Tables: ", OUT_DIR)
message("  RDS:    ", RDS_DIR)
message("  Plots:  ", PLOT_DIR)
message(paste(rep("=", 60), collapse = ""))
