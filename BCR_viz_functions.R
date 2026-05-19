# =============================================================================
# BCR_viz_functions.R
# Helper functions for BCR clonality visualization and phylogenetic trees
#
# These functions are sourced by BCR_GEX_Tutorial_Part3.Rmd -- do not run this
# file directly. All required packages are loaded by the tutorial script.
#
# INPUT: The annotated BCR data frame produced in Part 3 after cell_type
# metadata has been joined from the integrated Seurat object. This data frame
# must contain one row per BCR contig (heavy and light chains), with the
# following key columns:
#   - cell_id       : barcode with sample suffix (e.g. "ACGT..._P1_LN")
#   - clone_id      : informative clone ID from annotate_clone_ids()
#   - clone_count   : number of cells in this clone
#   - c_call        : isotype (e.g. "IGHG1", "IGHM")
#   - locus         : chain type ("IGH", "IGK", "IGL")
#   - cell_type     : cluster annotation from the Seurat object (NA for cells
#                     filtered out during GEX QC -- these are retained in the
#                     BCR data but excluded from cluster-level visualizations)
#   - sample_id     : sample of origin (e.g. "P1_LN", "P1_PT")
#   - sequence      : nucleotide sequence (required for trees)
#   - germline_alignment_d_mask : germline sequence (required for trees)
#
# Functions in this file (in order of use):
#   1. plot_combined_donut()       -- all clones in one donut, colored by isotype
#   2. plot_isotype_donuts()       -- one donut per isotype
#   3. plot_cluster_donuts()       -- one donut per cell type cluster
#   4. build_bcr_trees()           -- phylogenetic trees for top N clones or a
#                                    user-supplied list of clone IDs
#
# External dependency for trees:
#   IQ-TREE 2 must be installed and findable by dowser::getTrees().
#   See the Prerequisites section in the tutorial for installation instructions.
# =============================================================================


# =============================================================================
# Shared constants
# =============================================================================

# Standard isotype color palette used consistently across all plots.
# Override by passing your own named vector to the color_map argument.
BCR_ISOTYPE_COLORS <- c(
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

ISOTYPE_PALETTES <- c(
  IGHD  = "Reds",
  IGHM  = "Blues",
  IGHA1 = "BrBG",
  IGHA2 = "RdPu",
  IGHG1 = "Greens",
  IGHG2 = "Purples",
  IGHG3 = "Oranges",
  IGHG4 = "YlOrBr",
  IGHE  = "PuBuGn"
)


# =============================================================================
# Internal helper: build donut plot data from a clone/isotype data frame
# =============================================================================

.make_donut_data <- function(bcr_df) {
  # bcr_df must have columns: cell_id, clone_id, c_call (isotype)
  # Returns a data frame ready for geom_rect donut plotting, with columns:
  #   clone_id, clone_count, isotype, is_singleton, fraction, ymin, ymax, group

  # Filter to heavy chains only -- one row per cell with the correct isotype.
  # Using all contigs (heavy + light) causes each cell to appear twice with
  # different c_call values, inflating clone counts and corrupting the donut.
  df <- bcr_df %>%
    dplyr::filter(locus == "IGH", !is.na(clone_id), !is.na(c_call), c_call != "") %>%
    dplyr::select(cell_id, clone_id, isotype = c_call) %>%
    dplyr::distinct()

  if (nrow(df) == 0) return(NULL)

  expanded <- df %>%
    dplyr::group_by(clone_id, isotype) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::summarise(clone_count = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(is_singleton = FALSE)

  # All singletons collapsed into ONE slice -- keeping them separate produces
  # a white donut with visible slice lines for every singleton clone, which
  # is visually misleading (looks like expanded clones but they are all white)
  n_singletons <- df %>%
    dplyr::group_by(clone_id) %>%
    dplyr::filter(dplyr::n() == 1) %>%
    nrow()

  if (n_singletons > 0) {
    singleton <- data.frame(
      clone_id     = "Singleton",
      isotype      = "Singleton",
      clone_count  = n_singletons,
      is_singleton = TRUE,
      stringsAsFactors = FALSE
    )
  } else {
    singleton <- data.frame(
      clone_id = character(0), isotype = character(0),
      clone_count = integer(0), is_singleton = logical(0)
    )
  }

  out <- dplyr::bind_rows(expanded, singleton) %>%
    dplyr::arrange(is_singleton, dplyr::desc(clone_count)) %>%
    dplyr::mutate(
      group    = factor(paste0("clone", dplyr::row_number()),
                        levels = paste0("clone", seq_len(dplyr::n()))),
      fraction = clone_count / sum(clone_count),
      ymax     = cumsum(fraction),
      ymin     = c(0, head(ymax, -1))
    )

  return(out)
}


# =============================================================================
# Internal helper: render a donut ggplot from pre-built donut data
# =============================================================================

.render_donut <- function(donut_data, center_label, title = NULL,
                          fill_values = NULL) {
  # fill_values: named character vector of hex colors, names matching group levels.
  # If NULL, fills are drawn from the isotype palette (for combined donuts).

  if (is.null(fill_values)) {
    # Color by isotype: expanded clones get isotype color, singletons white
    fill_values <- setNames(
      ifelse(donut_data$is_singleton, "white",
             BCR_ISOTYPE_COLORS[donut_data$isotype]),
      donut_data$group
    )
    fill_values[is.na(fill_values)] <- "grey80"
  }

  p <- ggplot2::ggplot(
    donut_data,
    ggplot2::aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = group)
  ) +
    ggplot2::geom_rect(linewidth = 0.2, color = "black") +
    ggplot2::coord_polar(theta = "y", direction = 1, start = 1e-10) +
    ggplot2::xlim(c(2, 4)) +
    ggplot2::theme_void() +
    ggplot2::annotate("text", x = 2, y = 0, label = center_label,
                      size = 5, hjust = 0.5, vjust = 0.5) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_fill_manual(values = fill_values)

  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 11))
  }

  return(p)
}


# =============================================================================
# 1. plot_combined_donut
# =============================================================================

#' Plot a single donut showing all clones across all isotypes
#'
#' Each slice represents one clone (or the aggregated singleton pool).
#' Expanded clones are colored by their isotype using the standard palette;
#' singletons are shown in white.
#'
#' @param bcr_data Annotated BCR data frame. Must contain columns cell_id,
#'   clone_id, and c_call (isotype). Typically the long-format BCR data with
#'   cell_type joined from the Seurat object.
#' @param sample_name Character string used in plot titles and center labels.
#' @param color_map Optional named character vector of isotype colors. Defaults
#'   to the standard BCR_ISOTYPE_COLORS palette defined in this file.
#' @return Named list with one ggplot object: all_clones.
#'
#' @examples
#' donuts <- plot_combined_donut(bcr_data, sample_name = "P1_LN")
#' print(donuts$all_clones)

plot_combined_donut <- function(bcr_data, sample_name, color_map = NULL) {

  if (!is.null(color_map)) {
    BCR_ISOTYPE_COLORS[names(color_map)] <- color_map
  }

  message("Building combined donut for: ", sample_name)

  donut_data <- .make_donut_data(bcr_data)
  if (is.null(donut_data) || nrow(donut_data) == 0) {
    warning("No valid clone/isotype data found for sample: ", sample_name)
    return(NULL)
  }

  total_cells    <- sum(donut_data$clone_count)
  expanded_cells <- sum(donut_data$clone_count[!donut_data$is_singleton])
  n_expanded     <- sum(!donut_data$is_singleton)

  label_all <- paste0(sample_name, "\nTotal: ", total_cells,
                      "\nExpanded: ", expanded_cells)
  p_all <- .render_donut(donut_data, center_label = label_all,
                         title = paste0(sample_name, " -- All clones"))

  message("  Total cells: ", total_cells, " | Expanded clones: ", n_expanded)
  return(list(all_clones = p_all))
}


# =============================================================================
# 2. plot_isotype_donuts
# =============================================================================

#' Plot one donut per isotype, each showing clonal structure within that isotype
#'
#' For each isotype present in the data, generates a donut where each slice is
#' one expanded clone (colored by shade within the isotype palette) or the
#' aggregated singletons (white).
#'
#' @param bcr_data Annotated BCR data frame. Must contain cell_id, clone_id,
#'   and c_call.
#' @param sample_name Character string used in plot titles.
#' @param isotype_order Optional character vector specifying the order of
#'   isotypes to process. Defaults to the standard BCR isotype order.
#' @param min_cells Minimum number of cells required for an isotype to be
#'   plotted. Default: 5.
#' @return Named list. Each element is named by isotype and contains a nested
#'   list with one ggplot object: all_clones.
#'
#' @examples
#' iso_plots <- plot_isotype_donuts(bcr_data, sample_name = "P1_LN")
#' print(iso_plots$IGHG1$all_clones)

plot_isotype_donuts <- function(bcr_data, sample_name,
                                isotype_order = names(BCR_ISOTYPE_COLORS),
                                min_cells = 5) {

  message("Building per-isotype donuts for: ", sample_name)

  present_isotypes <- bcr_data %>%
    dplyr::filter(!is.na(c_call), c_call != "", locus == "IGH") %>%
    dplyr::pull(c_call) %>%
    unique()

  # Preserve standard order where possible
  ordered_iso <- c(
    intersect(isotype_order, present_isotypes),
    setdiff(present_isotypes, isotype_order)
  )

  results <- list()

  for (iso in ordered_iso) {

    iso_data <- bcr_data %>% dplyr::filter(c_call == iso)
    n_cells  <- nrow(iso_data)

    if (n_cells < min_cells) {
      message("  Skipping ", iso, " -- only ", n_cells, " cells (min_cells = ", min_cells, ")")
      next
    }

    donut_data <- .make_donut_data(iso_data)
    if (is.null(donut_data) || nrow(donut_data) == 0) next

    n_expanded <- sum(!donut_data$is_singleton)
    total      <- sum(donut_data$clone_count)

    # Colors: shades within the isotype palette for expanded clones, white for singletons
    palette_name <- ISOTYPE_PALETTES[iso]
    if (is.na(palette_name)) palette_name <- "Greys"

    if (palette_name == "BrBG") {
      exp_colors <- colorRampPalette(RColorBrewer::brewer.pal(11, "BrBG")[1:5])(max(1, n_expanded))
    } else {
      exp_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, palette_name)))(max(1, n_expanded))
    }

    all_colors <- c(exp_colors, rep("white", sum(donut_data$is_singleton)))
    # Use as.character() to avoid dropped factor levels after filtering
    fill_values <- setNames(all_colors, as.character(donut_data$group))

    label_all <- paste0(iso, "\n", total, " cells")
    p_all <- .render_donut(donut_data, center_label = label_all,
                           title = paste0(sample_name, " -- ", iso),
                           fill_values = fill_values)

    message("  ", iso, ": ", total, " cells | ", n_expanded, " expanded clones")
    results[[iso]] <- list(all_clones = p_all)
  }

  message("plot_isotype_donuts complete. Isotypes plotted: ",
          paste(names(results), collapse = ", "))
  return(results)
}


# =============================================================================
# 3. plot_cluster_donuts
# =============================================================================

#' Plot one donut per cell type cluster, showing clonal structure within each
#'
#' Loops through all cell type clusters present in the data and generates a
#' combined donut for each cluster using plot_combined_donut(). Clusters with
#' fewer than min_cells cells are skipped. Cells with NA cell_type (BCR cells
#' absent from the Seurat object due to GEX QC filtering) are excluded.
#'
#' @param bcr_data Annotated BCR data frame. Must contain cell_id, clone_id,
#'   c_call, and a cell type column (default: "cell_type").
#' @param sample_name Character string used in plot titles.
#' @param cell_type_col Name of the column containing cell type labels.
#'   Default: "cell_type".
#' @param min_cells Minimum number of cells required per cluster to generate
#'   a plot. Default: 10.
#' @return Named list. Each element is named by cluster label and contains a
#'   nested list with one ggplot object: all_clones.
#'
#' @examples
#' cluster_plots <- plot_cluster_donuts(bcr_data, sample_name = "P1_LN")
#' print(cluster_plots$GC$all_clones)

plot_cluster_donuts <- function(bcr_data, sample_name,
                                cell_type_col = "cell_type",
                                min_cells = 10) {

  if (!cell_type_col %in% colnames(bcr_data)) {
    stop("Column '", cell_type_col, "' not found in bcr_data. ",
         "Check that cell_type metadata has been joined from the Seurat object.")
  }

  message("Building per-cluster donuts for: ", sample_name)

  # Exclude NA cell_type rows -- these are BCR cells filtered out during GEX QC
  # and have no valid cluster assignment
  clusters <- bcr_data %>%
    dplyr::filter(!is.na(.data[[cell_type_col]])) %>%
    dplyr::pull(.data[[cell_type_col]]) %>%
    unique() %>%
    sort()

  message("  Clusters found: ", paste(clusters, collapse = ", "))

  results <- list()

  for (cl in clusters) {

    cl_data <- bcr_data %>%
      dplyr::filter(.data[[cell_type_col]] == cl)

    n_cells <- nrow(cl_data)

    if (n_cells < min_cells) {
      message("  Skipping '", cl, "' -- only ", n_cells,
              " cells (min_cells = ", min_cells, ")")
      next
    }

    tryCatch({
      plots <- plot_combined_donut(cl_data,
                                   sample_name = paste0(sample_name, " / ", cl))
      results[[cl]] <- list(all_clones = plots$all_clones)
      message("  Done: ", cl, " (", n_cells, " cells)")
    }, error = function(e) {
      message("  Error for cluster '", cl, "': ", e$message)
    })
  }

  message("plot_cluster_donuts complete. Clusters plotted: ",
          paste(names(results), collapse = ", "))
  return(results)
}



# =============================================================================
# Internal helper: rank-based tip size mapping for tree plots
# =============================================================================

# Maps seq_group_count values to point sizes by rank (not by count value).
# Sizes are evenly spaced from size_min to size_max across unique count ranks.
# Legend shows at most max_legend entries: all if <= max_legend unique counts,
# otherwise min, max, and max_legend-2 evenly spaced from actual values between.

.make_tip_sizes <- function(counts, size_min = 2, size_max = 6, max_legend = 5) {
  actual_counts <- sort(unique(counts))
  n_unique      <- length(actual_counts)

  if (n_unique == 1) {
    size_map <- setNames(size_min, as.character(actual_counts))
  } else {
    sizes    <- seq(size_min, size_max, length.out = n_unique)
    size_map <- setNames(sizes, as.character(actual_counts))
  }

  pt_size <- size_map[as.character(counts)]

  if (n_unique <= max_legend) {
    legend_counts <- actual_counts
  } else {
    inner_idx     <- round(seq(2, n_unique - 1, length.out = max_legend - 2))
    inner_idx     <- unique(pmax(2, pmin(n_unique - 1, inner_idx)))
    legend_counts <- actual_counts[c(1, inner_idx, n_unique)]
    legend_counts <- unique(legend_counts)
  }
  legend_sizes <- size_map[as.character(legend_counts)]

  list(
    pt_size       = unname(pt_size),
    legend_counts = legend_counts,
    legend_sizes  = unname(legend_sizes)
  )
}


# =============================================================================
# 4. build_bcr_trees
# =============================================================================

#' Build and plot BCR phylogenetic trees for expanded clones
#'
#' Uses the dowser package to build maximum likelihood phylogenetic trees
#' rooted on the unmutated germline sequence. Tips are colored by isotype and
#' annotated with cell type labels. Tip size reflects how many cells share the
#' same VDJ sequence (V-start to J-end, from AIRR coordinates).
#'
#' When top_n is used, a candidate pool of up to 3x top_n clones is tried in
#' size order. Skipped clones (identical VDJ, too small after filtering) are
#' replaced by the next available clone so that up to top_n plots are returned.
#'
#' Prerequisites -- IQ-TREE 2:
#'   dowser::getTrees() requires IQ-TREE 2 to be installed on your system.
#'   1. Download from: https://github.com/Cibiv/IQ-TREE/releases
#'   2. Unzip and note the full path to the executable
#'      (e.g. "C:/tools/iqtree2/bin/iqtree2.exe" on Windows,
#'             "/usr/local/bin/iqtree2" on macOS/Linux)
#'   3. Pass this path to exec, or add the bin folder to your system PATH.
#'
#' @param bcr_data Annotated BCR data frame (output of Part 3 join step).
#'   Must contain: clone_id, cell_id, c_call, locus, sequence,
#'   germline_alignment_d_mask, v_sequence_start, j_sequence_end, and the
#'   column specified by cell_type_col.
#' @param cluster_colors Named character vector mapping cell type labels to
#'   hex colors. Names must match cell_type_col values exactly.
#' @param top_n Integer. Target number of successful tree plots. Default: 5.
#'   Ignored if clone_ids is supplied.
#' @param clone_ids Optional character vector of specific clone IDs. If
#'   supplied, top_n and fill mode are ignored.
#' @param cell_type_col Name of the cell type column. Default: "cell_type".
#' @param exec Path to IQ-TREE 2 executable, or "iqtree2" if on system PATH.
#' @param min_clone_size Minimum paired cells required to attempt a tree.
#'   Default: 3.
#' @param size_range Numeric vector of length 2: c(min_size, max_size) for
#'   tip point sizes. Default: c(2, 5).
#' @return Named list of ggplot objects. Skipped clones are NULL.
#'
#' @examples
#' tree_plots <- build_bcr_trees(
#'   bcr_data       = bcr_annotated[["P1_LN"]],
#'   cluster_colors = cluster_colors,
#'   top_n          = 5,
#'   exec           = "iqtree2"
#' )
#' print(tree_plots[["P1_LN_4aFs_6_G2G4"]])

build_bcr_trees <- function(bcr_data,
                             cluster_colors,
                             top_n          = 5,
                             clone_ids      = NULL,
                             cell_type_col  = "cell_type",
                             exec           = "iqtree2",
                             min_clone_size = 3,
                             size_range     = c(2, 5)) {

  # ---- Package checks -------------------------------------------------------
  for (pkg in c("dowser", "ggtree", "ape")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required. Install with: ",
           if (pkg == "ggtree") "BiocManager::install('ggtree')"
           else paste0("install.packages('", pkg, "')"))
    }
  }

  # ---- Column checks --------------------------------------------------------
  required_cols <- c("clone_id", "cell_id", "c_call", "locus",
                     "sequence", "germline_alignment_d_mask", cell_type_col)
  missing_cols <- setdiff(required_cols, colnames(bcr_data))
  if (length(missing_cols) > 0) {
    stop("The following required columns are missing from bcr_data:\n  ",
         paste(missing_cols, collapse = ", "), "\n",
         "Make sure germline reconstruction (Part 1) and cell_type joining ",
         "(Part 3) have both been completed.")
  }

  has_vdj_coords <- all(c("v_sequence_start", "j_sequence_end") %in% colnames(bcr_data))
  if (!has_vdj_coords) {
    warning("Columns v_sequence_start and j_sequence_end not found. ",
            "Tip size grouping will use the full sequence column instead of VDJ region. ",
            "Re-run MakeDb.py (Step 00, see 00_Docker_Setup_and_VDJ_Assignment.md) to obtain VDJ coordinates.")
  }

  # ---- Select clones --------------------------------------------------------
  if (!is.null(clone_ids)) {
    not_found <- setdiff(clone_ids, unique(bcr_data$clone_id))
    if (length(not_found) > 0)
      warning("Clone IDs not found and will be skipped: ",
              paste(not_found, collapse = ", "))
    target_clones <- intersect(clone_ids, unique(bcr_data$clone_id))
    fill_mode     <- FALSE
  } else {
    clone_sizes <- bcr_data %>%
      dplyr::filter(locus == "IGH") %>%
      dplyr::group_by(clone_id) %>%
      dplyr::summarise(n_cells = dplyr::n_distinct(cell_id), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(n_cells))
    pool_size     <- min(nrow(clone_sizes), top_n * 3)
    target_clones <- clone_sizes$clone_id[seq_len(pool_size)]
    fill_mode     <- TRUE
    message("Attempting top ", top_n, " trees from a pool of ",
            pool_size, " candidates (skipped clones will be replaced):")
    for (i in seq_len(min(pool_size, 10))) {
      sz <- clone_sizes$n_cells[i]
      message("  ", i, ". ", target_clones[i], "  (", sz, " cells)")
    }
    if (pool_size > 10) message("  ... and ", pool_size - 10, " more")
  }

  if (length(target_clones) == 0) stop("No valid clones to build trees for.")

  # ---- Build one tree per clone ---------------------------------------------
  tree_plots <- list()

  for (cid in target_clones) {

    message("\n", paste(rep("-", 60), collapse = ""))
    message("Building tree for clone: ", cid)

    tryCatch({

      df <- bcr_data %>% dplyr::filter(clone_id == cid)

      n_cells <- dplyr::n_distinct(df$cell_id[df$locus == "IGH"])
      if (n_cells < min_clone_size) {
        message("  Skipping -- only ", n_cells, " heavy chain cells")
        tree_plots[[cid]] <- NULL
        next
      }

      # Rename c_call -> c_gene; add cell_type_display
      df <- df %>%
        dplyr::rename(c_gene = c_call) %>%
        dplyr::mutate(
          cell_type_display = dplyr::case_when(
            is.na(.data[[cell_type_col]]) ~ NA_character_,
            .data[[cell_type_col]] == ""  ~ NA_character_,
            TRUE ~ as.character(.data[[cell_type_col]])
          )
        )

      # Correct locus from c_gene (more reliable than original locus column)
      df <- df %>%
        dplyr::mutate(locus = dplyr::case_when(
          stringr::str_detect(c_gene, "^IGH") ~ "IGH",
          stringr::str_detect(c_gene, "^IGK") ~ "IGK",
          stringr::str_detect(c_gene, "^IGL") ~ "IGL",
          TRUE ~ locus
        ))

      # Remove contigs with no isotype call
      df <- df %>% dplyr::filter(!is.na(c_gene), c_gene != "")

      has_heavy   <- any(df$locus == "IGH")
      has_light   <- any(df$locus %in% c("IGK", "IGL"))
      chain_param <- if (has_heavy && has_light) "HL" else if (has_heavy) "H" else "L"
      message("  Chain type: ", chain_param, " | Sequences: ", nrow(df))

      # For HL trees: keep only cells with both heavy AND light chain present
      if (chain_param == "HL") {
        cells_h    <- df$cell_id[df$locus == "IGH"]
        cells_l    <- df$cell_id[df$locus %in% c("IGK", "IGL")]
        paired     <- intersect(cells_h, cells_l)
        n_unpaired <- length(cells_h) - length(paired)
        if (n_unpaired > 0)
          message("  Removed ", n_unpaired,
                  " cells: heavy chain present but no paired light chain")
        df <- df %>% dplyr::filter(cell_id %in% paired)
      }

      # Re-check size after all filtering
      n_cells_post <- dplyr::n_distinct(df$cell_id[df$locus == "IGH"])
      if (n_cells_post < min_clone_size) {
        message("  Skipping -- only ", n_cells_post,
                " cells remain after filtering (min_clone_size = ",
                min_clone_size, ")")
        tree_plots[[cid]] <- NULL
        next
      }
      message("  Building tree on ", n_cells_post, " cells",
              if (n_cells_post < n_cells) paste0(" (", n_cells - n_cells_post,
              " removed: no isotype call or unpaired light chain)") else "")

      # Compute VDJ_DNA_sequence and seq_group_count
      # VDJ_DNA_sequence: V-start to J-end trim (AIRR coordinates from MakeDb.py)
      # seq_group: VDJ sequence + isotype + cell_type -- cells sharing all three
      # collapse to one tree tip; tip size = count of cells in that group.
      igh_df <- df %>%
        dplyr::filter(locus == "IGH") %>%
        dplyr::mutate(
          VDJ_DNA_sequence = if (has_vdj_coords) {
            dplyr::if_else(
              !is.na(v_sequence_start) & !is.na(j_sequence_end) &
                v_sequence_start > 0   & j_sequence_end > 0,
              substr(sequence, v_sequence_start, j_sequence_end),
              sequence
            )
          } else { sequence },
          # Collapsing key: VDJ sequence + isotype only.
          # Cell type is NOT part of the key -- two cells with identical VDJ
          # and identical isotype are making the same antibody regardless of
          # their GEX cluster. The dominant cell type across cells sharing a
          # seq_group is assigned as the tip label below.
          seq_group = paste(VDJ_DNA_sequence, c_gene, sep = "_")
        )

      # Skip clones where all heavy VDJ sequences are identical --
      # IQ-TREE cannot build a tree from one unique sequence.
      # This represents clonal expansion with no further SHM (biologically
      # valid finding, but there is no phylogenetic signal to visualize).
      n_unique_heavy <- dplyr::n_distinct(igh_df$VDJ_DNA_sequence)
      if (n_unique_heavy < 2) {
        message("  Skipping -- all ", n_cells_post, " cells share the same ",
                "heavy chain VDJ sequence. No phylogenetic signal to build ",
                "a tree from (clonal expansion with no further SHM).")
        tree_plots[[cid]] <- NULL
        next
      }

      # For each seq_group, count cells and find the dominant cell type.
      # Dominant = most frequent non-NA cell_type among cells in that group.
      # This becomes the tip label. If all cells in a group have NA cell_type,
      # the tip is shown without a label rectangle.
      seq_group_summary <- igh_df %>%
        dplyr::group_by(seq_group) %>%
        dplyr::summarise(
          seq_group_count   = dplyr::n(),
          dominant_celltype = {
            ct <- cell_type_display[!is.na(cell_type_display)]
            if (length(ct) == 0) NA_character_
            else names(sort(table(ct), decreasing = TRUE))[1]
          },
          .groups = "drop"
        )

      # Collapse to one representative cell per seq_group BEFORE formatClones.
      # Without this, formatClones receives multiple cells with identical VDJ
      # sequences and internally keeps an arbitrary subset as tips. When those
      # tips are joined back to seq_counts by sequence_id, multiple tips from
      # the same seq_group each receive the full group count, making the plot
      # show inflated tip sizes (e.g. two size-4 dots instead of one).
      # Solution: feed formatClones exactly one cell per seq_group. Each tip in
      # the resulting tree then has a unique, unambiguous seq_group_count.
      rep_cells <- igh_df %>%
        dplyr::left_join(
          seq_group_summary %>% dplyr::select(seq_group, seq_group_count,
                                              dominant_celltype),
          by = "seq_group"
        ) %>%
        dplyr::group_by(seq_group) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup() %>%
        dplyr::select(cell_id, seq_group_count, dominant_celltype)

      # Keep all loci rows (heavy + light) for representative cells only
      df_fmt <- df %>%
        dplyr::filter(cell_id %in% rep_cells$cell_id)

      # seq_counts: one row per representative sequence_id -> clean 1-to-1 join
      seq_counts <- df_fmt %>%
        dplyr::filter(locus == "IGH") %>%
        dplyr::left_join(rep_cells, by = "cell_id") %>%
        dplyr::select(sequence_id, seq_group_count, dominant_celltype)

      # Format clones for dowser using collapsed data
      has_cell_type <- any(!is.na(df_fmt$cell_type_display))
      trait_cols    <- if (has_cell_type) c("c_gene", "cell_type_display") else "c_gene"

      if (chain_param == "HL") {
        df_fmt <- df_fmt %>%
          dplyr::mutate(subgroup = dplyr::case_when(
            locus == "IGH" ~ "heavy",
            locus == "IGK" ~ "kappa",
            locus == "IGL" ~ "lambda",
            TRUE           ~ "unknown"
          ))
        clones_fmt <- dowser::formatClones(
          df_fmt, traits = trait_cols, locus = "locus",
          chain = chain_param, subgroup = "subgroup", minseq = 1
        )
      } else {
        clones_fmt <- dowser::formatClones(
          df_fmt, traits = trait_cols, locus = "locus",
          chain = chain_param, minseq = 1
        )
      }

      if (is.null(clones_fmt) || nrow(clones_fmt) == 0) {
        message("  Skipping -- formatClones returned no data.")
        tree_plots[[cid]] <- NULL
        next
      }

      # Build trees with IQ-TREE 2
      message("  Running IQ-TREE 2...")
      trees <- dowser::getTrees(clones_fmt, exec = exec)
      trees <- dowser::scaleBranches(trees, edge_type = "mutations")
      message("  Tree built successfully.")

      # plotTrees wrapped separately -- degenerate topologies skip cleanly
      base_plot <- tryCatch(
        dowser::plotTrees(trees, tips = "c_gene", scale = FALSE)[[1]],
        error = function(e) {
          message("  plotTrees failed: ", e$message, " -- skipping.")
          NULL
        }
      )
      if (is.null(base_plot)) { tree_plots[[cid]] <- NULL; next }

      tree_data <- base_plot$data

      tip_data <- tree_data %>%
        dplyr::filter(isTip == TRUE) %>%
        dplyr::arrange(y) %>%
        dplyr::left_join(
          seq_counts %>%
            dplyr::select(sequence_id, seq_group_count, dominant_celltype) %>%
            dplyr::distinct(),
          by = c("label" = "sequence_id")
        ) %>%
        dplyr::mutate(seq_group_count = tidyr::replace_na(seq_group_count, 1L))

      # Rank-based tip sizes: evenly spaced from size_range[1] to size_range[2]
      # by count rank. Legend shows up to 5 entries using actual count values.
      sizing           <- .make_tip_sizes(tip_data$seq_group_count,
                                          size_min = size_range[1],
                                          size_max = size_range[2])
      tip_data$pt_size <- sizing$pt_size
      legend_counts    <- sizing$legend_counts
      legend_sizes     <- sizing$legend_sizes

      # Only draw cell type label rectangles for tips with known dominant
      # cell type and seq_group_count >= 2 (singleton tips get point only)
      tip_data_labeled <- tip_data %>%
        dplyr::filter(!is.na(dominant_celltype), seq_group_count >= 2)

      p <- ggplot2::ggplot(tree_data) +
        ggtree::geom_tree(ggplot2::aes(x = x, y = y),
                          linewidth = 0.5, color = "black") +
        ggplot2::geom_vline(xintercept = seq(20, 80, by = 20),
                            linetype = "dashed", color = "grey60",
                            linewidth = 0.3) +
        ggplot2::geom_point(
          data = tip_data,
          ggplot2::aes(x = x, y = y, color = c_gene, size = pt_size),
          stroke = 0.5
        ) +
        ggplot2::geom_rect(
          data = tip_data_labeled,
          ggplot2::aes(xmin = x + 0.5, xmax = x + 11,
                       ymin = y - 0.45, ymax = y + 0.45,
                       fill = dominant_celltype),
          color = "white", linewidth = 0.4
        ) +
        ggplot2::geom_text(
          data = tip_data_labeled,
          ggplot2::aes(x = x + 5.75, y = y, label = dominant_celltype),
          size = 3.2, color = "white", fontface = "bold"
        ) +
        ggplot2::scale_color_manual(
          values = BCR_ISOTYPE_COLORS,
          name   = "Isotype",
          guide  = ggplot2::guide_legend(override.aes = list(size = 4))
        ) +
        ggplot2::scale_fill_manual(
          values   = cluster_colors,
          name     = "Cell type
(dominant)",
          na.value = "grey85"
        ) +
        ggplot2::scale_size_identity(
          name   = "Cells (identical VDJ)",
          breaks = legend_sizes,
          labels = as.character(legend_counts),
          guide  = ggplot2::guide_legend(
            override.aes = list(color = "black")
          )
        ) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title    = ggplot2::element_text(hjust = 0.5, size = 13,
                                                face = "bold"),
          plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10,
                                                color = "grey40"),
          legend.position = "right",
          panel.grid      = ggplot2::element_blank(),
          axis.text.y     = ggplot2::element_blank(),
          axis.ticks.y    = ggplot2::element_blank(),
          axis.title.y    = ggplot2::element_blank(),
          axis.text.x     = ggplot2::element_text(size = 10),
          axis.title.x    = ggplot2::element_text(size = 12),
          plot.margin     = ggplot2::margin(10, 40, 10, 20)
        ) +
        ggplot2::ggtitle(
          label    = paste0("Clone: ", cid),
          subtitle = paste0(n_cells_post, " cells in tree | chain: ", chain_param)
        ) +
        ggplot2::labs(
          x       = "Mutations from germline root",
          caption = "Tip size = cells with identical VDJ sequence | Tips colored by isotype"
        )

      tree_plots[[cid]] <- p
      message("  Plot created.")

      # In fill_mode, stop as soon as we have top_n successful plots
      if (fill_mode && sum(!sapply(tree_plots, is.null)) >= top_n) {
        message("  Reached target of ", top_n, " successful trees. Stopping.")
        break
      }

    }, error = function(e) {
      message("  ERROR for clone '", cid, "': ", e$message)
      tree_plots[[cid]] <<- NULL
    })
  }

  n_success <- sum(!sapply(tree_plots, is.null))
  message("\n", paste(rep("=", 60), collapse = ""))
  message("build_bcr_trees complete.")
  message("Trees built successfully: ", n_success, " / ",
          if (fill_mode) top_n else length(target_clones))
  skipped <- names(tree_plots)[sapply(tree_plots, is.null)]
  if (length(skipped) > 0)
    message("Skipped: ", paste(skipped, collapse = ", "))

  # Return only non-NULL plots
  tree_plots <- Filter(Negate(is.null), tree_plots)
  return(tree_plots)
}
