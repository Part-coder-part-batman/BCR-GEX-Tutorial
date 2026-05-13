# =============================================================================
# 01_BCR_Pipeline.R
# BCR QC, clonal assignment, and somatic hypermutation analysis
# =============================================================================
#
# HOW TO USE:
#   1. Fill in the file paths in the USER SETTINGS section below
#   2. Run the script (Ctrl+Shift+Enter in RStudio to source the whole file)
#
# INPUT FILES (two per sample):
#   - *_airr_db-pass.tsv         output from step 00 (Immcantation / MakeDb.py)
#   - filtered_contig_annotations.csv  from Cell Ranger, found at:
#       per_sample_outs > <sample> > vdj_b > filtered_contig_annotations.csv
#
# OUTPUT:
#   - <SAMPLE_ID>_bcr_data.tsv   cleaned BCR table ready for GEX integration
#
# Run this script once per sample. The tutorial uses two samples:
#   P1_LN  (lymph node)
#   P1_PT  (primary tumor)
#
# IMGT GERMLINE DATABASE NOTE:
#   This pipeline requires a local copy of the IMGT germline reference database
#   for germline reconstruction (Step 5). The database is NOT bundled with any
#   R package and must be obtained separately. Two options:
#
#   Option A — copy from the Docker container you used in step 00:
#     docker cp <container_id>:/usr/local/share/germlines/imgt/human/vdj C:/Users/YourName/Documents/Immcantation/imgt/human/vdj
#
#   Option B — download via the Immcantation fetch script:
#     git clone https://github.com/immcantation/immcantation
#     bash immcantation/scripts/fetch_imgtdb.sh -o germlines
#
#   Point IMGT_REF_DIR below to the human/vdj folder you obtained.
# =============================================================================


# =============================================================================
# !! ONLY EDIT THIS SECTION !!
# =============================================================================

# Path to BCR_functions.R — update to wherever you saved it
BCR_FUNCTIONS_FILE <- "C:/Users/YourName/Documents/BCR-GEX-Tutorial/scripts/BCR_functions.R"

# Sample name — used in output file names and encoded into clone IDs
# Run this script once per sample, changing SAMPLE_ID each time
SAMPLE_ID <- "MySample"

# AIRR TSV from step 00 — the file ending in _airr_db-pass.tsv
AIRR_FILE <- "C:/Users/YourName/Documents/Immcantation/MySample/results/MySample_airr_db-pass.tsv"

# Cell Ranger annotations CSV
# Found at: per_sample_outs > <sample> > vdj_b > filtered_contig_annotations.csv
ANNOTATIONS_FILE <- "C:/Users/YourName/Documents/CellRanger/MySample/filtered_contig_annotations.csv"

# Output folder — will be created if it does not exist
OUTPUT_DIR <- "C:/Users/YourName/Documents/Immcantation/MySample/Output"

# IMGT germline reference folder (see note at top of this script)
IMGT_REF_DIR <- "C:/Users/YourName/Documents/Immcantation/imgt/human/vdj"

# Number of CPU cores for germline reconstruction and mutation analysis
# Keep at 1 unless you know your machine can handle more
NPROC <- 1

# =============================================================================
# END OF EDITABLE SECTION — do not change anything below
# =============================================================================


# ---- Load functions ---------------------------------------------------------
if (!file.exists(BCR_FUNCTIONS_FILE)) {
  stop(
    "Cannot find BCR_functions.R at:\n  ", BCR_FUNCTIONS_FILE,
    "\nPlease update BCR_FUNCTIONS_FILE at the top of this script."
  )
}
source(BCR_FUNCTIONS_FILE)

suppressPackageStartupMessages({
  library(airr)
  library(alakazam)
  library(dplyr)
  library(dowser)
  library(ggplot2)
  library(readr)
  library(scoper)
  library(shazam)
  library(stringr)
})


# ---- Check all input files exist --------------------------------------------
to_check <- c(
  "AIRR TSV"           = AIRR_FILE,
  "CR annotations CSV" = ANNOTATIONS_FILE,
  "IMGT reference dir" = IMGT_REF_DIR
)
missing_files <- to_check[!file.exists(to_check)]
if (length(missing_files) > 0) {
  stop(
    "The following paths could not be found:\n  ",
    paste(names(missing_files), "->", missing_files, collapse = "\n  "),
    "\nPlease check the paths at the top of this script."
  )
}

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)


# =============================================================================
# PIPELINE — runs top to bottom when you source this file
# =============================================================================

# ---- Step 1: Load and filter ------------------------------------------------
# Reads the AIRR TSV and Cell Ranger annotations, merges UMI counts,
# appends sample ID to all barcodes, then removes non-productive sequences,
# multi-heavy-chain cells, unpaired light chains, and duplicate light chains.
message("\n>> Step 1: Loading and filtering BCR data for sample: ", SAMPLE_ID)

bcr_data <- process_bcr_data(
  airr_file        = AIRR_FILE,
  annotations_file = ANNOTATIONS_FILE,
  sample_id        = SAMPLE_ID
)
bcr_data <- filter_bcr_data(bcr_data)


# ---- Step 2: Clonal threshold -----------------------------------------------
# Calculates pairwise Hamming distances between IGH sequences and fits a
# Gaussian mixture model to find the threshold that best separates clonally
# related from unrelated sequences. Inspect both plots before proceeding.
# If the auto threshold looks wrong, you can pass a manual value to
# define_clones_basic() in the next step instead.
message("\n>> Step 2: Computing clonal threshold...")

thr_res   <- plot_clonal_thresholds(bcr_data, user_defined_threshold = NULL)
threshold <- thr_res$auto_threshold
message("   Auto threshold: ", signif(threshold, 4))

print(thr_res$user_plot)
print(thr_res$shazam_plot)


# ---- Step 3: Clonal assignment ----------------------------------------------
# Assigns BCR sequences to clones using hierarchical clustering on IGH,
# with light chain splitting to resolve ambiguous cases.
# Then replaces numeric clone IDs with informative strings encoding
# sample ID, a unique random code, clone size, and isotype composition.
# Example clone ID: P1_LN_aB3x_14_G1G2M
message("\n>> Step 3: Assigning and annotating clones...")

bcr_data <- define_clones_basic(bcr_data, threshold = threshold, sample_id = SAMPLE_ID)
bcr_data <- annotate_clone_ids(bcr_data, sample_id = SAMPLE_ID)
qc_clone_assignment(bcr_data, bcr_data)


# ---- Step 4: Clone repertoire plots -----------------------------------------
# Rank-abundance curve, clone size distribution, and clonal diversity.
# These give an overview of the repertoire before looking at individual clones.
message("\n>> Step 4: Generating clone repertoire plots...")

clone_plots <- visualize_clones(bcr_data)
print(clone_plots$rank_abundance_plot)
print(clone_plots$clone_size_plot)
print(clone_plots$diversity_plot)


# ---- Step 5: Germline reconstruction and mutation analysis ------------------
# Infers the unmutated germline sequence for each clone using the IMGT
# reference, then calculates somatic hypermutation (SHM) frequency and
# raw count across the V gene region.
message("\n>> Step 5: Reconstructing germlines and computing mutation metrics...")

mut_outputs <- reconstruct_germlines_and_mutations(
  bcr_data,
  reference_dir = IMGT_REF_DIR,
  nproc         = NPROC
)
bcr_data <- mut_outputs$results_with_mut
print(mut_outputs$mut_histogram_plot)


# ---- Step 6: Mutation frequency plots ---------------------------------------
message("\n>> Step 6: Plotting mutation frequencies...")

mf_plots <- plot_mutation_frequencies(bcr_data)
print(mf_plots$histogram)
print(mf_plots$by_subject)
print(mf_plots$by_isotype)


# ---- Step 7: Save -----------------------------------------------------------
# Saves the cleaned, clone-annotated BCR table as a TSV.
# This file is the input for step 04 (Add BCR to Seurat), after GEX QC
# and doublet removal are complete.
bcr_tsv <- file.path(OUTPUT_DIR, paste0(SAMPLE_ID, "_bcr_data.tsv"))
message("\n>> Step 7: Saving BCR table to: ", bcr_tsv)

write.table(bcr_data, file = bcr_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

message("\nDone. Results saved to: ", OUTPUT_DIR)
