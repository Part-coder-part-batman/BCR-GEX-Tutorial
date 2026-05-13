# BCR + GEX Integration Tutorial

A three-part R tutorial for single-cell BCR sequence processing and integration with gene expression data, using the Immcantation toolkit and Seurat.

The dataset comes from Nathan et al., *Immunity* (2026): "Tumor-draining lymph nodes in ovarian cancer lack germinal centers but harbor tumor-reactive memory B cells clonally linked to intra-tumoral B cells." We use a paired lymph node and primary tumor sample from the same patient. Data are deposited on Zenodo: ([https://doi.org/10.5281/zenodo.20134219](https://zenodo.org/records/20134219))

---

## What the tutorial covers

**Part 1 -- BCR sequence processing (Immcantation)**
- Docker setup and V(D)J assignment
- BCR QC and filtering
- Clonal threshold
- Clone assignment
- Germline reconstruction and SHM
- Repertoire overview

**Part 2 -- GEX analysis (Seurat)**
- GEX QC
- Doublet removal
- BCR integration
- Solo GEX analysis
- RPCA integration
- Cluster annotation

**Part 3 -- Repertoire visualization (Immcantation + scRepertoire)**
- Clonal structure visualization
- Phylogenetic trees
- Shared clone analysis

Standalone `.R` scripts are also provided for each step, for users who prefer working outside of R Markdown:

- `00_Docker_Setup_and_VDJ_Assignment.md`
- `01_BCR_Pipeline.R`
- `02_GEX_QC.R`
- `03_Doublet_Removal.R`
- `04_Add_BCR.R`
- `05_Solo_GEX.R`
- `06_Integrate_Samples.R`
- `07_Annotate_Clusters.R`
- `08_BCR_Visualization.R`
- `09_Shared_Clone_Analysis.R`

---

## Tutorial parts

| Part | File | Covers |
|------|------|--------|
| Part 1 | `BCR_GEX_Tutorial_Part1.Rmd` | BCR processing with Immcantation: V(D)J assignment, QC, clonal assignment, SHM |
| Part 2 | `BCR_GEX_Tutorial_Part2.Rmd` | GEX analysis with Seurat: QC, doublet removal, BCR integration, RPCA, cluster annotation |
| Part 3 | `BCR_GEX_Tutorial_Part3.Rmd` | Repertoire visualization: donut plots, phylogenetic trees, shared clone analysis |

Each part has a Prerequisites section listing exactly which input files are needed. The Zenodo deposit includes intermediate output files so you can start from Part 2 or Part 3 without running earlier parts.

---

## Getting started

1. Download data files from [Zenodo](https://doi.org/10.5281/zenodo.20134219)
2. Clone this repository
3. Open `BCR_GEX_Tutorial_Part1.Rmd` in RStudio and update the paths chunk at the top
4. Knit each part in order, or start from a later part using the Zenodo entry-point files

---

## Requirements

- **R 4.3+** and RStudio
- **Seurat v5** (required for Part 2 -- v4 will not work with the `IntegrateLayers()` workflow)
- **Docker Desktop** (required for Part 1) -- Docker is used to run the Immcantation container, which handles V(D)J gene assignment via IgBLAST. No local IgBLAST installation is needed; everything runs inside the container. Step-by-step setup instructions are in `00_Docker_Setup_and_VDJ_Assignment.md`.
- **IQ-TREE 2** (required for phylogenetic trees in Part 3, Section 8.4 only)

R package installation instructions are at the top of each Rmd file.

---

## External resources

The following resources are required by the tutorial but must be obtained separately:

| Resource | Used in | Where to get it |
|----------|---------|-----------------|
| Immcantation Docker image (`immcantation/suite:4.5.0`) | Part 1 | `docker pull immcantation/suite:4.5.0` -- see `00_Docker_Setup_and_VDJ_Assignment.md` |
| IMGT human VDJ germline reference | Part 1 (germline reconstruction) | Download via `fetch_imgtdb.sh` from the [Immcantation repository](https://github.com/immcantation/immcantation) -- instructions in the Part 1 Prerequisites section |
| IQ-TREE 2 binary | Part 3 (phylogenetic trees only) | [github.com/Cibiv/IQ-TREE/releases](https://github.com/Cibiv/IQ-TREE/releases) -- download the binary for your OS |

---

## Citation

Nathan et al., "Tumor-draining lymph nodes in ovarian cancer lack germinal centers but harbor tumor-reactive memory B cells clonally linked to intra-tumoral B cells," *Immunity* (2026).

