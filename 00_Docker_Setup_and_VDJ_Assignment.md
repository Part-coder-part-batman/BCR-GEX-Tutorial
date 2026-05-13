# Docker Setup and V(D)J Gene Assignment with Immcantation

This guide walks you through installing Docker, setting up the Immcantation container, and running V(D)J gene assignment on BCR sequences from Cell Ranger output.

> **Platform note:** This tutorial was developed and tested on Windows. Mac users can follow the same steps, but terminal commands and file paths may differ slightly.

> **Alternative for HPC users:** If you are working on a high-performance computing cluster, consider [nf-core/airrflow](https://nf-co.re/airrflow), a Nextflow pipeline that wraps the full Immcantation workflow and is designed for scalable, reproducible analysis across many samples. This tutorial takes a hands-on Docker approach that is more accessible for researchers working on a local machine.

---

## Part 1: Install Docker Desktop

1. Go to [https://www.docker.com/products/docker-desktop/](https://www.docker.com/products/docker-desktop/) and download the installer for your operating system.

2. Run the installer and follow the setup steps.

3. After installation, find the Docker icon in the system tray (bottom-right corner of your screen). Right-click it and check the option — if it says **Switch to Linux Containers**, click it. If it says **Switch to Desktop Containers**, you are already set up correctly.

4. Restart your PC to make sure all changes are applied.

---

## Part 2: Test Your Docker Installation

Open a terminal (Command Prompt or PowerShell on Windows) and run:

```bash
docker --version
docker info
```

If Docker is installed correctly, both commands will return version and system information without errors.

---

## Part 3: Pull the Immcantation Docker Image

Download the Immcantation suite image:

```bash
docker pull immcantation/suite:4.5.0
```

This may take a few minutes depending on your internet connection. You only need to do this once.

---

## Part 4: Prepare Your Data

You will need the Cell Ranger VDJ output file `filtered_contig.fasta`. This file is located inside your Cell Ranger output folder at:

```
per_sample_outs > vdj_b > filtered_contig.fasta
```

Create a dedicated project folder on your computer where all Immcantation input and output will be stored. For example:

```
C:\Users\YourName\Documents\Immcantation\ProjectName
```

Copy your `filtered_contig.fasta` file into this folder.

---

## Part 5: Run the Docker Container

Run the following command to start the container and mount your project folder into it:

```bash
docker run -it -v "C:/Users/YourName/Documents/Immcantation/ProjectName:/data:z" immcantation/suite:4.5.0 bash
```

Replace `C:/Users/YourName/Documents/Immcantation/ProjectName` with the actual path to your project folder. This mounts your folder to `/data` inside the container, meaning any files you write to `/data` inside the container will appear in your project folder on your PC.

Once inside the container, verify your file is accessible:

```bash
ls /data
head /data/filtered_contig.fasta
```

---

## Part 6: V(D)J Gene Assignment with IgBLAST

Run the following command to assign V, D, and J genes to each contig:

```bash
AssignGenes.py igblast \
   -s /data/filtered_contig.fasta \
   -b /usr/local/share/igblast \
   --organism human \
   --loci ig \
   --format blast \
   --outdir /data/results \
   --outname Project_BCR
```

**Argument explanation:**

| Argument | Description |
|----------|-------------|
| `-s` | Path to the input FASTA file |
| `-b` | Path to the IgBLAST database inside the container |
| `--organism` | Organism (use `human` for human samples) |
| `--loci` | Locus to analyze (`ig` for B cell receptor) |
| `--format` | Output format (use `blast`) |
| `--outdir` | Directory where results will be saved |
| `--outname` | Prefix for all output files |

---

## Part 7: Convert to AIRR Format

Convert the IgBLAST output into a standardized AIRR-compliant TSV table with one row per contig:

```bash
MakeDb.py igblast \
   -i /data/results/Project_BCR_igblast.fmt7 \
   -s /data/filtered_contig.fasta \
   -r /usr/local/share/germlines/imgt/human/vdj/ \
   --outdir /data/results \
   --outname Project_BCR_airr
```

**Argument explanation:**

| Argument | Description |
|----------|-------------|
| `-i` | IgBLAST output file from the previous step |
| `-s` | Original input FASTA file |
| `-r` | Path to IMGT germline reference database inside the container |
| `--outdir` | Directory where results will be saved |
| `--outname` | Prefix for all output files |

---

## Part 8: Verify the Output

Check that the output files were created:

```bash
ls /data/results
```

Inspect the first few lines of the AIRR table:

```bash
head /data/results/Project_BCR_airr_db-pass.tsv
```

You should see a tab-separated table with one contig per row and columns for sequence ID, V/D/J gene calls, junction sequence, and other annotations.

---

## Part 9: Access Your Results

Exit the container when done:

```bash
exit
```

Your results are saved in your project folder on your PC:

```
C:\Users\YourName\Documents\Immcantation\ProjectName\results\
```

The key output file is:

- `Project_BCR_airr_db-pass.tsv` — AIRR-formatted table ready for downstream BCR QC and clonal analysis in R

---

## Next Step

With the AIRR table in hand, you are ready to move to R for BCR quality control, clonal assignment, and integration with gene expression data. See `01_BCR_Pipeline.R` for the next step.
