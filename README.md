<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/github_username/repo_name">
    <img src="images/logo.png" alt="Logo" width="80" height="80">
  </a>

  <h3 align="center">Single Cell Multiome Analysis</h3>

  <p align="center">
    TF-induced expresssion profile
    <br />
  </p>
</div>

## Table of Contents
<details open>
<summary><b>Catalogue</b></summary>

1. [Description](#description)
2. [Dependencies](#dependencies)
3. [File Format](#file-format)
    - [Sample Sheet](#sample-sheet)
    - [TF Barcodes](#tf-barcodes)
4. [Usage](#usage)
    - [Run](#run)
    - [Options](#options)
5. [Outputs](#outputs)
    - [Structure](#structure)
    - [File Description](#file-description)
6. [Note](#note)
    - [Junction Classification](#junction-classification)
</details>

<!-- Description-->
## Description
Synthetic Lineage Project -- by Connor Rogerson

To learn sequence to expression rules, we need to figure out how TFs interact with DNA and affect gene expression.

The pipeline is designed for single-cell multiome analysis, including scRNA-seq, scATAC-seq, and transcription factor profiling.

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- Dependencies -->
### Dependencies
<details>
<summary><b>Python Libraries</b></summary>

    pandas = 2.2.2
    numpy = 2.1.2
    biopython = 1.85
    polars = 1.25.2
</details>

<details>
<summary><b>R Packages</b></summary>

    optparse = 1.7.4
    tidyverse = 2.0.0
    data.table = 1.15.4
    ggplot2 = 3.5.2
    ggrepel = 0.9.5
    corrplot = 0.92
    GenomicRanges = 1.54.1
    Seurat = 5.3.0
    Signac = 1.12.0
    SeuratWrappers = 0.3.4
    SoupX = 1.6.2
    DoubletFinder = 2.0.6
    GenomeInfoDb = 1.38.8
    scDblFinder = 1.16.0
</details>

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- File Format -->
## File Format

<!-- Sample Sheet -->
### Sample Sheet

| sample_id | rep_id | dir_cellranger_arc | r1_tf_barcodes | r2_tf_barcodes | tf_barcodes |
|-|-|-|-|-|-|
| morf10 | B1 | /path/of/cellranger_arc/ | /path/of/R1.fastq.gz | /path/of/R2.fastq.gz | /path/of/tf_barcodes.tsv |
| morf10 | B2 | /path/of/cellranger_arc/ | /path/of/R1.fastq.gz | /path/of/R2.fastq.gz | /path/of/tf_barcodes.tsv |

> [!IMPORTANT]
> 1. The header must be like above in the example

<!-- TF Barcodes -->
### TF Barcodes

| barcode | tf_name |
|-|-|
| AGTCAAGACCCTCGGGCTCTGTGG | HIF3A-1 |
| CAATTACACCACGTCTGCCTACTA | HIF3A-2 |
| TTCAGACGTTTCGCGCCTGGAGCT | HIF3A-3 |

> [!IMPORTANT]
> 1. The header must be like above in the example
> 2. It includes TF names and barcodes

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- Usage -->
## Usage

<!-- Run -->
### Run
Please submit the bash script below

```bash
#!/bin/bash
#BSUB -o %J.o
#BSUB -e %J.e
#BSUB -R "select[mem>1000] rusage[mem=1000]"
#BSUB -M 1000
#BSUB -q long
#BSUB -J tf_pipeline

# modules
module load HGI/common/nextflow/24.10.4
module load HGI/softpack/users/fs18/sc_analysis

PIPELINE=/path/of/pipeline/dir/

SAMPLE_SHEET=/path/of/sample_sheet.tsv
OUTPUT_DIR=/path/of/output/dir
                
nextflow run -resume $PIPELINE/main.nf --sanger_module true \
                                       --sample_sheet  $SAMPLE_SHEET \
                                       --outdir        $OUTPUT_DIR
```

<!-- Options -->
### Options
#### Mandatory arguments
    --sample_sheet                path of the sample sheet
    --outdir                      the directory path of output results, default: the current directory

#### Optional arguments
    SC data QC:
    --n_rna_count                 minimum number of RNA counts, default: 1000
    --n_gene_feature              minimum number of gene features, default: 200
    --pct_mito                    maximum percentage of mitochondrial genes, default: 10
    --n_atac_count                minimum number of ATAC counts, default: 1000
    --tss_enrichment              minimum TSS enrichment score, default: 2
    --del_ambient                 whether to remove ambient RNA, default: false
    --mark_doublet                whether to mark doublets, default: false

    TF barcode QC:
    --has_umi                     whether the TF barcode has UMI, default: false
    --tf_barcode_len              length of TF barcode, default: 24
    --marker_seq                  marker sequence for TF barcode, default: "GAAAGGACGA"
    --marker_start                start position of marker sequence, default: 25
    --marker_end                  end position of marker sequence, default: 50
    --max_mismatch                maximum mismatch allowed for match sequence, default: 1
    --top_n                       keep top N TFs if 0 keep all, default: 10
    --tf_cutoff                   minimum cutoff for TF barcode filtering, default: 0
                                  if set to 0, kmeans clustering will be performed to determine the cutoff


<p align="right">(<a href="#top">back to top</a>)</p>

<!-- Outputs -->
### Outputs
### Structure

```bash
📁 output_directory
    ├─── 📁 qc_outputs
    │       ├─── 📁 morf10_B1_py_inputs
    │       │       ├─── 📁 atac
    │       │       │       ├─── 📄 barcodes.tsv.gz
    │       │       │       ├─── 📄 features.tsv.gz
    │       │       │       └─── 📄 matrix.mtx.gz
    │       │       ├─── 📁 rna
    │       │       └─── 📁 tf  
    │       ├─── 📄 morf10_B1_qc_obj_tf.rds
    │       └─── 📄 morf10_B1_qc_doublet_status.tsv
    ├─── 📁 qc_report
    │       └─── 📄 morf10_B1_qc_summary.html
    ├─── 📁 qc_sc
    │       └─── 📁 morf10_B1
    │               ├─── 📄 rds
    │               ├─── 📄 tsv
    │               └─── 📄 png
    └─── 📁 qc_tf
            └─── 📁 morf10_B1
                    ├─── 📄 tsv
                    └─── 📄 png

```

### File Description

<p align="right">(<a href="#top">back to top</a>)</p>