# 18S_metabarcoding_pipeline
A bioinformatics pipeline for analyzing long-read 18S rRNA metabarcoding data from Oxford Nanopore Technologies (ONT). It processes raw sequencing reads through primer trimming, 18S region selection, optimized quality-based clustering, consensus polishing, and taxonomic assignment against the [PR² database](https://pr2-database.org/).

This pipeline was developed and validated for assessing soil protist diversity, demonstrating that long-read metabarcoding with optimized bioinformatic analysis outperforms short-read metabarcoding for this application.
 
> **Status:** the article describing and applying this pipeline is currently in preparation / under review. Citation details will be added here once it is published.

## Table of Contents

- [Overview](#overview)
- [Pipeline Workflow](#pipeline-workflow)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
  - [Step 1: Read processing and individual taxonomic assignment](#step-1:-read-processing-and-individual-taxonomic-assignment)
  - [Step 2: Optimized clustering and consensus sequence taxonomic assigment](#step-2:-optimized-clustering-and-consensus-sequence-taxonomic-assignment)
- [Key Output Files](#key-output-files)
- [License](#license)
- [Contact](#contact)

## Overview
 
This pipeline takes raw ONT `fastq_pass` reads (organised by barcode) and processes them through primer trimming, 18S rRNA region extraction, quality/length filtering, per-read taxonomic assignment, sequence clustering, consensus polishing, and a final taxonomic assignment of consensus (cluster) sequences against the PR2 reference database.
 
It is composed of two consecutive bash scripts, each calling a set of custom Python scripts:
 
1. **`18S_metabarcoding_pipeline_WP2_article.sh`** per-barcode preprocessing, primer trimming, 18S extraction, and individual taxonomic assignment.
3. **`18S_metabarcoding_clustering_WP2_article.sh`** optimized clustering of the pooled 18S sequences, consensus polishing, and taxonomic assignment of consensus sequences.
   
## Pipeline workflow
 
```
Raw fastq_pass reads (per barcode)
        │
        ▼
Step 1.1 Merge & convert to FASTA, log raw read counts, BLAST primer sequences against reads (optional check)
        │
        ▼
Step 1.2 Primer trimming (cutadapt)
        │
        ▼
Step 1.3 rRNA gene detection & 18S region extraction (barrnap + custom script)
        │
        ▼
Step 1.4 Pooling 18S sequences 
        │
        ▼
Step 1.5 Quality & length filtering, Per-read taxonomic assignment vs. PR2 
        │
        ▼
Step 2.1 Sort sequences by quality, filter by length
        │
        ▼
Step 2.2 Cluster sequences (vsearch --cluster_smallmem, 98% identity)
        │
        ▼
Step 2.3 Select clusters ≥ 5 sequences, Consensus polishing per cluster (medaka)
        │
        ▼
Step 2.4 Taxonomic assignment of consensus sequences vs. PR2
        │
        ▼
Final grouped taxonomic OTU table
```

## Requirements

### Software Dependencies

| Tool | Version | Conda Environment |
|------|---------|-------------------|
| [NCBI BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) | 2.16.0+ | `blast_env` |
| [Cutadapt](https://cutadapt.readthedocs.io/) | 5.0 | `cutadapt_env` |
| [Barrnap](https://github.com/tseemann/barrnap) | 0.9 | `barrnap_env` |
| [VSEARCH](https://github.com/torognes/vsearch) | 2.29.4 | `vsearch_env` |
| [Medaka](https://github.com/nanoporetech/medaka) | 2.0.1 | `medaka_env` |
| [Biopython](https://biopython.org/) | 1.85 | (Python scripts) |

### Reference Database
This pipeline was used with the PR2 database. 
| Database | Version | File Used |
|----------|---------|-----------|
| [PR²](https://pr2database.appspot.com/) | 5.1.0 | `pr2_version_5.1.0_SSU_mothur.fasta` + `.tax` |

#### Use the pipeline with another database

To use the pipeline with another database, the script can be adapted. The easiest is to replace the names of the .fasta and .tax files with another database files in `18S_metabarcoding_pipeline/18S_metabarcoding_clustering_WP2_article.sh` lines 129 and 151 and in `18S_metabarcoding_pipeline/18S_metabarcoding_pipeline_WP2_article.sh` line 233.

### Primers

The pipeline currently uses hardcoded primer sequences:

| Primer | Sequence |
|--------|----------|
| Forward | `GGCAAGTCTGGTGCCAG` |
| Reverse variant 1 | `AAGGTAGCCAAATGCCTCGTC` |
| Reverse variant 2 | `AYTWGTGACGYGCATGAATGG` |

> ⚠️ If your study uses different primers, you will need to modify the `cutadapt` commands in the script accordingly.

## Installation
### 1. Clone the repository

```bash
git clone https://github.com/Aline-Git/18S_metabarcoding_pipeline.git
cd 18S_metabarcoding_pipeline
chmod a+x *.sh
```

### 2. Set up conda environments

Create the required conda environments for each dependency:

```bash
conda create -n blast_env -c bioconda blast=2.16.0 
conda create -n cutadapt_env -c bioconda cutadapt=5.0 
conda create -n barrnap_env -c bioconda barrnap=0.9 
conda create -n vsearch_env -c bioconda vsearch=2.29.4 
conda create -n medaka_env -c bioconda medaka=2.0.1
```

Install Biopython in the environment used to run the Python helper scripts:

```bash
pip install biopython==1.85
```

### 3. Download the PR² database

Download PR² version 5.1.0 from the [PR² database website](https://pr2database.appspot.com/) and place the Mothur-formatted files in your database directory:

```
<path_to_database_dir>/ 
                  ├── pr2_version_5.1.0_SSU_mothur.fasta 
                  └── pr2_version_5.1.0_SSU_mothur.tax
```

### 4. Make Python scripts executable:

```bash
chmod a+x ~/scripts/.py ~/scripts/extract_18S_sequences_with_barcodeID/.py
```

## Usage

### Input Data Structure

Organize your ONT basecalled output in the standard `fastq_pass` directory format:

```
workdir/ 
 └── fastq_pass/ 
       ├── barcode01/
       │   ├── batch_1.fastq.gz
       │   ├── batch_2.fastq.gz 
       │   └── ... 
       ├── barcode02/ 
       │   └── ... 
       └── ...
```

### Step 1: Read processing and individual taxonomic assignment

```bash
./18S_metabarcoding_pipeline/18S_metabarcoding_pipeline_WP2_article.sh \
  <path_to_workdir> \
  <path_to_primer_file> \
  <path_to_database_dir>
```

| Argument | Description |
|----------|-------------|
| `path_to_workdir` | Path to the working directory containing `fastq_pass/` |
| `path_to_primer_file` | Path to the FASTA file containing the sequences of the primers used to amplify the DNA sequenced |
| `path_to_database_dir` | Path to the directory containing the database files (pr2_version_5.1.0_SSU_mothur.fasta and pr2_version_5.1.0_SSU_mothur.tax) |

**Example:**

```bash
./18S_metabarcoding_pipeline/18S_metabarcoding_pipeline_WP2_article.sh \
  /data/runs/soil_protist_2024 \
  /data/primers/my_primers.fasta \
  /data/references/PR2_db_v5_1_0
```

The selected 18S sequences from each barcode folder will be concatenated in grouped file at `$WORKDIR/clustering/input_sequences/grouped_18S_selection.fastq`, preparing the input for Step 2.

### Step 2: Clustering and consensus sequence taxonomic assigment

Once all barcodes have been processed, run the clustering script:

```bash
./18S_metabarcoding_pipeline/18S_metabarcoding_clustering_WP2_article.sh \
  <path_to_workdir> \
  <clustering_dir_name> \
  <path_to_database_dir>
```

| Argument | Description |
|----------|-------------|
| `path_to_workdir` | Same working directory used in Step 1 |
| `clustering_dir_name` | Name for the clustering output directory (e.g., `clustering`) |
| `path_to_database_dir` | Path to the directory containing the database files (pr2_version_5.1.0_SSU_mothur.fasta and pr2_version_5.1.0_SSU_mothur.tax) |

**Example:**

```bash
./18S_metabarcoding_pipeline/18S_metabarcoding_clustering_WP2_article.sh \
  /data/runs/soil_protist_2024 \
  clustering
```

### Key Output Files

| File | Description |
|------|-------------|
| `logs/<barcode>_recap_count_sequences.csv` | Read counts at each processing stage sample <barcode> |
| `<barcode>_18S_selection.fastq` | Reads identified as 18S rRNA |
| `\usearch_18S_selection_pr2v51\<barcode>\<barcode>_usearch_output_query_cov_90_best_strand.txt` | Per-read taxonomic assignment  (Step 1) |
| `<clustering_dir>/polishing/centroids_consensus.fasta` | Polished consensus sequences for all clusters ≥ 5 reads |
| `<clustering_dir>/usearch_global_cov90_entroids_18S_clusters_size5/usearch_output_query_cov_90_best_strand.txt` | taxonomic assignment of consensus sequences (Step 2)|
| Final taxonomic table | Grouped OTU table with taxonomy (from `make_grouped_taxonomic_table_from_usearch_output.py`) |


## License

This project is licensed under the terms specified in the [LICENSE](LICENSE) file.

## Contact

For questions or support, please open an issue on [GitHub](https://github.com/Aline-Git/18S_metabarcoding_pipeline/issues) or contact:

📧 **aline.git@etik.com**



*This pipeline was developed for use with Oxford Nanopore Technologies long-read sequencing data.*
