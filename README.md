# 18S_metabarcoding_pipeline
A bioinformatics pipeline for analyzing long-read 18S rRNA metabarcoding data from Oxford Nanopore Technologies (ONT). It processes raw sequencing reads through primer trimming, 18S region selection, quality-based clustering, consensus polishing, and taxonomic assignment against the [PR² database](https://pr2-database.org/).

This pipeline was developed and validated for assessing soil protist diversity, demonstrating that long-read metabarcoding with optimized bioinformatic analysis outperforms short-read metabarcoding for this application.
 
> **Status:** the article describing and applying this pipeline is currently in preparation / under review. Citation details will be added here once it is published.

## Table of Contents

- [Overview](#overview)
- [Pipeline Workflow](#pipeline-workflow)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
  - [Step 1: Read processing and individual taxonomic assignment](#step-1:-read-processing-and-individual-taxonomic-assignment)
  - [Step 2: Clustering and consensus sequence taxonomic assigment](#step-2:-clustering-and-consensus-sequence-taxonomic-assignment)
- [Output Structure](#output-structure)
- [Helper Scripts](#helper-scripts)
- [License](#license)
- [Contact](#contact)


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

| Database | Version | File Used |
|----------|---------|-----------|
| [PR²](https://pr2database.appspot.com/) | 5.1.0 | `pr2_version_5.1.0_SSU_mothur.fasta` + `.tax` |

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
git clone https://github.com/Aline-Git/18S_metabarcoding_pipeline.git cd 18S_metabarcoding_pipeline chmod a+x *.sh
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
       │ ├── batch_1.fastq.gz
       │ ├── batch_2.fastq.gz 
       │ └── ... 
       ├── barcode02/ 
       │ └── ... 
       └── ...
```

### Step 1: Read processing and individual taxonomic assignment

```bash
./18S_metabarcoding_pipeline/18S_metabarcoding_pipeline_WP2_article.sh \
  <path_to_workdir> \
  <path_to_primer_file> \
  <path_to_database_folder>
```

| Argument | Description |
|----------|-------------|
| `path_to_workdir` | Path to the working directory containing `fastq_pass/` |
| `path_to_primer_file` | Path to the FASTA file containing the sequences of the primers used to amplify the DNA sequenced |
| `path_to_database_folder` | Path to the folder containing the database files (pr2_version_5.1.0_SSU_mothur.fasta and pr2_version_5.1.0_SSU_mothur.tax) |

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
  <clustering_dir_name>
```

| Argument | Description |
|----------|-------------|
| `path_to_workdir` | Same working directory used in Step 1 |
| `clustering_dir_name` | Name for the clustering output directory (e.g., `clustering`) |

**Example:**

```bash
./18S_metabarcoding_pipeline/18S_metabarcoding_clustering_WP2_article.sh \
  /data/runs/soil_protist_2024 \
  clustering
```

## License

This project is licensed under the terms specified in the [LICENSE](LICENSE) file.

## Contact

For questions or support, please open an issue on [GitHub](https://github.com/Aline-Git/18S_metabarcoding_pipeline/issues) or contact:

📧 **aline.git@etik.com**



*This pipeline was developed for use with Oxford Nanopore Technologies long-read sequencing data.*
