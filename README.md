# ChIP-Seq Analysis Pipeline

This pipeline is designed to perform comprehensive ChIP-Seq analysis, including quality control, alignment, peak calling, blacklist filtering, annotation, motif analysis, and visualization. The pipeline is implemented using Nextflow DSL2 and supports paired-end and single-end sequencing data.

## Pipeline Overview

This ChIP-Seq analysis pipeline performs the following steps:

1. **Download Blacklist Regions**: Downloads the genome-specific blacklist file.
2. **Quality Control**: Performs read quality control using FastQC.
3. **Read Trimming**: Uses TrimGalore to trim low-quality reads and adapters.
4. **Alignment**: Aligns reads to the reference genome using Bowtie2.
5. **Post-Alignment Processing**:
   - Sorting and indexing of BAM files.
   - Removing duplicate reads using Picard.
   - Generating alignment QC metrics.
6. **Peak Calling**: Calls peaks using MACS2.
7. **Blacklist Filtering**: Removes peaks that overlap with blacklist regions.
8. **Peak Annotation**: Annotates peaks using HOMER.
9. **Motif Analysis**: Identifies enriched motifs using HOMER.
10. **Coverage Analysis**:
    - Generates BigWig files for visualization.
    - Computes coverage matrices and plots coverage profiles.
11. **Correlation Analysis**: Generates correlation matrices and plots.

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/KavyaBanerj/ChIP-Seq-Nexflow-Pipeline/tree/main
   cd <repository-directory>
   ```
2. Install Nextflow:
   ```bash
   curl -s https://get.nextflow.io | bash
   mv nextflow ~/bin/
   ```
3. Ensure the necessary software and tools are installed (see Requirements).

## Requirements

The following tools are required for running the pipeline:

- Nextflow
- FastQC
- TrimGalore
- Bowtie2
- SAMtools
- Picard
- MACS2
- HOMER
- bedtools
- deepTools
- MultiQC
- wget

Ensure these tools are available in your system PATH or in the container used for the pipeline.

## Parameters

| Parameter           | Description                                               | Default Value                             |
|---------------------|-----------------------------------------------------------|-------------------------------------------|
| `reads`             | Location of input reads (supports glob patterns).          | `"$PWD/data/reads/*{1,2}.fastq.gz"`      |
| `outdir`            | Output directory.                                         | `"$PWD/results"`                         |
| `genome`            | Path to the genome FASTA file.                            | `"$PWD/data/refGenome/mm10.fa"`          |
| `gtf`               | Path to the GTF file (optional for ChIP-Seq).             | `"$PWD/data/refGenome/mm10.gtf"`         |
| `blacklist_url`     | URL to download blacklist regions.                        | `"https://raw.githubusercontent.com/..."`|
| `blacklist_path`    | Local path to the blacklist file.                         | `"$PWD/resources/blacklist/mm10-blacklist.v2.bed"` |
| `genome_size`       | Genome size for MACS2 peak calling.                       | `"mm"`                                   |
| `bowtie2_threads`   | Number of threads for Bowtie2 alignment.                  | `4`                                       |
| `bowtie2_index`     | Bowtie2 index directory.                                  | `"$PWD/results/bowtie2_index"`           |
| `keep_dup`          | MACS2 keep-dup parameter.                                 | `"auto"`                                  |
| `skip_alignment`    | Skip alignment step if set to true.                       | `false`                                   |
| `test_mode`         | Run a test process to verify setup.                       | `false`                                   |
| `read_type`         | Specify read type: `paired` or `single`.                  | `"paired"`                                |

## Pipeline Workflow

The pipeline is divided into several processes, each handling a specific task:

### 1. Download Blacklist Regions
Downloads the genome-specific blacklist file from the specified URL.

### 2. Quality Control (FastQC)
Runs FastQC to generate quality control reports for the input reads.

### 3. Read Trimming (TrimGalore)
Trims low-quality bases and adapters from the reads using TrimGalore.

### 4. Alignment (Bowtie2)
Aligns the trimmed reads to the reference genome using Bowtie2 and converts the output to BAM format using SAMtools.

### 5. Post-Alignment Processing
- **Sorting and Indexing**: Sorts and indexes the aligned BAM files.
- **Duplicate Removal**: Removes duplicate reads using Picard.
- **Alignment QC**: Generates alignment statistics and indices.

### 6. Peak Calling (MACS2)
Identifies enriched regions (peaks) using MACS2 with the specified genome size.

### 7. Blacklist Filtering (bedtools)
Filters out peaks overlapping with blacklist regions.

### 8. Peak Annotation (HOMER)
Annotates peaks using HOMER, providing information on genomic features.

### 9. Motif Analysis (HOMER)
Performs motif analysis to identify enriched motifs in the peak regions.

### 10. Coverage Analysis (deepTools)
- **BigWig Generation**: Creates BigWig files for visualization.
- **Compute Matrix**: Computes coverage matrices.
- **Plot Coverage Profile**: Generates coverage profile plots.

### 11. Correlation Analysis (deepTools)
Generates correlation matrices and plots based on the coverage data.

## Output Structure

The pipeline generates the following output directories:

```
results/
├── aligned/            # Aligned BAM files and indices
├── annotated_peaks/    # Annotated peak files
├── blacklist/          # Blacklist regions
├── bigwig/             # BigWig files for coverage visualization
├── correlation/        # Correlation matrices and plots
├── filtered_peaks/     # Blacklist-filtered peak files
├── matrix/             # Coverage matrices
├── motifs/             # HOMER motif analysis results
├── peaks/              # MACS2 peak calling results
├── qc/                 # FastQC quality control reports
├── trimmed/            # Trimmed reads
```

## Usage

Run the pipeline using the following command:

```bash
nextflow run main.nf --reads "path/to/reads/*{1,2}.fastq.gz" --genome "path/to/genome.fa" --outdir "path/to/output"
```

To run in test mode:

```bash
nextflow run main.nf --test_mode true
```

## Customizing the Pipeline

You can customize the pipeline by modifying the parameters in the `main.nf` file or by specifying them at runtime using the `--` prefix.

Example:

```bash
nextflow run main.nf --reads "data/*.fastq.gz"  --read_type paired
```

## Test Mode

The pipeline includes a test process that can be run independently to validate the setup and ensure that the required tools are available.

## Dependencies

The pipeline uses the following software containers from [Biocontainers](https://biocontainers.pro/). Please refer to the nextflow.config file for the Docker containers used. Ensure these containers are available in your Docker environment.


## **Future Scope**

This ChIP-Seq analysis pipeline is designed to be modular, scalable, and adaptable. Future enhancements and extensions planned for the pipeline include:

1. **Spike-in Normalization**  
   Add support for spike-in controls to normalize ChIP-Seq signals and ensure accurate comparison across samples.

2. **Cloud Deployment**  
   Improve cloud compatibility by creating profiles for AWS to facilitate large-scale data processing and reproducibility in cloud environments.



