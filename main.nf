#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// -----------------------------
// 1. Define Pipeline Parameters
// -----------------------------

params.reads = params.reads ?: "$PWD/data/reads/*{1,2}.fastq.gz" // Location of input reads
params.outdir = params.outdir ?: "$PWD/results" // Output directory
params.genome = params.genome ?: "$PWD/data/refGenome/mm10.fa" // Path to the mouse genome FASTA file
params.gtf = params.gtf ?: "$PWD/data/refGenome/mm10.gtf" // Path to the GTF file (optional for ChIP-seq)
params.blacklist_url = params.blacklist_url ?: "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz" // URL to mouse blacklist file
params.blacklist_path = params.blacklist_path ?: "$PWD/resources/blacklist/mm10-blacklist.v2.bed" // Local path to blacklist file
params.genome_size = params.genome_size ?: "mm" // Genome size for MACS2 (mm for mouse)
params.bowtie2_threads = params.bowtie2_threads ?: 4 // Number of threads for Bowtie2
params.bowtie2_index = params.bowtie2_index ?: "$PWD/results/bowtie2_index" // Bowtie2 index directory
params.skip_alignment  = params.skip_alignment ?: false
params.test_mode       = params.test_mode ?: false

// -----------------------------
// 2. Validate Parameters
// -----------------------------

if (!params.containsKey('reads')) {
    error "Parameter 'reads' is required"
}

if (!file(params.genome).exists()) {
    error "Genome file not found: ${params.genome}"
}
if (!file(params.gtf).exists()) {
    log.warn "GTF file not found: ${params.gtf}. GTF is optional for ChIP-seq."
}


// Log Pipeline Information
log.info """\
ChIP-Seq Analysis Pipeline - FastQC, TrimGalore, Bowtie2Index, Bowtie2Align, Peak Calling (MACS2), Blacklist Filtering, Peak Annotation (ChIPseeker), and MultiQC
Reads: ${params.reads}
Output directory: ${params.outdir}
Genome: ${params.genome}
GTF: ${params.gtf}
Blacklist URL: ${params.blacklist_url}
Blacklist Path: ${params.blacklist_path}
Genome Size: ${params.genome_size}
Bowtie2 Threads: ${params.bowtie2_threads}
Bowtie2 Index: ${params.bowtie2_index}
Test mode : ${params.test_mode}
"""

// Create output directories if they do not exist
if (!file(params.outdir).exists()) {
    file(params.outdir).mkdirs()
}

// -----------------------------
// 3. Define Processes
// -----------------------------

// 3.1 Download Blacklist Regions
process DownloadBlacklist {
    tag "DownloadBlacklist"
    container 'biocontainers/wget:1.21.3--hda2e75f_0'
    publishDir "${params.outdir}/blacklist", mode: 'copy'

    output:
    path "mm10-blacklist.v2.bed" into blacklist_channel

    script:
    """
    mkdir -p blacklist
    wget -O blacklist/mm10-blacklist.v2.bed.gz ${params.blacklist_url}
    gunzip blacklist/mm10-blacklist.v2.bed.gz
    mv blacklist/mm10-blacklist.v2.bed blacklist/
    """
}

// 3.2 FastQC Process for Quality Control of Reads
process FastQC {
    tag "FastQC on ${sample_id}"
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_fastqc.zip"), path("*_fastqc.html"), emit: reports

    script:
    """
    mkdir -p qc
    fastqc -o qc ${reads[0]} ${reads[1]}
    echo "FASTQC completed for ${sample_id}"
    """
}

// 3.3 TrimGalore Process for Trimming Low-Quality Reads
process TrimGalore {
    tag "TrimGalore on ${sample_id}"
    container 'quay.io/biocontainers/trim-galore:0.6.6--hdfd78af_1'
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_1_trimmed.fq.gz"), path("${sample_id}_2_trimmed.fq.gz"), emit: trimmed_reads

    script:
    """
    mkdir -p trimmed
    trim_galore --paired --quality 20 --length 36 --gzip --output_dir trimmed ${reads[0]} ${reads[1]}
    mv trimmed/${sample_id}_1_val_1.fq.gz trimmed/${sample_id}_1_trimmed.fq.gz
    mv trimmed/${sample_id}_2_val_2.fq.gz trimmed/${sample_id}_2_trimmed.fq.gz
    echo "Trimmed Reads completed for ${sample_id}"
    """
}

// 3.4 Creating Bowtie2 Index for Alignment
process Bowtie2Index {
    tag "Bowtie2Index"
    container 'biocontainers/bowtie2:2.4.5--hdfd78af_0'
    publishDir "${params.outdir}/bowtie2_index", mode: 'copy'

    input:
    path genome_fasta

    output:
    path("bowtie2_index/*"), emit: bowtie2_index

    script:
    """
    mkdir -p bowtie2_index
    bowtie2-build ${genome_fasta} bowtie2_index/genome
    """
}

// 3.5 Bowtie2 Alignment Process for Creating Aligned BAM Files
process Bowtie2Align {
    tag "Bowtie2Align on ${sample_id}"
    container 'biocontainers/bowtie2:2.4.5--hdfd78af_0'
    publishDir "${params.outdir}/aligned", mode: 'copy'

    input:
    tuple val(sample_id), path(trimmed_reads), path(bowtie2_index)

    output:
    tuple val(sample_id), path("aligned/${sample_id}.bam"), emit: aligned_bam
    tuple val(sample_id), path("aligned/${sample_id}.bam.bai"), emit: bam_index

    script:
    """
    mkdir -p aligned
    bowtie2 -p ${params.bowtie2_threads} -x bowtie2_index/genome -1 ${trimmed_reads[0]} -2 ${trimmed_reads[1]} \
        | samtools view -bS - > aligned/${sample_id}.bam
    samtools sort aligned/${sample_id}.bam -o aligned/${sample_id}.sorted.bam
    samtools index aligned/${sample_id}.sorted.bam
    mv aligned/${sample_id}.sorted.bam aligned/${sample_id}.bam
    rm aligned/${sample_id}.sorted.bam
    echo "Bowtie2 alignment completed for ${sample_id}"
    """
}

// 3.6 Post-Alignment Processing: Sorting and Removing Duplicates
process PostAlignChIPseq {
    tag "PostAlignChIPseq on ${sample_id}"
    container 'biocontainers/samtools:1.15--hdfd78af_0'
    publishDir "${params.outdir}/aligned", mode: 'copy'

    input:
    tuple val(sample_id), path(bam) from aligned_bam

    output:
    tuple val(sample_id), path("aligned/${sample_id}.dedup.bam") emit: dedup_bam

    script:
    """
    samtools sort -o aligned/${sample_id}.sorted.bam ${bam}
    samtools index aligned/${sample_id}.sorted.bam
    picard MarkDuplicates I=aligned/${sample_id}.sorted.bam O=aligned/${sample_id}.dedup.bam M=aligned/${sample_id}.dedup_metrics.txt REMOVE_DUPLICATES=true
    samtools index aligned/${sample_id}.dedup.bam
    rm aligned/${sample_id}.sorted.bam aligned/${sample_id}.sorted.bam.bai
    echo "Post-alignment processing completed for ${sample_id}"
    """
}

// 3.7 Peak Calling with MACS2
process PeakCalling {
    tag "MACS2 Peak Calling on ${sample_id}"
    container 'biocontainers/macs2:2.2.7.1--pyhdfd78af_0'
    publishDir "${params.outdir}/peaks", mode: 'copy'

    input:
    tuple val(sample_id), path(bam) from dedup_bam

    // If you have control samples, include them here
    // For simplicity, we'll assume no control samples

    output:
    path("peaks/${sample_id}_peaks.narrowPeak") emit: peaks

    script:
    """
    mkdir -p peaks
    macs2 callpeak -t ${bam} -f BAM -g ${params.genome_size} -n ${sample_id} --outdir peaks --keep-dup all
    mv peaks/${sample_id}_peaks.narrowPeak peaks/
    echo "MACS2 peak calling completed for ${sample_id}"
    """
}

// 3.8 Blacklist Filtering of Peaks
process FilterBlacklistedPeaks {
    tag "Blacklist Filtering for ${sample_id}"
    container 'biocontainers/bedtools:2.30.0--hdfd78af_0'
    publishDir "${params.outdir}/filtered_peaks", mode: 'copy'

    input:
    tuple val(sample_id), path(peak_files) from peaks
    path blacklist from blacklist_channel

    output:
    path("filtered_peaks/${sample_id}_filtered.narrowPeak") emit: filtered_peaks

    script:
    """
    mkdir -p filtered_peaks
    bedtools intersect -v -a ${peak_files} -b ${blacklist} > filtered_peaks/${sample_id}_filtered.narrowPeak
    echo "Blacklist filtering completed for ${sample_id}"
    """
}

// 3.9 Peak Annotation with ChIPseeker
process AnnotatePeaks {
    tag "Peak Annotation for ${sample_id}"
    container 'biocontainers/r-base:4.1.0--hdfd78af_0'
    publishDir "${params.outdir}/annotated_peaks", mode: 'copy'

    input:
    tuple val(sample_id), path(filtered_peak) from filtered_peaks

    output:
    path("annotated_peaks/${sample_id}_annotated.txt") emit: annotated_peaks

    script:
    """
    Rscript annotate_peaks.R ${filtered_peak} annotated_peaks/
    echo "Peak annotation completed for ${sample_id}"
    """
}

// 3.10 MultiQC for Aggregating Reports
process MultiQC {
    tag "MultiQC"
    container 'quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path qc_reports from reports
    path "peaks/*_peaks.narrowPeak" from peaks
    path "filtered_peaks/*_filtered.narrowPeak" from filtered_peaks
    path "annotated_peaks/*.annotated.txt" from annotated_peaks

    output:
    path("multiqc_report.html") emit: multiqc_report

    script:
    """
    mkdir -p multiqc
    multiqc . -o multiqc
    echo "MultiQC report generated"
    """
}

// 3.11 Test Process for Validation (Optional)
process TestProcess {
    tag 'TestProcess'
    publishDir "${params.outdir}/test_output", mode: 'copy'

    script:
    """
    set -e
    echo "This is a test process to verify pipeline setup." > test_output.txt
    """
}

// -----------------------------
// 4. Define Workflow
// -----------------------------

workflow {
    // Define channels
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true, size: 2)
    genome_ch = Channel.value(file(params.genome))
    gtf_ch = Channel.value(file(params.gtf))
    blacklist_ch = DownloadBlacklist()

    // Bowtie2 Indexing
    bowtie2_index_ch = Bowtie2Index(genome_ch)

    // FastQC
    fastqc_ch = reads_ch.map { sample_id, files -> tuple(sample_id, files) }
    FastQC(fastqc_ch)

    // Trim Galore
    trimmed_ch = TrimGalore(reads_ch)

    // Bowtie2 Alignment
    bowtie2_align_input = trimmed_ch.combine(bowtie2_index_ch).map { sample_id, files, bowtie2_index -> tuple(sample_id, files, bowtie2_index) }
    aligned_output = Bowtie2Align(bowtie2_align_input)

    // Post-Alignment Processing
    dedup_bam = PostAlignChIPseq(aligned_output)

    // Peak Calling
    peaks = PeakCalling(dedup_bam)

    // Blacklist Filtering
    filtered_peaks = FilterBlacklistedPeaks(peaks, blacklist_ch)

    // Peak Annotation
    // annotated_peaks = AnnotatePeaks(filtered_peaks)

    // // MultiQC
    // MultiQC(fastqc_ch, peaks, filtered_peaks, annotated_peaks)

    // Optional Test Process
    if (params.test_mode) {
        TestProcess()
    }
}
