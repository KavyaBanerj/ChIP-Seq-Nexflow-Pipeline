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
params.keep_dup = params.keep_dup ?: "auto" // MACS2 keep-dup parameter
params.skip_alignment = params.skip_alignment ?: false
params.test_mode = params.test_mode ?: false
params.read_type = params.read_type ?: 'paired' // 'paired' or 'single'

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
log.info """
ChIP-Seq Analysis Pipeline - FastQC, TrimGalore, Bowtie2Index, Bowtie2Align, Peak Calling (MACS2), Blacklist Filtering, Peak Annotation (ChIPseeker), and MultiQC
Reads: ${params.reads}
Read type: ${params.read_type}
Output directory: ${params.outdir}
Genome: ${params.genome}
GTF: ${params.gtf}
Blacklist URL: ${params.blacklist_url}
Blacklist Path: ${params.blacklist_path}
Genome Size: ${params.genome_size}
Bowtie2 Threads: ${params.bowtie2_threads}
Bowtie2 Index: ${params.bowtie2_index}
Keep Duplication: ${params.keep_dup}
Skip Alignment: ${params.skip_alignment}
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
    publishDir "${params.outdir}/blacklist", mode: 'copy'

    output:
    path "mm10-blacklist.v2.bed", emit: blacklist_channel

    script:
    """
    wget -O mm10-blacklist.v2.bed.gz ${params.blacklist_url}
    gunzip mm10-blacklist.v2.bed.gz
    """
}

// 3.2 FastQC Process for Quality Control of Reads
process FastQC {
    tag "FastQC on ${sample_id}"
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.html", emit: reports_html
    path "*.zip", emit: reports_zip

    script:
    """
    fastqc -t ${task.cpus} -o . ${reads.join(' ')}
    echo "FASTQC completed for ${sample_id}"
    """
}

// 3.3 TrimGalore Process for Trimming Low-Quality Reads
process TrimGalore {
    tag "TrimGalore on ${sample_id}"
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_trimmed.fq.gz"), emit: trimmed_reads

    script:
    if (params.read_type == 'paired') {
        """
        trim_galore --paired --cores ${task.cpus} --quality 20 --length 36 --gzip --output_dir . ${reads[0]} ${reads[1]}
        mv *_val_1.fq.gz ${sample_id}_1_trimmed.fq.gz
        mv *_val_2.fq.gz ${sample_id}_2_trimmed.fq.gz
        echo "Trimmed paired-end reads for ${sample_id}"
        """
    } else {
        """
        trim_galore --cores ${task.cpus} --quality 20 --length 36 --gzip --output_dir . ${reads[0]}
        echo "Trimmed single-end reads for ${sample_id}"
        """
    
    }
}  

// 3.4 Bowtie2 Index Process
process Bowtie2Index {
    tag "Bowtie2Index"
    publishDir "${params.outdir}/bowtie2_index", mode: 'copy'

    input:
    path genome_fasta

    output:
    path("genome.*"), emit: bowtie2_index

    script:
    """
    bowtie2-build --threads ${task.cpus} ${genome_fasta} genome
    echo "Genome build completed for ${genome_fasta}"
    """
}

// 3.5 Bowtie2 Alignment Process
process Bowtie2Align {
    tag "Bowtie2Align on ${sample_id}"
    publishDir "${params.outdir}/aligned", mode: 'copy'

    input:
    tuple val(sample_id), path(trimmed_reads), val(bowtie2_index_prefix)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam_output

    script:
    if (params.read_type == 'paired') {
        """
        bowtie2 -p ${params.bowtie2_threads} -x ${bowtie2_index_prefix} \
                -1 ${trimmed_reads[0]} -2 ${trimmed_reads[1]} \
            | samtools view -bS - > ${sample_id}.bam
        echo "Bowtie2 alignment completed for ${sample_id}"
        """
    } else {
        """
        bowtie2 -p ${params.bowtie2_threads} -x ${bowtie2_index_prefix} \
                -U ${trimmed_reads[0]} \
            | samtools view -bS - > ${sample_id}.bam
        echo "Bowtie2 alignment completed for ${sample_id}"
        """
    }
}

// 3.6 Sort & Index Bam Process
process SortIndexBam {
    tag "SortIndexBam on ${sample_id}"
    publishDir "${params.outdir}/aligned", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_output)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: sorted_bam

    script:
    """
    samtools sort -@ ${task.cpus} -m 1G -l 9 -o ${sample_id}.sorted.bam ${bam_output}
    samtools index ${sample_id}.sorted.bam
    echo "Samtools sorting and indexing completed for ${sample_id}"
    """
}


// 3.7 Post-Alignment Processing: Removing Duplicates
process PostAlignChIPseq {
    tag "PostAlignChIPseq on ${sample_id}"
    publishDir "${params.outdir}/aligned", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam"), emit: dedup_bam

    script:
    """
    picard MarkDuplicates I=${bam} O=${sample_id}.dedup.bam M=${sample_id}.dedup_metrics.txt REMOVE_DUPLICATES=true
    echo "Post-alignment processing completed for ${sample_id}"
    """
}


// 3.8 Alignment QC Metrics & BAM Indexing
process AlignmentQCAndIndex {
    tag "Alignment QC & BAM Indexing on ${sample_id} dedup reads"
    publishDir "${params.outdir}/aligned", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    path("${sample_id}_flagstat.txt"), emit: qc_metrics
    path("${sample_id}.dedup.bam.bai"), emit: bam_index

    script:
    """
    samtools index ${bam}
    samtools flagstat ${bam} > ${sample_id}_flagstat.txt
    echo "BAM indexing and QC metrics generated for ${sample_id}"
    """
}

// 3.9 Peak Calling with MACS2
// to do: figure out the control/spike-in sample command
process PeakCalling {
    tag "MACS2 Peak Calling on ${sample_id}"
    publishDir "${params.outdir}/peaks", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_peaks.narrowPeak"), emit: peaks
    tuple val(sample_id), path("${sample_id}_peaks.xls")
    tuple val(sample_id), path("${sample_id}_summits.bed")

    script:
    """
    macs2 callpeak -t ${bam} -f BAM -g ${params.genome_size} -n ${sample_id} --outdir . --keep-dup ${params.keep_dup} --nolambda --nomodel -q 0.05
    echo "MACS2 peak calling completed for ${sample_id}"
    """
}

// 3.10 Blacklist Filtering of Peaks
process FilterBlacklistedPeaks {
    tag "Blacklist Filtering for ${sample_id}"
    publishDir "${params.outdir}/filtered_peaks", mode: 'copy'

    input:
    tuple val(sample_id), path(peak_files)
    path blacklist

    output:
    tuple val(sample_id), path("${sample_id}_filtered.narrowPeak"), emit: filtered_peaks

    script:
    """
    bedtools intersect -v -a ${peak_files} -b ${blacklist} > ${sample_id}_filtered.narrowPeak
    echo "Blacklist filtering completed for ${sample_id}"
    """
}

// 3.11 Peak Annotation with HOMER
process AnnotatePeaks {
    tag "Annotate Peaks with HOMER for ${sample_id}"
    publishDir "${params.outdir}/annotated_peaks", mode: 'copy'

    input:
    tuple val(sample_id), path(peak_files)
    path genome_fasta

    output:
    path("${sample_id}_annotated_peaks.txt"), emit: homer_annotated

    script:
    """
    annotatePeaks.pl ${peak_files} ${genome_fasta} > ${sample_id}_annotated_peaks.txt
    echo "HOMER annotation completed for ${sample_id}"
    """
}

// 3.13  Custom Genome reference for Motif Analysis with HOMER - to do, requires a custom container with homer and samtools
// process CreateHomerGenome {
//     tag "CreateHomerGenome for ${params.genome}"
//     publishDir "${params.outdir}/homer_genomes", mode: 'copy'

//     input:
//     path genome_fasta

//     output:
//     path "homer_genome", emit: homer_genome_dir

//     script:
//     """
//     mkdir -p homer_genome
//     cp ${genome_fasta} homer_genome/genome.fa
//     samtools faidx homer_genome/genome.fa
//     cut -f1,2 homer_genome/genome.fa.fai > homer_genome/chrom.sizes

//     perl /usr/local/share/homer/configureHomer.pl -install homer_genome
//     echo "HOMER genome directory created for ${genome_fasta}"
//     """
// }

// 3.13  Motif Analysis with HOMER
process MotifAnalysis {
    tag "MotifAnalysis on ${sample_id}"
    publishDir "${params.outdir}/motifs", mode: 'copy'

    input:
    tuple val(sample_id), path(peaks)

    output:
    path "${sample_id}_motifs/knownResults.html", emit: motif_results

    script:
    """
    mkdir -p ${sample_id}_motifs
    perl /usr/local/share/homer/.//configureHomer.pl -install mm10
    findMotifsGenome.pl ${peaks} mm10 ${sample_id}_motifs -size 200 -mask
    echo "Motif analysis completed for ${sample_id}"
    """
}


// 3.13 Bam Coverage (BigWig file generation) with deepTools
process BamCoverage {
    tag "bamCoverage on ${sample_id}"
    publishDir "${params.outdir}/bigwig", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.bw"), emit: bigwig

    script:
    """
    bamCoverage -b ${bam} -o ${sample_id}.bw --normalizeUsing CPM --numberOfProcessors ${task.cpus}
    echo "BigWig file generated for ${sample_id}"
    """
}


// 3.14 Multi Bigwig Summary with deepTools
process MultiBigwigSummary {
    tag "multiBigwigSummary"
    publishDir "${params.outdir}/correlation", mode: 'copy'

    input:
    path bigwigs from bigwig.collect()

    output:
    path "correlation_matrix.npz", emit: summary_matrix

    script:
    """
    multiBigwigSummary bins -b ${bigwigs.join(' ')} -o correlation_matrix.npz --numberOfProcessors ${task.cpus}
    """
}

// 3.15 Correlation Analysis with deepTools
process PlotCorrelation {
    tag "plotCorrelation"
    publishDir "${params.outdir}/correlation", mode: 'copy'

    input:
    path summary_matrix from MultiBigwigSummary.summary_matrix

    output:
    path "correlation_heatmap.png", emit: correlation_plot
    path "correlation_matrix.tsv", emit: correlation_matrix

    script:
    """
    plotCorrelation --corData ${summary_matrix} --corMethod pearson --skipZeros \
                    --whatToPlot heatmap --plotFile correlation_heatmap.png \
                    --outFileCorMatrix correlation_matrix.tsv
    """
}


// 3.16 MultiQC for Aggregating Reports
process MultiQC {
    tag "MultiQC"
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path qc_reports from reports
    path "peaks/*_peaks.narrowPeak" from peaks
    path "filtered_peaks/*_filtered.narrowPeak" from filtered_peaks
    path "annotated_peaks/*.annotated.txt" from annotated_peaks

    output:
    path("multiqc_report.html"), emit: multiqc_report

    script:
    """
    multiqc . -o .
    echo "MultiQC report generated"
    """
}

// 3.17 Test Process for Validation (Optional)
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
    // Define channels for input reads and genome
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true, size: params.read_type == 'paired' ? 2 : 1)
    genome_ch = Channel.value(file(params.genome))
    gtf_ch = Channel.value(file(params.gtf))
    blacklist_ch = DownloadBlacklist()

    // Step 1: Bowtie2 Indexing
    bowtie2_index_ch = Bowtie2Index(genome_ch).map { index_paths -> index_paths[0].parent + '/genome' }

    // Step 2: FastQC
    fastqc_ch = reads_ch.map { sample_id, files -> tuple(sample_id, files) }
    FastQC(fastqc_ch)

    // Step 3: Trim Galore
    trimmed_ch = TrimGalore(reads_ch)

    // Step 4: Bowtie2 Alignment
    bowtie2_align_input = trimmed_ch.combine(bowtie2_index_ch)
    bam_output = Bowtie2Align(bowtie2_align_input)

    // Step 5: Sorting and Indexing BAM
    sorted_bam = SortIndexBam(bam_output)

    // Step 6: Post-Alignment Processing (Removing Duplicates)
    dedup_bam = PostAlignChIPseq(sorted_bam)

    // Step 7: Alignment QC and BAM Indexing
    qc_metrics, bam_index = AlignmentQCAndIndex(dedup_bam)

    // Step 8: Peak Calling with MACS2
    peaks_ch = PeakCalling(dedup_bam)

    // Step 9: Blacklist Filtering
    filtered_peaks = FilterBlacklistedPeaks(peaks_ch.peaks, blacklist_ch)

    // Step 10: Annotate peaks with HOMER
    annotated_peaks = AnnotatePeaks(filtered_peaks, genome_ch)

    // Step 11: Motif Analysis with HOMER
    motif_results = MotifAnalysis(filtered_peaks)

    // Step 12: BigWig Generation
    bigwig = BamCoverage(bam_index)

    // Optional Correlation Analysis
    // summary_matrix = MultiBigwigSummary(bigwig.collect())
    // correlation_plot = PlotCorrelation(summary_matrix)

    // Optional MultiQC
    // MultiQC(fastqc_ch, peaks, filtered_peaks, annotated_peaks)

    // Optional Test Process
    if (params.test_mode) {
        TestProcess()
    }
}