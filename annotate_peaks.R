#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ChIPseeker)
    library(GenomicFeatures)
    library(TxDb.Mmusculus.UCSC.mm10.knownGene) # Mouse TxDb
    library(org.Mm.eg.db) # Mouse OrgDb
})

args <- commandArgs(trailingOnly = TRUE)
peak_file <- args[1]
output_dir <- args[2]

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene # Mouse TxDb

peak <- readPeakFile(peak_file)
peakAnno <- annotatePeak(peak, TxDb=txdb, tssRegion=c(-3000, 3000), annoDb="org.Mm.eg.db")

# Save annotation
output_file <- file.path(output_dir, paste0(basename(peak_file), ".annotated.txt"))
write.table(as.data.frame(peakAnno), file=output_file, sep="\t", quote=FALSE, row.names=FALSE)
