#!/usr/bin/env Rscript

# Which SNPs overlap ATAC-seq peaks?

library("GenomicRanges")
library("readr")
library("rtracklayer")
library("tidyr")

# Define genome annotation
genome <- Seqinfo(genome = "hg19")

# Format SNPs
bim <- read_tsv("../data/plink/cardioqtl.bim", col_names = FALSE)
colnames(bim) <- c("chr", "rs", "X3", "ps", "allele1", "allele0")
bim$chr <- paste0("chr", bim$chr)
gr_bim <- makeGRangesFromDataFrame(bim,
                                   keep.extra.columns = TRUE,
                                   seqinfo = genome,
                                   start.field = "ps", end.field = "ps")

# Format ATAC-seq peaks
f_atac <- file.path("../data/atac/",
                    c("CM_ATACpeaks.bed.gz",
                      "IPS_ATACpeaks.bed.gz",
                      "LCL_ATACpeaks.bed.gz"))
names(f_atac) <- c("cm", "ips", "lcl")
l_atac <- List(lapply(f_atac, import, genome = genome))
gr_atac <- stack(l_atac, "cell")

# Intersect the SNPs with the ATAC-seq peaks
pairs <- findOverlapPairs(gr_bim, gr_atac, ignore.strand = TRUE)

# Format into 0 or 1 for overlap or not
snps_in_peaks <- data.frame(rs = mcols(first(pairs))$rs,
                            cell = mcols(second(pairs))$cell,
                            detect = 1)
snps_in_peaks_wide <- spread(snps_in_peaks, key = cell, value = detect)
stopifnot(nrow(snps_in_peaks_wide) == length(unique(snps_in_peaks$rs)))
for (cell in c("cm", "ips", "lcl")) {
  snps_in_peaks_wide[, cell] <- ifelse(is.na(snps_in_peaks_wide[, cell]),
                                       0, snps_in_peaks_wide[, cell])
}

write_tsv(snps_in_peaks_wide, "../data/atac/intersect.txt")
