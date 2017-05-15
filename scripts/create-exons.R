#!/usr/bin/env Rscript

# Create exons file for for mapping reads to genes with Subread featureCounts.

# Usage:
#   Rscript create-exons.R <archive> > out.saf
#
# archive - URL of Ensembl archive, e.g. mar2017.archive.ensembl.org
#
# Notes:
# + Output is in Simplified Annotation Format (SAF)
#     + Columns: GeneID, Chr, Start, End, Strand
#     + Coordinates are 1-based, inclusive on both ends
# + Contains duplicate and overlapping exons (featureCounts handles this)

suppressMessages(library("biomaRt"))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
archive <- args[1]

ensembl <- useMart(host = archive,
                   biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "hsapiens_gene_ensembl")

exons_all <- getBM(attributes = c("ensembl_gene_id", "ensembl_exon_id",
                                  "chromosome_name", "exon_chrom_start",
                                  "exon_chrom_end", "strand",
                                  "external_gene_name",
                                  "gene_biotype"),
                   mart = ensembl)
exons_final <- exons_all[exons_all$chromosome_name %in% c(1:22, "X", "Y", "MT"),
                         c("ensembl_gene_id", "chromosome_name",
                           "exon_chrom_start", "exon_chrom_end",
                           "strand", "external_gene_name",
                           "gene_biotype")]
colnames(exons_final) <- c("GeneID", "Chr", "Start", "End", "Strand",
                           "Name", "Biotype")
# Sort by Ensembl gene ID, then start and end positins
exons_final <- exons_final[order(exons_final$GeneID,
                                 exons_final$Start,
                                 exons_final$End), ]
# Fix strand
exons_final$Strand <- ifelse(exons_final$Strand == 1, "+", "-")

# Save as tab-separated file in Simplified Annotation Format (SAF)
write.table(exons_final, "", quote = FALSE, sep = "\t",
            row.names = FALSE)
