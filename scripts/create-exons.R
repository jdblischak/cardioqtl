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
# For testing:
# archive <- "mar2017.archive.ensembl.org"
# archive <- "feb2014.archive.ensembl.org"

ensembl <- useMart(host = archive,
                   biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "hsapiens_gene_ensembl")

# The name of the attribute to obtain the gene name has changed over different
# Ensembl versions (I didn't try to pinpoint the exact switch point). Thus need
# to obtain this dynamically based on the input.
#
# Ensembl 75 (GRCh37): external_gene_id "Associated Gene Name"
#
# Ensembl 88 (GRCh38): external_gene_name "Gene name"
atts <- listAttributes(ensembl, page = "feature_page")
potential_gene_name <- union(match("Associated Gene Name", atts$description),
                             match("Gene name", atts$description))
potential_gene_name <- potential_gene_name[!is.na(potential_gene_name)]
att_gene_name <- atts[potential_gene_name, "name"]
stopifnot(is.character(att_gene_name), length(att_gene_name) == 1)

exons_all <- getBM(attributes = c("ensembl_gene_id", "ensembl_exon_id",
                                  "chromosome_name", "exon_chrom_start",
                                  "exon_chrom_end", "strand",
                                  att_gene_name,
                                  "gene_biotype"),
                   mart = ensembl)
exons_final <- exons_all[exons_all$chromosome_name %in% c(1:22, "X", "Y", "MT"),
                         c("ensembl_gene_id", "chromosome_name",
                           "exon_chrom_start", "exon_chrom_end",
                           "strand", att_gene_name,
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
