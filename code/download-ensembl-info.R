#!/usr/bin/env Rscript

# Create Ensembl transcript mappings to genes.
#
# Usage:
#
# Rscript download-ensembl-info.R <archive> > out.txt
#
# archive - URL of Ensembl archive, e.g. mar2015.archive.ensembl.org
#
# A plain-text, tab-delimited file is written to stdout.

suppressMessages(library("biomaRt"))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
archive <- args[1]

ensembl <- useMart(host = archive,
                   biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "hsapiens_gene_ensembl")
tx2gene <- getBM(attributes = c("ensembl_transcript_id",
                                "ensembl_gene_id",
                                "chromosome_name",
                                "external_gene_name",
                                "gene_biotype",
                                "source"),
                 mart = ensembl)

write.table(tx2gene, file = "", quote = FALSE, sep = "\t",
            row.names = FALSE)
