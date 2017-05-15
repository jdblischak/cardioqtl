#!/usr/bin/env Rscript

# Obtain the most upstream TSS per gene.
#
# Usage:
#   Rscript get-tss.R <archive> input > output.bed
#
# archive - URL of Ensembl archive, e.g. mar2015.archive.ensembl.org
# input - gene expression matrix. Only rownames are used
# output.bed - BED format with one TSS per gene

suppressMessages(library("biomaRt"))
library("stringr")

# Input ------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
archive <- args[1]
input <- args[2]
# For testing:
# archive <- "mar2015.archive.ensembl.org"
# input <- "../data/counts-clean.txt"
stopifnot(file.exists(input))
raw <- read.delim(input, check.names = FALSE)
genes <- rownames(raw)
stopifnot(str_sub(genes, 1, 4) == "ENSG")

# Download TSS from Ensembl ----------------------------------------------------
ensembl <- useMart(host = archive,
                   biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "hsapiens_gene_ensembl")

tss_all <- getBM(attributes = c("ensembl_gene_id",
                                "chromosome_name",
                                "transcription_start_site",
                                "start_position",
                                "end_position",
                                "strand",
                                "external_gene_name",
                                "gene_biotype"),
                 filters = "ensembl_gene_id",
                 values = genes,
                 mart = ensembl)

# Check download ---------------------------------------------------------------

# The Ensembl variables can be cryptic. Confirming my intuition in my SO answer:
# http://stackoverflow.com/questions/13012210/find-transcription-start-sites-with-biomart/35544659#comment58776614_25434411

is_unique <- function(x) length(unique(x)) == 1

# All transcripts of a gene should have the same start_position, i.e. the 5´
# most recorded position for any of the transcripts
stopifnot(tapply(tss_all$start_position, tss_all$ensembl_gene_id, is_unique))

# All transcripts of a gene should have the same end_position, i.e. the 3´
# most recorded position for any of the transcripts
stopifnot(tapply(tss_all$end_position, tss_all$ensembl_gene_id, is_unique))

# Therefore start_position should always be less than end_position, no matter
# the strand of the gene
stopifnot(mapply(function(x, y) x < y,
                 tss_all$start_position,
                 tss_all$end_position))

# Reduce to most upstream TSS per gene -----------------------------------------

# Detach biomaRt before loading dplyr because they both define `select`
detach("package:biomaRt")
suppressPackageStartupMessages(library("dplyr"))

tss_per_gene <- tss_all %>%
  group_by(ensembl_gene_id, chromosome_name,
           strand, external_gene_name, gene_biotype) %>%
  summarize(tss = if (all(strand == 1)) min(start_position) else max(end_position)) %>%
  ungroup() %>%
  arrange(ensembl_gene_id)

# Convert to BED format --------------------------------------------------------

tss_bed <- tss_per_gene %>%
  rename(chr = chromosome_name) %>%
  mutate(chr = ifelse(chr == "MT", "M", chr),
         chr = paste0("chr", chr),
         start = tss - 1,
         end = tss,
         strand = ifelse(strand == 1, "+", "-")) %>%
  select(chr, start, end, ensembl_gene_id, strand)

# Output -----------------------------------------------------------------------
write.table(tss_bed, "", quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)
