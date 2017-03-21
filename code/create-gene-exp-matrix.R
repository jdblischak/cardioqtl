#!/usr/bin/env Rscript

# Create gene expression matrix from kallisto results.
#
# Usage:
#
# Rscript create-gene-exp-matrix.R <tx2gene> <kallisto> > out.txt
#
# tx2gene - Path to tab-delimited file with ENST in col 1 and ENSG in col 2
# kallisto - Path to directory containing kallisto output
#
# Use tximport to convert kallisto transcript-level TPM to gene-level
# lengthScaledTPM. A plain-text, tab-delimited file is written to stdout.

library("tximport")
library("readr")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
tx2gene_path <- args[1]
stopifnot(file.exists(tx2gene_path))
kallisto <- args[2]
stopifnot(dir.exists(kallisto))

tx2gene <- read.delim(tx2gene_path, stringsAsFactors = FALSE)

files <- Sys.glob(file.path(kallisto, "*", "abundance.tsv"))
suppressMessages(
  tx_list <- tximport(files,
                      type = "kallisto",
                      countsFromAbundance = "lengthScaledTPM",
                      tx2gene = tx2gene,
                      reader = read_tsv)
)

# tximport discards the names. The sample name is the directory that contains
# the data file.
samples <- basename(dirname(files))
tx_list <- lapply(tx_list,
                  function(x) {if (is.matrix(x)) colnames(x) <- samples; x})

write.table(tx_list$counts, file = "", quote = FALSE, sep = "\t",
            row.names = TRUE)
