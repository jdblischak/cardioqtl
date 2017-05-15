#!/usr/bin/env Rscript

# Parse the results from GEMMA. Keep the SNP with the lowest p-value.
#
# Usage: Rscript parse-gemma.R gene gene.assoc.txt > gene.top.txt

# Input ------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
gene <- args[1]
f <- args[2]
# For testing:
# gene <- "ENSG00000070404"
# f <- "../data/gemma/ENSG00000070404.assoc.txt"
stopifnot(substr(gene, 1, 4) == "ENSG",
          file.exists(f))
gemma <- read.delim(f, stringsAsFactors = FALSE)

# Get top SNP ------------------------------------------------------------------

n_snps <- nrow(gemma)
top <- data.frame(gene, gemma[which.min(gemma$p_lrt), ], n_snps,
                  stringsAsFactors = FALSE)

# Output -----------------------------------------------------------------------
write.table(top, file = "", quote = FALSE, sep = "\t", row.names = FALSE)
