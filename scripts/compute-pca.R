#!/usr/bin/env Rscript

# Compute principal components.
#
# Usage: Rscript compute-pca.R n_pcs input > output
#
# n_pcs: The number of PCs to include in output file.
#
# input: A gene-by-sample matrix of normalized gene expression levels.
#
# output: A sample-by-pc matrix of covariates written to standard out. The first
# column is 1 for the intercept in the linear model (format for GEMMA). Does not
# contain a header row.

set.seed(12345)

# Functions --------------------------------------------------------------------

# PCA
#
# x - a gene-by-sample matrix of gene expression levels
run_pca <- function(x) prcomp(t(x), retx = TRUE, center = TRUE, scale. = TRUE)

# Specify the number of digits to output
# http://stackoverflow.com/questions/3443687/formatting-decimal-places-in-r/12135122#12135122
specify_decimal <- function(x, k) format(round(x, k), nsmall = k)

# Input ------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
n_pcs <- as.numeric(args[1])
input <- args[2]
stopifnot(file.exists(input))
# For testing:
# n_pcs <- 10
# input <- "../data/counts-normalized.txt"
raw <- read.delim(input, check.names = FALSE)
stopifnot(n_pcs >= 0, n_pcs <= ncol(raw))

# Run PCA ----------------------------------------------------------------------
if (n_pcs > 0) {
  pca <- run_pca(raw)
  out <- pca$x[, seq(n_pcs), drop = FALSE]
  stopifnot(rownames(out) == colnames(raw), is.matrix(out))
}

# Output -----------------------------------------------------------------------

# The first column needs to be a vector of 1's to represent the intercept for
# GEMMA. First convert to character matrix to specify the number of decimal
# places (6) used for the PCs.
if (n_pcs > 0) {
  out_character <- apply(out, 2, specify_decimal, k = 6)
  out_character <- cbind("1", out_character)
} else {
  out_character <- matrix("1", nrow = ncol(raw))
}

write.table(out_character, file = "",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
