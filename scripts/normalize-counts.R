#!/usr/bin/env Rscript

# Normalize the counts from Subread.
#
# 1. Normalize to N(0, 1), i.e. standard normal, within each sample (column).
# 2. Normalize to N(0, 1), i.e. standard normal, across each gene (row).
#
# Usage: Rscript normalize-counts.R input > output
#
# The input file is a tab-separated matrix of gene counts with a header row.
#
# The results are exported as tab-separated columns to standard out.

# Functions --------------------------------------------------------------------

# Normalize a vector to N(0, 1)
norm_to_standard <- function(x) {
  qq <- qqnorm(x, plot.it = FALSE)
  out <- qq$x
  stopifnot(abs(mean(out)) < 0.1,
            abs(var(out) - 1) < 0.1)
  return(out)
}

# Check if two vectors have the same ranks
compare_ranks <- function(x, y) {
  x_rank <- rank(x, ties.method = "first")
  y_rank <- rank(y, ties.method = "first")
  return(all(x_rank == y_rank))
}

# Specify the number of digits to output
# http://stackoverflow.com/questions/3443687/formatting-decimal-places-in-r/12135122#12135122
specify_decimal <- function(x, k) format(round(x, k), nsmall = k)

# Input ------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1, file.exists(args))
input <- args
# For testing:
# input <- "../data/counts-clean.txt"
raw <- read.delim(input, check.names = FALSE)

# Normalize within each sample (column) ----------------------------------------
n_within_sample <- apply(raw, 2, norm_to_standard)
n_within_sample <- as.data.frame(n_within_sample)

# Normalization should not have changed rank of genes within each sample
ranks_within_sample <- mapply(compare_ranks, raw, n_within_sample)
stopifnot(ranks_within_sample)

# Normalize across each gene (row) ----------------------------------------
n_across_gene <- apply(n_within_sample, 1, norm_to_standard)

# Fix row and column names
# Annoyingly, R transposes the output, and removes row and column names
n_across_gene <- t(n_across_gene)
n_across_gene <- as.data.frame(n_across_gene)
stopifnot(dim(n_across_gene) == dim(raw))
rownames(n_across_gene) <- rownames(raw)
colnames(n_across_gene) <- colnames(raw)

# Normalization should not have changed rank of genes across each gene
ranks_across_gene <- mapply(compare_ranks,
                            as.data.frame(t(n_within_sample)),
                            as.data.frame(t(n_across_gene)))
stopifnot(ranks_across_gene)

# Output -----------------------------------------------------------------------
# Truncate to 6 decimal places
write.table(specify_decimal(n_across_gene, 6), file = "",
            quote = FALSE, sep = "\t")
