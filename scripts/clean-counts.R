#!/usr/bin/env Rscript

# Clean the counts from Subread.
#
# 1. Remove invalid samples:
#      - Too few reads: 26302, 110232
#      - Sample swap: 160001
# 2. Remove genes with zero counts
# 2. Remove lowly expressed genes. Keep genes with median log2 cpm > 0
#
# Usage: Rscript clean-counts.R input > output
#
# The input file is a tab-separated matrix of gene counts with a header row.
#
# The results are exported as tab-separated columns to standard out.

suppressPackageStartupMessages(library("edgeR"))

# Input ------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1, file.exists(args))
input <- args
raw <- read.delim(input, check.names = FALSE)

# Remove outliers --------------------------------------------------------------
to_remove <- c("26302", "110232", "160001")
e <- raw[, !(colnames(raw) %in% to_remove)]

# Remove zeros -----------------------------------------------------------------
zeros <- apply(e, 1, function(x) all(x == 0))
# sum(zeros)
e_detect <- e[!zeros, ]

# Remove lowly expressed genes -------------------------------------------------
cpm_all <- cpm(e_detect, log = TRUE)
cpm_all_median <- apply(cpm_all, 1, median)
# plot(density(cpm_all_median))
# abline(0, 1, col = "red")
# sum(cpm_all_median > 0)
e_expr <- e_detect[cpm_all_median > 0, ]

# Output -----------------------------------------------------------------------
write.table(e_expr, file = "", quote = FALSE, sep = "\t")
