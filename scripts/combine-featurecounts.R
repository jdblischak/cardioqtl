#!/usr/bin/env Rscript

# Combine the featureCounts output per sample into one matrix.
#
# Usage: Rscript combine-featurecounts.R [file, file, ...] > output
#
# Each input file is the featureCounts results for one sample. The gene names
# are extracted from column 1 and the counts from column 7.
#
# The results are exported as tab-separated columns to standard out.

suppressPackageStartupMessages(library("edgeR"))
library("magrittr")
library("stringr")

# Input files ------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) > 0, file.exists(args))
f_counts <- args
# For testing:
# f_counts <- Sys.glob("/scratch/midway2/jdblischak/cardioqtl-counts/*txt")
raw <- readDGE(f_counts, columns = c(1, 7), comment.char = "#")
counts <- as.matrix(raw)

# Change column names from filenames to sample names ---------------------------
colnames(counts) <- colnames(counts) %>%
  str_split("/", simplify = TRUE) %>%
  `[`(, ncol(.)) %>%
  str_replace(".genecounts", "")

# Sort the rownames and colnames -----------------------------------------------
counts <- counts[order(rownames(counts)), order(as.numeric(colnames(counts)))]

# Write to standard out --------------------------------------------------------
# head(counts)
write.table(counts, file = "", quote = FALSE, sep = "\t")
