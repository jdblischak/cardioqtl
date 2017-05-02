#!/usr/bin/env Rscript

# Subset the individuals in the relatedness matrix.
#
# 1. Remove invalid samples:
#      - Too few reads: 26302, 110232
#      - Sample swap: 160001
# 2. Remove genes with zero counts
# 2. Remove lowly expressed genes. Keep genes with median log2 cpm > 0
#
# Usage: Rscript subset-relatedness.R exp-mat relatedness-all > relatedness-sub
#
# exp-mat - gene expression matrix. Extracts samples from header row
# relatedness-all - tab-separated input file with pairwise relatedness
#  measurements for all individuals.
# relatedness-sub - tab-separated output file with pairwise relatedness
#  measurements for individuals in this study.

# Input ------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2, file.exists(args))
exp_file <- args[1]
relat_file <- args[2]
exp_file <- "../data/counts-clean.txt"
relat_file <- "../data/relatedness-matrix-all.txt"
exp <- read.delim(exp_file, check.names = FALSE)
inds <- colnames(exp)
relat <- read.delim(relat_file, check.names = FALSE, stringsAsFactors = FALSE)

# Subset individuals -----------------------------------------------------------
ind1_keep <- relat$Ind1 %in% inds
ind2_keep <- relat$Ind2 %in% inds
relat_sub <- relat[ind1_keep & ind2_keep, ]

# Sanity checks ----------------------------------------------------------------

inds_sub <- c(relat_sub$Ind1, relat_sub$Ind2)
stopifnot(inds %in% unique(inds_sub),
          table(inds_sub) == ncol(exp) + 1, # Plus 1 for comparison to self
          relat_sub[relat_sub$Ind1 == relat_sub$Ind2, 3] > 1) # not sure why self-self isn't exactly 1

# Output -----------------------------------------------------------------------
write.table(relat_sub, file = "", quote = FALSE, sep = "\t",
            row.names = FALSE)
