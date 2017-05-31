#!/usr/bin/Rscript

# Filter the phenotype file (IID is column 1) to remove any individuals used in
# the eQTL experiment, which are contained in the file to filter individuals
# from the PLINK data for the eQTL analysis (IID is column 2). Also sorts output
# by IID.
#
# usage: Rscript scripts/filter-gwas.R {input.raw} {input.filter} >
# {output.clean}

# Input
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2,
          file.exists(args))
f_raw <- args[1]
f_filter <- args[2]
# For testing:
# f_raw <- "../data/gwas/ykl/raw-ykl.txt"
# f_filter <- "../data/plink/filter-individuals.txt"
raw <- read.delim(f_raw, stringsAsFactors = FALSE)
filter <- read.table(f_filter, stringsAsFactors = FALSE)
colnames(filter) <- c("fid", "iid")

# Filter eQTL individuals
clean <- raw[!(raw$findiv %in% filter$iid), ]
stopifnot(nrow(clean) <= nrow(raw),
          nrow(clean) >= nrow(raw) - nrow(filter),
          ncol(clean) == ncol(raw))

# Sort by IID
clean <- clean[order(clean$findiv), ]

# Output
write.table(clean, file = "", quote = FALSE, sep = "\t", row.names = FALSE)
