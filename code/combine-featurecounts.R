#!/usr/bin/env Rscript

library("edgeR")
library("magrittr")
library("stringr")

f_counts <- Sys.glob("/scratch/midway2/jdblischak/cardioqtl-counts/*txt")
raw <- readDGE(f_counts, columns = c(1, 7), comment.char = "#")
counts <- as.matrix(raw)
colnames(counts) <- colnames(counts) %>%
  str_split("/", simplify = TRUE) %>%
  `[`(, ncol(.)) %>%
  str_replace(".genecounts", "")

head(counts)
write.table(counts, file = "../data/counts-subread.txt", quote = FALSE, sep = "\t")
