#!/usr/bin/Rscript

library("magrittr")
library("stringr")

# Input ------------------------------------------------------------------------

# List of fastq files
fq_in <- Sys.glob("/project2/gilad/jdblischak/dox/fastq/*-0.000*.fastq.gz")
stopifnot(length(fq_in) == 46 * 2)

# Conversion table between cell line and ID
cell2id <- read.delim("/home/jdblischak/dox/data/cell-line-to-findiv.txt")

# Process filenames ------------------------------------------------------------

fq_in_df <- fq_in %>% 
  basename() %>%
  str_replace(".fastq.gz", "") %>%
  str_split_fixed("-", n = 5) %>%
  as.data.frame(stringsAsFactors = FALSE)

stopifnot(is.data.frame(fq_in_df), !is.na(fq_in_df))
colnames(fq_in_df) <- c("sample", "cell", "conc", "flow_cell", "lane")
stopifnot(fq_in_df$conc == "0.000")

fq_in_df$cell_line <- fq_in_df$cell %>%
  str_replace("c", "") %>%
  str_split("\\.") %>%
  sapply(function(x) x[1]) %>%
  as.numeric()

fq_in_df_id <- merge(fq_in_df, cell2id, all.x = TRUE, sort = FALSE)
stopifnot(fq_in_df_id[, c("cell", "flow_cell", "lane")] == 
          fq_in_df[, c("cell", "flow_cell", "lane")])

fq_out <- paste(fq_in_df_id$findiv, fq_in_df_id$flow_cell, fq_in_df_id$lane,
                sep = "-")
fq_out <- paste0("/project2/gilad/jdblischak/cardioqtl/data/fastq/", fq_out,
                 ".fastq.gz")

# Create symlinks --------------------------------------------------------------

dir.create("/project2/gilad/jdblischak/cardioqtl/data/fastq/", recursive = TRUE)

file.symlink(from = fq_in, to = fq_out)
