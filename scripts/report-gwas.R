#!/usr/bin/env Rscript

# Report the GWAS results for a given phenotype by running it through a
# parameterized R Markdown document.
#
# usage: Rscript report-gwas.R <template> <trait> <outfile>
#
# template - parameterized R Markdown document
#
# trait - 3-digit code for trait
#
# outfile - File name to save the HTML

library("magrittr")
library("rmarkdown")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
template <- args[1]
trait <- args[2]
outfile <- args[3]
stopifnot(file.exists(template),
          nchar(trait) == 3,
          dir.exists(dirname(outfile)))

# For testing:
# template <- "../scratch/gwas.Rmd"
# trait <- "tri"
# outfile <- "../data/gwas/tri/gwas-tri.html"

# render creates an intermediate Markdown file which always has the same name.
# This causes chaos when convert the Markdown to HTML with pandoc. The traits
# get randomized and some even fail because the file is deleted by another
# process. Need to create a unique temporary file for the template.
template_tmp <- template %>%
  tools::file_path_sans_ext() %>%
  paste0("-", trait, ".Rmd")

file.copy(from = template, to = template_tmp)

render(input = template_tmp,
       output_file = outfile,
       params = list(trait = trait))

file.remove(template_tmp)
