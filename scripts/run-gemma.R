#!/usr/bin/env Rscript

# Run GEMMA for eQTL analysis.
#
# Usage: Rscript run-gemma.R counts tss pca relat prefix_plink window dir_pheno
# dir_plink dir_gemma
#
# counts: gene-by-sample matrix of normalzied gene expression values
#
# tss: tab-separated file with transcription start site info for each gene in
# counts. No column names, but has the following columns: gene, chromosome, tss,
# strand, biotype.
#
# pca: tab-separated file with principal components to include as covariates.
# First column is 1 for intercept. No column names.
#
# relat: pairwise relatedness values in GEMMA format "id1 id2 value"
#
# prefix_plink: The prefix of the PLINK .bed, .bim, and .fam files.
#
# window: The window to search upstream and downstream of the TSS for SNPs.
#
# dir_pheno: Top-level directory to save phenotype file for each gene (required
# to subset PLINK data).
#
# dir_plink: Top-level directory to save PLINK files for each gene.
#
# dir_gemma: Top-level directory to store GEMMA results for each gene.
#
# The main output file is saved as dir_gemma/top-pca-{n_pcs}.txt, where n_pcs is
# the number of PCs included as covariates. It has the following columns:
# gene	chr	rs	ps	n_miss	allele1	allele0	af	l_mle	p_lrt	n_snps n_pcs

# Arguments---------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 9)
# counts
f_exp <- args[1]
# tss
f_tss <- args[2]
# pca
f_pca <- args[3]
# relat
f_relat <- args[4]
# prefix_plink
plink_prefix_all <- args[5]
# window
window <- as.numeric(args[6])
# dir_pheno
dir_pheno <- args[7]
# dir_plink
dir_plink <- args[8]
# dir_gemma
dir_gemma <- args[9]

# For testing:
# f_exp <- "../data/counts-normalized.txt"
# f_tss <- "../data/tss/tss.txt"
# f_pca <- "../data/pca/pca-5.txt"
# f_relat <- "../data/relatedness-matrix-all.txt"
# plink_prefix_all <- "../data/plink/cardioqtl"
# window <- 10^6
# dir_pheno <- "../data/pheno/"
# dir_plink <- "../data/plink/"
# dir_gemma <- "../data/gemma/"

f_plink_bed <- paste0(plink_prefix_all, ".bed")
f_plink_bim <- paste0(plink_prefix_all, ".bim")
f_plink_fam <- paste0(plink_prefix_all, ".fam")

# Check input ------------------------------------------------------------------

stopifnot(file.exists(f_exp, f_tss, f_pca, f_relat,
                      f_plink_bed, f_plink_bim, f_plink_fam))

# Import data files
exp <- read.delim(f_exp, check.names = FALSE)
tss <- read.table(f_tss, stringsAsFactors = FALSE)
pca <- read.table(f_pca)
relat <- read.delim(f_relat, stringsAsFactors = FALSE)
plink_bim <- read.table(f_plink_bim, stringsAsFactors = FALSE)
plink_fam <- read.table(f_plink_fam, stringsAsFactors = FALSE)

stopifnot(nrow(exp) == nrow(tss),
          ncol(exp) == nrow(pca),
          rownames(exp) == tss[, 1])

colnames(tss) <- c("gene", "chr", "tss", "strand", "biotype")
if (ncol(pca) == 1) {
  colnames(pca) <- "intercept"
} else {
  colnames(pca) <- c("intercept", paste0("pc", seq_along(colnames(pca)[-1])))
}
rownames(pca) <- colnames(exp)

# Check that samples are in PLINK fam file and relatedness matrix
stopifnot(all(colnames(exp) %in% relat$Ind1),
          all(colnames(exp) %in% relat$Ind1),
          all(colnames(exp) %in% plink_fam[, 2]))

# Determine valid chromosomes for eQTL mapping from PLINK .bim file
chr_valid <- as.character(unique(plink_bim[, 1]))

# Determine number of PCs being regressed
n_pcs <- ncol(pca) - 1
stopifnot(n_pcs >= 0, n_pcs <= ncol(exp))

# Create output subdirectories, e.g.
#
#         -- pca0
# pheno ---- pca1
#         -- pcaX
#
#         -- pca0
# plink ---- pca1
#         -- pcaX
#
#         -- pca0
# gemma ---- pca1
#         -- pcaX
#
pca_subdir <- paste0("pca", n_pcs, "/")
dir_pheno_pca <- file.path(dir_pheno, pca_subdir)
dir_plink_pca <- file.path(dir_plink, pca_subdir)
dir_gemma_pca <- file.path(dir_gemma, pca_subdir)
dir.create(dir_pheno_pca, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_plink_pca, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_gemma_pca, showWarnings = FALSE, recursive = TRUE)

# Create global output file
f_top <- file.path(dir_gemma, paste0("top-pca-", n_pcs, ".txt"))
f_top_colnames <- c("gene", "chr", "rs", "ps", "n_miss", "allele1", "allele0",
                      "af", "l_mle", "p_lrt", "n_snps", "n_pcs")
cat(f_top_colnames, file = f_top, sep = "\t")
cat("\n", file = f_top, sep = "", append = TRUE)

# Functions for eQTL analysis --------------------------------------------------

get_tss_gene <- function(chr, tss, window) {
  stopifnot(is.character(chr), is.numeric(tss), is.numeric(window))
  start <- tss - window
  if (start < 0)
    start <- 0
  end = tss + window
  out <- c(chr, start, end)
  names(out) <- c("chr", "start", "end")
  return(out)
}

# Format is:
# Column 1/2 = FID/IID
# Column 3 = Phenotype (i.e. gene expression levels)
create_pheno_file <- function(fid, iid, e, f = "") {
  stopifnot(length(fid) == 1 | length(fid) == length(e),
            length(iid) == 1 | length(iid) == length(e),
            is.numeric(e))
  d <- data.frame(fid, iid, e)
  write.table(d, file = f, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  return(invisible(f))
}

create_plink_g <- function(prefix_in, f_pheno, chr, from, to, prefix_out) {
  cmd <- sprintf("plink2 --bfile %s --make-bed --pheno %s --chr %s --from-bp %s --to-bp  %s --out %s",
                 prefix_in, f_pheno, chr, from, to, prefix_out)
  system(cmd)
  return(invisible(cmd))
}

run_gemma <- function(plink, relatedness, pca, out_prefix, outdir) {
  cmd <- sprintf("gemma -bfile %s -k %s -km 2 -c %s -lmm 2 -o %s",
                 plink, relatedness, pca, out_prefix)
  system(cmd)
  # Move to output directory
  cmd2 <- sprintf("mv output/%s* %s", out_prefix, outdir)
  system(cmd2)
  outfile <- paste0(outdir, out_prefix, ".assoc.txt")
  return(invisible(outfile))
}

parse_gemma <- function(f, gene, n_pcs, outfile = "") {
  gemma <- read.delim(f, stringsAsFactors = FALSE)
  n_snps <- nrow(gemma)
  top <- data.frame(gene, gemma[which.min(gemma$p_lrt), ], n_snps, n_pcs,
                    stringsAsFactors = FALSE)
  write.table(top, file = outfile, quote = FALSE, sep = "\t", row.names = FALSE)
  return(invisible(outfile))
}

# Per gene eQTL analysis -------------------------------------------------------

genes <- rownames(exp)
iid <- colnames(exp)
fid <- "HUTTERITES"

for (i in seq_along(genes)) {
  g <- rownames(exp)[i]
  message("\n\n####\ngene: ", g, "\n####\n\n")
  chr = tss$chr[i]
  if (!(chr %in% chr_valid)) {
    message("skipped:\tinvalid chromosome")
    next
  }
  tss_g <- get_tss_gene(chr = chr, tss = tss$tss[i], window = window)
  f_pheno_g <- paste0(dir_pheno_pca, g, ".pheno")
  create_pheno_file(fid = fid, iid = iid, e = as.numeric(exp[i, ]),
                    f = f_pheno_g)
  plink_prefix_g <- paste0(dir_plink_pca, g)
  message("Running PLINK")
  create_plink_g(prefix_in = plink_prefix_all,
                 f_pheno = f_pheno_g,
                 chr = tss_g["chr"],
                 from = tss_g["start"],
                 to = tss_g["end"],
                 prefix_out = plink_prefix_g)
  if (!file.exists(paste0(plink_prefix_g, ".bed"))) {
    message("\nskipped:\tNo PLINK output\n")
    next
  }
  message("Running GEMMA")
  f_gemma <- run_gemma(plink = plink_prefix_g, relatedness = f_relat, pca = f_pca,
                       out_prefix = paste0("pca", n_pcs, "-", g),
                       outdir = dir_gemma_pca)
  message("Parse GEMMA")
  f_gemma_top <- sub(".assoc.txt", ".top.txt", f_gemma)
  message(f_gemma_top)
  parse_gemma(f = f_gemma, gene = g, n_pcs = n_pcs, outfile = f_gemma_top)
  message(f_top)
  # Write to global results file
  system(sprintf("cat %s | sed -e '1d' >> %s", f_gemma_top, f_top))
}
