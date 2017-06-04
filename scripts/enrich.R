# Are cardiomyocyte eQTLs enriched in signal from GWAS?

library("Hmisc")
library("readr")

# Functions --------------------------------------------------------------------

#' Main enrichment function.
#'
#' For each interval, calculate the enrichment compared to the background
#' enrichment.
#'
#' x - the variable that is split into intervals
#'
#' y - the variable that is being tested for enrichment
#'
#' cutoff - the value of y at which to separate genes when cacluting enrichment
#'
#' m - the average number of genes in each interval
#'
#' x_direction - With increasing values of x, should genes be included if they
#' are greater than or less than the metric that defines that interval. Default
#' "greater", anything else with do "lesser".
#'
#' cutoff_direction - When testing for enrichment of y at the given interval of
#' x, are genes counted if they are greater or lesser than the cuottf.
#'
#' @export
enrich <- function(x, y, cutoff, m = 50,
                   x_direction = "greater",
                   cutoff_direction = "greater") {
  intervals <- Hmisc::cut2(x, m = m, onlycuts = TRUE)
  enrichment <- numeric(length = length(intervals))
  sizes <- numeric(length = length(intervals))
  if (cutoff_direction == "greater") {
    background_enrich <- sum(y > cutoff) / length(y)
  } else {
    background_enrich <- sum(y < cutoff) / length(y)
  }
  for (i in seq_along(intervals)) {
    if (x_direction == "greater") {
      y_sub <- y[x > intervals[i]]
    } else {
      y_sub <- y[x < intervals[i]]
    }
    sizes[i] <- length(y_sub)
    # browser()
    if (cutoff_direction == "greater") {
      enrichment[i] <- sum(y_sub > cutoff) / sizes[i] / background_enrich
    } else {
      enrichment[i] <- sum(y_sub < cutoff) / sizes[i] / background_enrich
    }
  }
  return(data.frame(enrichment = enrichment, intervals = intervals, sizes = sizes))
}

plot_enrich <- function(x, ...) {
  plot(x$enrichment, type = "l", col = "red", xaxt = "n", ...)
  n <- length(x$enrichment)
  spacing <- seq(1, n, by = 50)
  axis(side = 1, at = (1:n)[spacing], padj = 0.5,
       labels = paste0(x$sizes, "\n(", round(x$intervals, digits = 2),
                       ")")[spacing])
}

bonferroni <- function(p, n) {
  # p - nominal p-value(s)
  # n - number(s) of tests
  stopifnot(is.numeric(p), p >= 0, p <= 1,
            is.numeric(n), n > 0)
  n_valid_len <- length(n) == 1 || length(n) == length(p)
  if (!n_valid_len)
    stop("n must be length 1 or the same length as p")
  p_b <- pmin(1, p * n)
  stopifnot(p_b <= 1)
  return(p_b)
}

# Data -------------------------------------------------------------------------

tri <- read_tsv("../data/gwas/tri/gemma-tri.assoc.txt")
neu <- read_tsv("../data/gwas/neu/gemma-neu.assoc.txt")
ykl <- read_tsv("../data/gwas/ykl/gemma-ykl.assoc.txt")
chi <- read_tsv("../data/gwas/chi/gemma-chi.assoc.txt")

eqtl <- read_delim("../data/gemma/top-pca-0.txt", skip = 1, delim = "\t",
                   col_names = FALSE)
colnames(eqtl) <- c("gene", "chr", "rs", "ps", "n_miss",
                    "allele1", "allele0", "af", "beta",
                    "se", "l_remle", "l_mle", "p_wald",
                    "p_lrt", "p_score", "n_snps", "n_pcs")


eqtl$p_bonf <- bonferroni(eqtl$p_lrt, eqtl$n_snps)
eqtl$p_bh = p.adjust(eqtl$p_bonf, method = "BH")

# Enrichment--------------------------------------------------------------------

eqtl_tri <- merge(eqtl, tri,
                  by = c("chr", "rs", "ps", "allele1", "allele0"),
                  suffixes = c(".eqtl", ".tri"))

tri_enrich <- enrich(x = eqtl_tri$p_lrt.eqtl,
                     y = eqtl_tri$p_lrt.tri,
                     cutoff = 0.05,
                     m = 100,
                     x_direction = "lesser",
                     cutoff_direction = "lesser")

plot_enrich(tri_enrich, main = "tri")

eqtl_neu <- merge(eqtl, neu,
                  by = c("chr", "rs", "ps", "allele1", "allele0"),
                  suffixes = c(".eqtl", ".neu"))

neu_enrich <- enrich(x = eqtl_neu$p_lrt.eqtl,
            y = eqtl_neu$p_lrt.neu,
            cutoff = 0.05,
            m = 100,
            x_direction = "lesser",
            cutoff_direction = "lesser")

plot_enrich(neu_enrich, main = "neu")

eqtl_ykl <- merge(eqtl, ykl,
                  by = c("chr", "rs", "ps", "allele1", "allele0"),
                  suffixes = c(".eqtl", ".ykl"))

ykl_enrich <- enrich(x = eqtl_ykl$p_lrt.eqtl,
                     y = eqtl_ykl$p_lrt.ykl,
                     cutoff = 0.05,
                     m = 100,
                     x_direction = "lesser",
                     cutoff_direction = "lesser")

plot_enrich(ykl_enrich, main = "ykl")

eqtl_chi <- merge(eqtl, chi,
                  by = c("chr", "rs", "ps", "allele1", "allele0"),
                  suffixes = c(".eqtl", ".chi"))

chi_enrich <- enrich(x = eqtl_chi$p_lrt.eqtl,
                     y = eqtl_chi$p_lrt.chi,
                     cutoff = 0.05,
                     m = 100,
                     x_direction = "lesser",
                     cutoff_direction = "lesser")

plot_enrich(chi_enrich, main = "chi")
