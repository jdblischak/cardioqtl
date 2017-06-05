# Are cardiomyocyte eQTLs enriched in signal from GWAS?

library("cowplot")
library("flux")
library("ggplot2")
library("Hmisc")
library("magrittr")
library("plyr"); library("dplyr")
library("RColorBrewer")
library("readr")
library("stringr")
theme_set(theme_cowplot())

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


eqtl <- read_delim("../data/gemma/top-pca-0.txt", skip = 1, delim = "\t",
                   col_names = FALSE)
colnames(eqtl) <- c("gene", "chr", "rs", "ps", "n_miss",
                    "allele1", "allele0", "af", "beta",
                    "se", "l_remle", "l_mle", "p_wald",
                    "p_lrt", "p_score", "n_snps", "n_pcs")


eqtl$p_bonf <- bonferroni(eqtl$p_lrt, eqtl$n_snps)
eqtl$p_bh = p.adjust(eqtl$p_bonf, method = "BH")

f_gwas <- Sys.glob("../data/gwas/*/gemma-*.assoc.txt")
eqtl_gwas <- eqtl
for (f in f_gwas) {
  n <- f %>% basename %>%
    str_replace("gemma-", "") %>%
    str_replace(".assoc.txt", "")
  d <- read_tsv(f)
  p <- d$p_lrt
  names(p) <- d$rs
  eqtl_gwas <- cbind(eqtl_gwas, p[eqtl_gwas$rs])
  colnames(eqtl_gwas)[ncol(eqtl_gwas)] <- n
}

eqtl_gwas_final <- na.omit(eqtl_gwas)

# Enrichment--------------------------------------------------------------------

l_enrich <- list()
for (i in 20:39) {
  n <- colnames(eqtl_gwas_final)[i]
  l_enrich[[n]] <- enrich(x = -log10(eqtl_gwas_final$p_lrt),
                          y = eqtl_gwas_final[, i],
                          cutoff = 0.05,
                          m = 100,
                          x_direction = "greater",
                          cutoff_direction = "lesser")
  plot_enrich(l_enrich[[n]], main = n, ylim = c(0, 2))
}

# Plot with base ---------------------------------------------------------------

heart <- c("cho", "cim", "dbp", "hdl", "lav", "ldl", "lvm", "sbp", "tri")
immune <- c("chi", "eos", "ige", "lym", "mon", "neu", "ykl")
lung <- c("bri", "eno", "fev", "fvc")
stopifnot(length(c(heart, immune, lung)) == 20)

# Heart
plot_enrich(l_enrich[[heart[1]]], ylim = c(0, 2), main = "Cardiovascular")
for (pheno in heart[-1]) {
  lines(l_enrich[[pheno]]$enrichment)
}

# Immune
plot_enrich(l_enrich[[immune[1]]], ylim = c(0, 2), main = "Immunological")
for (pheno in immune[-1]) {
  lines(l_enrich[[pheno]]$enrichment)
}

# Lung
plot_enrich(l_enrich[[lung[1]]], ylim = c(0, 2), main = "Pulmonary")
for (pheno in lung[-1]) {
  lines(l_enrich[[pheno]]$enrichment)
}

# Plot with ggplot2 ------------------------------------------------------------

d_enrich <- ldply(l_enrich, .id = "pheno")
# Remove NaN for final bin with zero eQTLs
d_enrich <- d_enrich[!is.nan(d_enrich$enrichment), ]
len <- sum(d_enrich$pheno == "bri")
d_enrich$index <- seq(len)
d_enrich$label <- sprintf("%.2f\n(%d)", d_enrich$intervals, d_enrich$sizes)
breaks <- seq(1, len, length.out = 10)
labels <- d_enrich$label[d_enrich$pheno == "bri"][breaks]

p_heart <- ggplot(d_enrich[d_enrich$pheno %in% heart, ],
                  aes(x = index, y = enrichment, color = pheno)) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_x_continuous(breaks = breaks, labels = labels) +
  ylim(min(d_enrich$enrichment),
       max(d_enrich$enrichment)) +
  scale_color_brewer(name = NULL, palette = "Reds") +
  theme(legend.position = "none") +
  labs(x = "-log10 P cutoff for eQTL\n(Number of eQTLs)",
       y = "Fold enrichment of GWAS P < 0.05",
       title = "Cardiovascular traits")

p_immune <- p_heart %+% d_enrich[d_enrich$pheno %in% immune, ] +
  scale_color_brewer(name = NULL, palette = "Blues") +
  labs(title = "Immunological traits")

p_lung <- p_heart %+% d_enrich[d_enrich$pheno %in% lung, ] +
  scale_color_brewer(name = NULL, palette = "Purples") +
  labs(title = "Pulmonary traits")

plot_grid(p_heart +
            annotate("text", x = rep(len + 1, length(heart)),
                     y = d_enrich$enrichment[d_enrich$pheno %in% heart &
                                             d_enrich$index == len],
                     label = heart),
          p_immune +
            annotate("text", x = rep(len + 1, length(immune)),
                     y = d_enrich$enrichment[d_enrich$pheno %in% immune &
                                             d_enrich$index == len],
                     label = immune),
          p_lung +
            annotate("text", x = rep(len + 1, length(lung)),
                     y = d_enrich$enrichment[d_enrich$pheno %in% lung &
                                             d_enrich$index == len],
                     label = lung))

# Calculate AUC ----------------------------------------------------------------

background <- auc(x = 1:len, y = rep(1, len))
d_auc <- d_enrich %>%
  group_by(pheno) %>%
  summarize(auc_raw = auc(x = index, y = enrichment),
            auc_std = auc_raw - background,
            enr_max = max(enrichment),
            enr_final = enrichment[length(enrichment)])
set1 <- brewer.pal(n = 4, "Set1")
names(set1) <- c("red", "blue", "green", "purple")
pheno_colors <- character(length = length(l_enrich))
for (i in seq_along(pheno_colors)) {
  if (names(l_enrich)[i] %in% heart) {
    pheno_colors[i] <- set1["red"]
  } else if (names(l_enrich)[i] %in% immune) {
    pheno_colors[i] <- set1["blue"]
  } else {
    pheno_colors[i] <- set1["purple"]
  }
}

p_auc <- ggplot(d_auc,
       aes(x = reorder(pheno, auc_std), y = auc_std, fill = pheno)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = pheno_colors) +
  theme(legend.position = "none") +
  labs(x = "Trait", y = "Enrichment above background\n(area under the curve)")

p_max <- p_auc %+% aes(x = reorder(pheno, enr_max), y = enr_max) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "Trait", y = "Maximum enrichment")

p_final <- p_max %+% aes(x = reorder(pheno, enr_final), y = enr_final) +
  labs(x = "Trait", y = "Final enrichment")

plot_grid(p_auc + geom_hline(yintercept = 0, linetype = "dashed"),
          p_max,
          p_final)
