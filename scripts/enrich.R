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
theme_set(theme_cowplot(font_size = 10))

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

atac <- read_tsv("../data/atac/intersect.txt")

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

# 1000 permutations of a GWAS data set
set.seed(12345)
perms <- 100
l_permute <- vector(mode = "list", length = perms)
# All the GWAS p-values are uniformly distributed, so it shouldn't matter which
# one I permute.
for (i in 1:perms) {
  l_permute[[i]] <- enrich(x = -log10(eqtl_gwas_final$p_lrt),
                           y = sample(eqtl_gwas_final[, "bri"]),
                           cutoff = 0.05,
                           m = 100,
                           x_direction = "greater",
                           cutoff_direction = "lesser")
  if (i %% 10 == 0)
    plot_enrich(l_permute[[i]], main = i, ylim = c(0, 2))
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

# Permutations
plot_enrich(l_permute[[1]], ylim = c(0, 3), main = "Permutations")
for (i in 2:perms) {
  lines(l_permute[[i]]$enrichment)
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

# Calculate 95% confidence intervals based on the permutations
names(l_permute) <- paste0("perm", 1:perms)
d_permute <- ldply(l_permute, .id = "permutation")
# Remove NaN for final bin with zero eQTLs
d_permute <- d_permute[!is.nan(d_permute$enrichment), ]
len_permute <- sum(d_permute$permutation == "perm1")
stopifnot(len_permute == len)
d_permute$index <- seq(len_permute)
d_permute_ci <- d_permute %>%
  group_by(index) %>%
  summarize(lower_bound = quantile(enrichment, 0.05),
            upper_bound = quantile(enrichment, 0.95))

p_heart <- ggplot(d_enrich[d_enrich$pheno %in% heart, ],
                  aes(x = index, y = enrichment, color = pheno)) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_x_continuous(breaks = breaks, labels = labels) +
  ylim(min(d_enrich$enrichment),
       max(d_enrich$enrichment)) +
  scale_color_brewer(name = NULL, palette = "Reds") +
  geom_line(mapping = aes(x = index, y = lower_bound),
            data = d_permute_ci,
            col = "black", linetype = "dashed", size = 1.25) +
  geom_line(mapping = aes(x = index, y = upper_bound),
            data = d_permute_ci,
            col = "black", linetype = "dashed", size = 1.25) +
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
                     label = lung),
          nrow = 1, labels = LETTERS[1:3])
ggsave("enrich-phenos.png", width = 18, height = 6)

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

d_permute_auc <- d_permute %>%
  group_by(permutation) %>%
  summarize(auc_raw = auc(x = index, y = enrichment),
            auc_std = auc_raw - background,
            enr_max = max(enrichment),
            enr_final = enrichment[length(enrichment)])
boxplot(d_permute_auc$auc_std)
lapply(d_auc$auc_std, points, col = "red", pch = 19)
quantile(d_permute_auc$auc_std, probs = c(0.05, 0.95))

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
          p_final,
          nrow = 1, labels = LETTERS[1:3])
ggsave("enrich-phenos-auc.png", width = 18, height = 6)

# Enrichment of ATAC-seq peaks -------------------------------------------------

eqtl_atac <- merge(eqtl, atac, by = "rs", all.x = TRUE)
stopifnot(nrow(eqtl_atac) == nrow(eqtl))
# Convert NA to 0
for (cell in c("cm", "ips", "lcl")) {
  eqtl_atac[, cell] <- ifelse(is.na(eqtl_atac[, cell]),
                              0, eqtl_atac[, cell])
}
sum(eqtl_atac$cm)
sum(eqtl_atac$ips)
sum(eqtl_atac$lcl)

l_enrich_atac <- list()
for (i in 20:22) {
  n <- colnames(eqtl_atac)[i]
  l_enrich_atac[[n]] <- enrich(x = -log10(eqtl_atac$p_lrt),
                          y = eqtl_atac[, i],
                          cutoff = 0.5, # Just a hack to differentiate between
                          m = 100,      # 1 and 0 (i.e. presence or absence)
                          x_direction = "greater",
                          cutoff_direction = "greater")
  plot_enrich(l_enrich_atac[[n]], main = n, ylim = c(0, 2))
}

d_enrich_atac <- ldply(l_enrich_atac, .id = "cell")
# Remove NaN for final bin with zero eQTLs
d_enrich_atac <- d_enrich_atac[!is.nan(d_enrich_atac$enrichment), ]
len_atac <- sum(d_enrich_atac$cell == "cm")
d_enrich_atac$index <- seq(len_atac)
d_enrich_atac$label <- sprintf("%.2f\n(%d)", d_enrich_atac$intervals, d_enrich_atac$sizes)
breaks <- seq(1, len_atac, length.out = 10)
labels <- d_enrich_atac$label[d_enrich_atac$cell == "cm"][breaks]

cell_colors <- set1[c("red", "green", "blue")]
# ggplot2 doesn't accept manual color vectors with a names attribute
names(cell_colors) <- NULL

# Calculate permutations for ATAC
set.seed(12345)
perms_atac <- 100
l_permute_atac <- vector(mode = "list", length = perms_atac)
for (i in 1:perms_atac) {
  l_permute_atac[[i]] <- enrich(x = -log10(eqtl_atac$p_lrt),
                                y = sample(eqtl_atac[, "cm"]),
                                cutoff = 0.5, # Just a hack to differentiate between
                                m = 100,      # 1 and 0 (i.e. presence or absence)
                                x_direction = "greater",
                                cutoff_direction = "greater")
  if (i %% 10 == 0)
    plot_enrich(l_permute_atac[[i]], main = i, ylim = c(0, 2))
}

# Calculate 95% confidence intervals based on the permutations
names(l_permute_atac) <- paste0("perm", 1:perms)
d_permute_atac <- ldply(l_permute_atac, .id = "permutation")
# Remove NaN for final bin with zero eQTLs
d_permute_atac <- d_permute_atac[!is.nan(d_permute_atac$enrichment), ]
len_permute_atac <- sum(d_permute_atac$permutation == "perm1")
stopifnot(len_permute_atac == len_atac)
d_permute_atac$index <- seq(len_permute_atac)
d_permute_atac_ci <- d_permute_atac %>%
  group_by(index) %>%
  summarize(lower_bound = quantile(enrichment, 0.05),
            upper_bound = quantile(enrichment, 0.95))

p_atac <- ggplot(d_enrich_atac,
                  aes(x = index, y = enrichment, color = cell)) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_x_continuous(breaks = breaks, labels = labels) +
  ylim(min(d_enrich_atac$enrichment),
       max(d_enrich_atac$enrichment)) +
  scale_color_manual(name = "Cell type", values = cell_colors,
                     labels = c("CM", "iPSC", "LCL")) +
  geom_line(mapping = aes(x = index, y = lower_bound),
            data = d_permute_atac_ci,
            col = "black", linetype = "dashed", size = 1.25) +
  geom_line(mapping = aes(x = index, y = upper_bound),
            data = d_permute_atac_ci,
            col = "black", linetype = "dashed", size = 1.25) +
  labs(x = "-log10 P cutoff for eQTL\n(Number of eQTLs)",
       y = "Fold enrichment of eQTLs in ATAC-seq peaks",
       title = "Enrichment of eQTLs in open chromatin")

# Calculate AUC
background_atac <- auc(x = 1:len_atac, y = rep(1, len_atac))
d_auc_atac <- d_enrich_atac %>%
  group_by(cell) %>%
  summarize(auc_raw = auc(x = index, y = enrichment),
            auc_std = auc_raw - background_atac,
            enr_max = max(enrichment),
            enr_final = enrichment[length(enrichment)])
# Calculate AUC of ATAC permutations
d_permute_atac_auc <- d_permute_atac %>%
  group_by(permutation) %>%
  summarize(auc_raw = auc(x = index, y = enrichment),
            auc_std = auc_raw - background,
            enr_max = max(enrichment),
            enr_final = enrichment[length(enrichment)])
boxplot(d_permute_atac_auc$auc_std)
lapply(d_auc_atac$auc_std, points, col = "red", pch = 19)
quantile(d_permute_atac_auc$auc_std, probs = c(0.05, 0.95))

d_auc_atac$cell <- factor(d_auc_atac$cell, levels = c("cm", "ips", "lcl"),
                          labels = c("CM", "iPSC", "LCL"))
p_auc_atac <- p_auc %+% d_auc_atac %+%
  aes(x = reorder(cell, auc_std), fill = cell) +
  scale_fill_manual(values = cell_colors) +
  labs(x = "Cell type",
       title = "Quantification of enrichment of eQTLs in ATAC-seq peaks")

plot_grid(p_atac, p_auc_atac, labels = LETTERS[1:2])
ggsave("enrich-atac.png", width = 12, height = 6)
