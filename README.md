# cardioqtl [![Anaconda-Server Badge](https://anaconda.org/jdblischak/cardioqtl/badges/version.svg)](https://anaconda.org/jdblischak/cardioqtl)

## Installation

I used [conda][] and [bioconda][] to manage software dependencies. To replicate
the computing environment, you will need to complete the following 4 steps. Note
that this is only guaranteed to work on a Linux-64 based architecture.

1. Download and install Miniconda ([instructions](https://conda.io/miniconda.html))
2. Download environment file. If you cloned the Git repo, you already have the environment file.
    - `conda install anaconda-client; anaconda download jdblischak/cardioqtl`
3. Create conda environment:
    - `conda env create -n cardioqtl --file environment.yaml`
4. Activate conda environment:
    - To activate: `source activate cardioqtl`
    - To deactivate: `source deactivate cardioqtl`

[conda]: https://conda.io/docs/
[bioconda]: https://bioconda.github.io

## Code

The code is freely available for reuse with attribution via the [MIT
license][mit].

* `Snakefile` - Implements the analysis pipeline

* `submit-snakemake.sh` - Submits individual jobs produced in
  `Snakefile` to [Slurm][]. If your cluster uses a different job
  scheduler, you'll need to edit this file and `cluster.json`.

* `scripts/` - R scripts called by `Snakefile`

* `scratch/` - Exploratory analyses written in R Markdown

[mit]: https://choosealicense.com/licenses/mit/
[Slurm]: https://slurm.schedmd.com/overview.html

## Data

* `data/counts-subread.txt` - Gene counts after mapping to GRCh37 with
  Subjunc and summing counts per gene with featureCounts (Subread
  1.5.0p3). Includes all genes in Ensembl release 75 (i.e. protein
  coding plus all other biotypes; see `scripts/create-exons.R`).

* `data/counts-clean.txt` - Gene counts after removing samples 26302,
  110232, and 160001 and removing genes with log2 cpm less than 0 (see
  `scripts/clean-counts.R`).

* `data/counts-normalized.txt` - Gene counts after normalizing to
  N(0,1) within each sample follwed by normalizing to N(0,1) within
  each gene. Used the R function `qqnorm` (see
  `scripts/normalize-counts.R`).
