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
