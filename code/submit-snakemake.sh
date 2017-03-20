#!/bin/bash

snakemake \
    -kp \
    --ri \
    -j 500 \
    --cluster-config config-rcc.json \
    -c "sbatch \
        --mem={cluster.mem} \
        --nodes={cluster.n} \
        --tasks-per-node={cluster.tasks} \
        --partition=broadwl" \
    $*
