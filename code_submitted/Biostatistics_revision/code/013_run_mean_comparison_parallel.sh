#!/bin/bash

set -euo pipefail

# Run this script from the Biostatistics_revision/code directory.
number_of_workers="${1:-8}"

simulation_task_manifest_file="../data/intermediate/011_mean_comparison_task_manifest.csv"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

tail -n +2 "$simulation_task_manifest_file" |
  cut -d "," -f 1 |
  xargs \
    -P "$number_of_workers" \
    -n 1 \
    Rscript "012_run_mean_comparison_task.R"
