#!/bin/bash

set -euo pipefail

# Run this script from the Biostatistics_revision/code directory.
number_of_workers="${1:-8}"

Rscript "011_prepare_mean_comparison_tasks.R"

# Script 013 calls script 012 once for every simulation task.
bash "013_run_mean_comparison_parallel.sh" "$number_of_workers"

Rscript "014_collect_mean_comparison_results.R"

