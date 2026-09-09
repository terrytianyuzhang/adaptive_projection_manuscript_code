#!/bin/bash

set -euo pipefail

# Run this script from the Biostatistics_revision/code directory.
number_of_workers="${1:-8}"

Rscript "031_prepare_mean_comparison_tasks.R"

# Script 033 calls script 032 once for every simulation task.
bash "033_run_mean_comparison_parallel.sh" "$number_of_workers"

Rscript "034_collect_mean_comparison_results.R"

