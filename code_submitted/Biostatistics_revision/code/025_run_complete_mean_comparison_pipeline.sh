#!/bin/bash

set -euo pipefail

# Run this script from the Biostatistics_revision/code directory.
number_of_workers="${1:-8}"

Rscript "021_prepare_mean_comparison_tasks.R"

# Script 023 calls script 022 once for every simulation task.
bash "023_run_mean_comparison_parallel.sh" "$number_of_workers"

Rscript "024_collect_mean_comparison_results.R"
