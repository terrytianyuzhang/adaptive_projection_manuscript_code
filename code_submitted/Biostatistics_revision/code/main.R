# Required source scripts:
# src/001_find_gap_covariance.R
# src/002_threshold_sparse_covariance.R
# src/031_mean_comparison_simulation_config.R

# find proper covariance matrices for the simulation
file.path(".", "001_find_sparse_covariance_matrix.R")

# Run the simulation with 16 cores.
# bash 035_run_complete_mean_comparison_pipeline.sh 16

# Create type-I and power plot using the output from 035_
file.path(".", "036_plot_mean_comparison_results.R") 

# Create publciation plot
file.path(".", "037_publication_plot.R")



