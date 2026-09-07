library(here)

work_directory <- file.path(
  here(),
  "code_submitted",
  "Biostatistics_revision",
  "code"
)
source_directory <- file.path(
  here(),
  "code_submitted",
  "Biostatistics_revision",
  "src"
)
data_directory <- file.path(
  here(),
  "code_submitted",
  "Biostatistics_revision",
  "data"
)

file.path(work_directory, "001_find_sparse_covariance_matrix.R") |> source()
file.path(work_directory, "002_validate_sparse_covariance_estimator.R") |> source()
file.path(work_directory, "003_run_mean_comparison_simulation.R") |> source()
file.path(work_directory, "004_plot_mean_comparison_simulation.R") |> source()
