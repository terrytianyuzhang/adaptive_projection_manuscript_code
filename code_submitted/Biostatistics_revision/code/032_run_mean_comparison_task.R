# Run as: Rscript 032_run_mean_comparison_task.R task_id

# File input and output ------------------------------------------------------

simulation_task_table_input_file <- file.path(
  "..",
  "data",
  "intermediate",
  "031_mean_comparison_task_table.csv"
)
simulation_config_source_file <- file.path(
  "..",
  "src",
  "031_mean_comparison_simulation_config.R"
)
mean_comparison_source_file <- file.path(
  "..",
  "src",
  "002_threshold_sparse_covariance.R"
)
simulation_task_result_directory <- file.path(
  "..",
  "data",
  "intermediate",
  "032_mean_comparison_task_result"
)


# Simulation constants ------------------------------------------------------

repetition_number_per_batch <- 10L
simulation_seed_offset <- 20260829L
cross_fitting_fold_number <- 5L
threshold_multiplier <- 1.25
mean_difference_threshold_multiplier <- 1.25


# Read task -----------------------------------------------------------------

task_id <- as.integer(commandArgs(trailingOnly = TRUE)[1L])
simulation_task_table <- read.csv(
  simulation_task_table_input_file,
  stringsAsFactors = FALSE
)
task_specification <- simulation_task_table[
  simulation_task_table$id == task_id,
  ,
  drop = FALSE
]

method <- task_specification$method
setting_id <- task_specification$setting_id
batch_number <- task_specification$batch_number

source(simulation_config_source_file)
source(mean_comparison_source_file)


# Generate and analyze repetitions -----------------------------------------

population_covariance_cholesky_factor <- chol(
  simulation_config$covariance_matrix
)
simulation_results <- vector(
  "list",
  repetition_number_per_batch
)

for (repetition_number in 1:repetition_number_per_batch) {
  repeat_index <- (
    (batch_number - 1L) * repetition_number_per_batch +
      repetition_number
  )
  simulation_seed <- (
    simulation_seed_offset +
      setting_id * 10000L +
      repeat_index
  )
  set.seed(simulation_seed)

  treatment_sample <- (
    matrix(
      stats::rnorm(
        simulation_config$treatment_sample_size *
          simulation_config$dimension
      ),
      nrow = simulation_config$treatment_sample_size,
      ncol = simulation_config$dimension
    ) %*%
      population_covariance_cholesky_factor
  )
  treatment_sample <- sweep(
    treatment_sample,
    MARGIN = 2,
    STATS = simulation_config$treatment_mean,
    FUN = "+"
  )

  control_sample <- (
    matrix(
      stats::rnorm(
        simulation_config$control_sample_size *
          simulation_config$dimension
      ),
      nrow = simulation_config$control_sample_size,
      ncol = simulation_config$dimension
    ) %*%
      population_covariance_cholesky_factor
  )
  control_sample <- sweep(
    control_sample,
    MARGIN = 2,
    STATS = simulation_config$control_mean,
    FUN = "+"
  )

  if (method == "debiased_pc") {
    cross_fitted_result <- cross_fitted_pc_mean_comparison(
      group_1_sample = treatment_sample,
      group_2_sample = control_sample,
      fold_number = cross_fitting_fold_number,
      threshold_multiplier = threshold_multiplier,
      mean_difference_threshold_multiplier = (
        mean_difference_threshold_multiplier
      ),
      cross_fitting_seed = simulation_seed + 1L
    )
    testing_result <- cross_fitted_result$debiased
  } else if (method == "oracle_pc") {
    testing_result <- oracle_pc_mean_comparison(
      group_1_sample = treatment_sample,
      group_2_sample = control_sample,
      population_leading_eigenvector = leading_eigenvector
    )
  } else if (method == "plug_in_pc") {
    cross_fitted_result <- cross_fitted_pc_mean_comparison(
      group_1_sample = treatment_sample,
      group_2_sample = control_sample,
      fold_number = cross_fitting_fold_number,
      threshold_multiplier = threshold_multiplier,
      mean_difference_threshold_multiplier = (
        mean_difference_threshold_multiplier
      ),
      cross_fitting_seed = simulation_seed + 1L
    )
    testing_result <- cross_fitted_result$plug_in
  }

  simulation_results[[repetition_number]] <- list(
    repeat_index = repeat_index,
    simulation_seed = simulation_seed,
    setting_id = setting_id,
    method = method,
    parameter_estimate = testing_result$parameter_estimate,
    standard_error = testing_result$standard_error,
    test_statistic = testing_result$test_statistic,
    p_value = testing_result$p_value
  )
}


# Save results --------------------------------------------------------------

simulation_task_result_file <- file.path(
  simulation_task_result_directory,
  paste0(
    "032_mean_comparison_task_",
    task_id,
    ".rds"
  )
)
dir.create(
  simulation_task_result_directory,
  recursive = TRUE,
  showWarnings = FALSE
)
saveRDS(
  simulation_results,
  file = simulation_task_result_file
)
