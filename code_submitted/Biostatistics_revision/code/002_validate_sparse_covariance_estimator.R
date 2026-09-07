# File input and output -------------------------------------------------------

threshold_covariance_source_file <- file.path(
  source_directory,
  "002_threshold_sparse_covariance.R"
)

population_covariance_input_files <- c(
  feature_100 = file.path(
    data_directory,
    "raw",
    "001_random_sparse_covariance_p100.rds"
  ),
  feature_1000 = file.path(
    data_directory,
    "raw",
    "001_random_sparse_covariance_p1000.rds"
  )
)

operator_norm_diagnostic_output_file <- file.path(
  data_directory,
  "intermediate",
  "002_threshold_covariance_operator_norm_diagnostics.rds"
)


# Source computational functions ---------------------------------------------

threshold_covariance_source_file |> source()


# Diagnostic settings --------------------------------------------------------

sample_sizes_per_group <- c(50, 100, 200, 400, 800)
covariance_estimation_repeat_number <- 5
threshold_multipliers <- c(1.0, 1.25, 1.5)
covariance_estimation_seed <- 20260828


# Evaluate covariance estimation ---------------------------------------------

set.seed(covariance_estimation_seed)
operator_norm_diagnostics <- data.frame()

for (population_covariance_input_file in population_covariance_input_files) {
  population_covariance_result <- readRDS(
    population_covariance_input_file
  )
  population_covariance_matrix <- as.matrix(
    population_covariance_result$covariance_matrix
  )
  feature_number <- nrow(population_covariance_matrix)
  population_covariance_cholesky_factor <- chol(
    population_covariance_matrix
  )
  maximum_sample_size_per_group <- max(sample_sizes_per_group)

  for (
    repeat_index in seq_len(covariance_estimation_repeat_number)
  ) {
    maximum_group_1_sample <- (
      matrix(
        stats::rnorm(maximum_sample_size_per_group * feature_number),
        nrow = maximum_sample_size_per_group,
        ncol = feature_number
      ) %*%
        population_covariance_cholesky_factor
    )
    maximum_group_2_sample <- (
      matrix(
        stats::rnorm(maximum_sample_size_per_group * feature_number),
        nrow = maximum_sample_size_per_group,
        ncol = feature_number
      ) %*%
        population_covariance_cholesky_factor
    )

    for (sample_size_per_group in sample_sizes_per_group) {
      group_1_sample <- maximum_group_1_sample[
        seq_len(sample_size_per_group),
        ,
        drop = FALSE
      ]
      group_2_sample <- maximum_group_2_sample[
        seq_len(sample_size_per_group),
        ,
        drop = FALSE
      ]

      for (threshold_multiplier in threshold_multipliers) {
        covariance_estimation_result <- estimate_thresholded_common_covariance(
          group_1_sample = group_1_sample,
          group_2_sample = group_2_sample,
          threshold_multiplier = threshold_multiplier
        )
        covariance_estimation_error <- (
          covariance_estimation_result$estimated_covariance_matrix -
            population_covariance_result$covariance_matrix
        )
        covariance_operator_norm_error <- calculate_symmetric_operator_norm(
          covariance_estimation_error
        )
        estimated_leading_eigenpairs <- calculate_leading_symmetric_eigenpairs(
          covariance_estimation_result$estimated_covariance_matrix,
          eigenpair_number = 2
        )
        leading_eigenvector_alignment <- abs(
          crossprod(
            population_covariance_result$leading_eigenvectors[, 1],
            estimated_leading_eigenpairs$vectors[, 1]
          )
        )

        one_diagnostic_result <- data.frame(
          feature_number = feature_number,
          sample_size_per_group = sample_size_per_group,
          total_sample_size = 2 * sample_size_per_group,
          repeat_index = repeat_index,
          threshold_multiplier = threshold_multiplier,
          covariance_operator_norm_error = covariance_operator_norm_error,
          estimated_leading_eigenvalue_gap = (
            estimated_leading_eigenpairs$values[1] -
              estimated_leading_eigenpairs$values[2]
          ),
          leading_eigenvector_alignment = leading_eigenvector_alignment,
          retained_off_diagonal_entry_number = (
            covariance_estimation_result$diagnostics$
              retained_off_diagonal_entry_number
          ),
          positive_semidefinite_shift = (
            covariance_estimation_result$diagnostics$
              positive_semidefinite_shift
          )
        )

        operator_norm_diagnostics <- rbind(
          operator_norm_diagnostics,
          one_diagnostic_result
        )
      }
    }
  }
}

dir.create(
  dirname(operator_norm_diagnostic_output_file),
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  operator_norm_diagnostics,
  file = operator_norm_diagnostic_output_file
)

operator_norm_summary <- aggregate(
  covariance_operator_norm_error ~
    feature_number + sample_size_per_group + threshold_multiplier,
  data = operator_norm_diagnostics,
  FUN = mean
)

print(operator_norm_summary)
