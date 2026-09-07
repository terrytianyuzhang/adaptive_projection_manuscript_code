# File input and output -------------------------------------------------------

find_gap_covariance_source_file <- file.path(
  source_directory,
  "001_find_gap_covariance.R"
)

sparse_covariance_output_files <- c(
  feature_100 = file.path(
    data_directory,
    "raw",
    "random_sparse_covariance_p100.rds"
  ),
  feature_1000 = file.path(
    data_directory,
    "raw",
    "random_sparse_covariance_p1000.rds"
  )
)


# Source computational functions ---------------------------------------------

find_gap_covariance_source_file |> source()


# Covariance-generation settings ---------------------------------------------

feature_numbers <- c(100, 1000)
average_nonzero_off_diagonal_number_per_row_by_feature_number <- c(
  feature_100 = 8,
  feature_1000 = 8
)
covariance_perturbation_size <- 0.90
minimum_leading_eigenvalue_gap <- 0.20
maximum_search_attempt_number <- 500
covariance_seed <- 20260825


# Find and save one covariance matrix ----------------------------------------

all_sparse_covariance_diagnostics <- data.frame()

for (feature_number in feature_numbers) {
  average_nonzero_off_diagonal_number_per_row <- (
    unname(
      average_nonzero_off_diagonal_number_per_row_by_feature_number[
        paste0("feature_", feature_number)
      ]
    )
  )

  sparse_covariance_result <- find_sparse_covariance_matrix(
    feature_number = feature_number,
    average_nonzero_off_diagonal_number_per_row = (
      average_nonzero_off_diagonal_number_per_row
    ),
    covariance_perturbation_size = covariance_perturbation_size,
    minimum_leading_eigenvalue_gap = minimum_leading_eigenvalue_gap,
    maximum_search_attempt_number = maximum_search_attempt_number,
    covariance_seed = covariance_seed
  )

  sparse_covariance_output_file <- sparse_covariance_output_files[
    paste0("feature_", feature_number)
  ]

  dir.create(
    dirname(sparse_covariance_output_file),
    recursive = TRUE,
    showWarnings = FALSE
  )

  saveRDS(
    sparse_covariance_result,
    file = sparse_covariance_output_file
  )

  all_sparse_covariance_diagnostics <- rbind(
    all_sparse_covariance_diagnostics,
    sparse_covariance_result$diagnostics
  )
}

print(all_sparse_covariance_diagnostics)
