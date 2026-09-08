# File input ----------------------------------------------------------------

population_covariance_input_files <- c(
  dimension_100 = file.path(
    "..",
    "data",
    "raw",
    "001_random_sparse_covariance_p100.rds"
  ),
  dimension_1000 = file.path(
    "..",
    "data",
    "raw",
    "001_random_sparse_covariance_p1000.rds"
  )
)

# Sample size allocation ----------------------------------------------------

if (setting_id %in% c(1L, 6L, 11L, 16L)) {
  control_sample_size <- 50L
  treatment_sample_size <- 50L
} else if (setting_id %in% c(2L, 7L, 12L, 17L)) {
  control_sample_size <- 200L
  treatment_sample_size <- 200L
} else if (setting_id %in% c(3L, 8L, 13L, 18L)) {
  control_sample_size <- 400L
  treatment_sample_size <- 400L
} else if (setting_id %in% c(4L, 9L, 14L, 19L)) {
  control_sample_size <- 600L
  treatment_sample_size <- 600L
} else {
  control_sample_size <- 800L
  treatment_sample_size <- 800L
}


# Dimension allocation ------------------------------------------------------

if (setting_id %in% c(
  1L, 2L, 3L, 4L, 5L,
  11L, 12L, 13L, 14L, 15L
)) {
  dimension <- 100L
  population_covariance_input_file <- (
    population_covariance_input_files[["dimension_100"]]
  )
} else {
  dimension <- 1000L
  population_covariance_input_file <- (
    population_covariance_input_files[["dimension_1000"]]
  )
}

population_covariance_result <- readRDS(
  population_covariance_input_file
)
covariance_matrix <- as.matrix(
  population_covariance_result$covariance_matrix
)
leading_eigenvector <- (
  population_covariance_result$leading_eigenvectors[, 1L]
)
largest_loading_index <- which.max(abs(leading_eigenvector))
if (leading_eigenvector[largest_loading_index] < 0) {
  leading_eigenvector <- -leading_eigenvector
}


# Mean allocation -----------------------------------------------------------

active_feature_indices <- order(
  abs(leading_eigenvector),
  decreasing = TRUE
)[1:10]
projected_null_mean_difference <- rep(0, dimension)

for (pair_index in 1:5) {
  first_active_position <- 2L * pair_index - 1L
  second_active_position <- 2L * pair_index
  first_feature_index <- active_feature_indices[first_active_position]
  second_feature_index <- active_feature_indices[second_active_position]

  projected_null_mean_difference[first_feature_index] <- (
    leading_eigenvector[second_feature_index]
  )
  projected_null_mean_difference[second_feature_index] <- (
    -leading_eigenvector[first_feature_index]
  )
}
projected_null_mean_difference <- (
  projected_null_mean_difference /
    sqrt(sum(projected_null_mean_difference^2))
)

if (setting_id %in% c(
  1L, 2L, 3L, 4L, 5L,
  6L, 7L, 8L, 9L, 10L
)) {
  treatment_minus_control_mean <- projected_null_mean_difference
} else {
  truncated_leading_eigenvector <- rep(0, dimension)
  truncated_leading_eigenvector[active_feature_indices] <- (
    leading_eigenvector[active_feature_indices]
  )
  projected_alternative_addition <- (
    0.25 * truncated_leading_eigenvector /
      as.numeric(
        crossprod(
          leading_eigenvector,
          truncated_leading_eigenvector
        )
      )
  )
  treatment_minus_control_mean <- (
    projected_null_mean_difference +
      projected_alternative_addition
  )
}

control_mean <- -treatment_minus_control_mean / 2
treatment_mean <- treatment_minus_control_mean / 2


simulation_config <- list(
  setting_id = setting_id,
  control_mean = control_mean,
  treatment_mean = treatment_mean,
  control_sample_size = control_sample_size,
  treatment_sample_size = treatment_sample_size,
  covariance_matrix = covariance_matrix,
  dimension = dimension,
  distribution = "normal"
)
