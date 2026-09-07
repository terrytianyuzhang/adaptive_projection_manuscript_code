# File input and output -------------------------------------------------------

mean_comparison_source_file <- file.path(
  source_directory,
  "002_threshold_sparse_covariance.R"
)

population_covariance_input_file <- file.path(
  data_directory,
  "raw",
  "001_random_sparse_covariance_p100.rds"
)

mean_comparison_simulation_output_file <- file.path(
  data_directory,
  "intermediate",
  "003_mean_comparison_simulation_p100.rds"
)

mean_comparison_simulation_checkpoint_file <- file.path(
  data_directory,
  "intermediate",
  "003_mean_comparison_simulation_checkpoint_p100.rds"
)


# Source computational functions ---------------------------------------------

mean_comparison_source_file |> source()


# Simulation settings --------------------------------------------------------

sample_sizes_per_group <- c(50, 200, 400, 600, 800)
simulation_repeat_number <- 150
cross_fitting_fold_number <- 5
threshold_multiplier <- 1.25
mean_difference_threshold_multiplier <- 1.25
projected_null_mean_difference_size <- 1.00
projected_alternative_mean_difference_size <- 0.25
projected_alternative_support_size <- 10
significance_level <- 0.05
mean_comparison_simulation_seed <- 20260829

simulation_settings <- list(
  sample_sizes_per_group = sample_sizes_per_group,
  simulation_repeat_number = simulation_repeat_number,
  cross_fitting_fold_number = cross_fitting_fold_number,
  threshold_multiplier = threshold_multiplier,
  mean_difference_threshold_multiplier = (
    mean_difference_threshold_multiplier
  ),
  projected_null_mean_difference_size = (
    projected_null_mean_difference_size
  ),
  projected_alternative_mean_difference_size = (
    projected_alternative_mean_difference_size
  ),
  projected_alternative_support_size = (
    projected_alternative_support_size
  ),
  significance_level = significance_level,
  mean_comparison_simulation_seed = mean_comparison_simulation_seed
)


# Population quantities ------------------------------------------------------

population_covariance_result <- readRDS(
  population_covariance_input_file
)
population_covariance_matrix <- as.matrix(
  population_covariance_result$covariance_matrix
)
population_eigendecomposition <- eigen(
  population_covariance_matrix,
  symmetric = TRUE
)
population_leading_eigenvector <- (
  population_eigendecomposition$vectors[, 1]
)
global_null_mean_difference <- rep(
  0,
  nrow(population_covariance_matrix)
)
largest_leading_eigenvector_loading_indices <- order(
  abs(population_leading_eigenvector),
  decreasing = TRUE
)[1:2]
projected_null_direction <- rep(
  0,
  length(population_leading_eigenvector)
)
projected_null_direction[
  largest_leading_eigenvector_loading_indices[1]
] <- population_leading_eigenvector[
  largest_leading_eigenvector_loading_indices[2]
]
projected_null_direction[
  largest_leading_eigenvector_loading_indices[2]
] <- -population_leading_eigenvector[
  largest_leading_eigenvector_loading_indices[1]
]
projected_null_direction <- (
  projected_null_direction /
    sqrt(sum(projected_null_direction^2))
)
projected_null_mean_difference <- (
  projected_null_mean_difference_size * projected_null_direction
)
projected_alternative_loading_indices <- order(
  abs(population_leading_eigenvector),
  decreasing = TRUE
)[seq_len(projected_alternative_support_size)]
projected_alternative_direction <- rep(
  0,
  length(population_leading_eigenvector)
)
projected_alternative_direction[
  projected_alternative_loading_indices
] <- population_leading_eigenvector[
  projected_alternative_loading_indices
]
projected_alternative_direction <- (
  projected_alternative_direction /
    as.numeric(
      crossprod(
        population_leading_eigenvector,
        projected_alternative_direction
      )
    )
)
projected_alternative_mean_difference <- (
  projected_null_mean_difference +
    projected_alternative_mean_difference_size *
      projected_alternative_direction
)


# Result helpers --------------------------------------------------------------

create_one_simulation_result <- function(
    sample_size_per_group,
    repeat_index,
    test_family,
    simulation_target,
    method,
    testing_result) {
  data.frame(
    feature_number = nrow(population_covariance_matrix),
    sample_size_per_group = sample_size_per_group,
    repeat_index = repeat_index,
    test_family = test_family,
    simulation_target = simulation_target,
    method = method,
    parameter_estimate = testing_result$parameter_estimate,
    standard_error = testing_result$standard_error,
    test_statistic = testing_result$test_statistic,
    p_value = testing_result$p_value,
    rejected = testing_result$p_value <= significance_level
  )
}


# Run simulation --------------------------------------------------------------

simulation_start_time <- proc.time()[["elapsed"]]
mean_comparison_simulation_results <- list()
simulation_result_index <- 1

if (file.exists(mean_comparison_simulation_checkpoint_file)) {
  existing_checkpoint <- readRDS(
    mean_comparison_simulation_checkpoint_file
  )

  existing_settings_for_comparison <- existing_checkpoint$settings
  current_settings_for_comparison <- simulation_settings
  existing_settings_for_comparison$simulation_repeat_number <- NULL
  current_settings_for_comparison$simulation_repeat_number <- NULL
  existing_settings_for_comparison$sample_sizes_per_group <- NULL
  current_settings_for_comparison$sample_sizes_per_group <- NULL

  if (
    identical(
      existing_settings_for_comparison,
      current_settings_for_comparison
    )
  ) {
    existing_checkpoint$simulation_results <- (
      existing_checkpoint$simulation_results[
        existing_checkpoint$simulation_results$repeat_index <=
          simulation_repeat_number &
          existing_checkpoint$simulation_results$sample_size_per_group %in%
          sample_sizes_per_group,
        ,
        drop = FALSE
      ]
    )
    mean_comparison_simulation_results <- split(
      existing_checkpoint$simulation_results,
      seq_len(nrow(existing_checkpoint$simulation_results))
    )
    simulation_result_index <- (
      length(mean_comparison_simulation_results) + 1
    )
    message(
      "Restored ",
      nrow(existing_checkpoint$simulation_results),
      " checkpointed result rows."
    )
  }
}

for (sample_size_per_group in sample_sizes_per_group) {
  for (repeat_index in seq_len(simulation_repeat_number)) {
    existing_result_indicator <- vapply(
      mean_comparison_simulation_results,
      function(existing_result) {
        existing_result$sample_size_per_group == sample_size_per_group &&
          existing_result$repeat_index == repeat_index
      },
      logical(1)
    )

    if (sum(existing_result_indicator) == 10) {
      next
    }

    set.seed(
      mean_comparison_simulation_seed +
        10000 * sample_size_per_group +
        repeat_index
    )

    global_null_sample <- generate_gaussian_two_sample(
      sample_size_per_group = sample_size_per_group,
      population_covariance_matrix = population_covariance_matrix,
      population_mean_difference = global_null_mean_difference
    )
    projected_null_sample <- generate_gaussian_two_sample(
      sample_size_per_group = sample_size_per_group,
      population_covariance_matrix = population_covariance_matrix,
      population_mean_difference = projected_null_mean_difference
    )
    projected_alternative_sample <- generate_gaussian_two_sample(
      sample_size_per_group = sample_size_per_group,
      population_covariance_matrix = population_covariance_matrix,
      population_mean_difference = projected_alternative_mean_difference
    )

    global_null_pc_results <- cross_fitted_pc_mean_comparison(
      group_1_sample = global_null_sample$group_1_sample,
      group_2_sample = global_null_sample$group_2_sample,
      fold_number = cross_fitting_fold_number,
      threshold_multiplier = threshold_multiplier,
      mean_difference_threshold_multiplier = (
        mean_difference_threshold_multiplier
      )
    )
    global_null_anchored_result <- anchored_lasso_mean_comparison(
      group_1_sample = global_null_sample$group_1_sample,
      group_2_sample = global_null_sample$group_2_sample,
      fold_number = cross_fitting_fold_number,
      threshold_multiplier = threshold_multiplier
    )
    projected_null_pc_results <- cross_fitted_pc_mean_comparison(
      group_1_sample = projected_null_sample$group_1_sample,
      group_2_sample = projected_null_sample$group_2_sample,
      fold_number = cross_fitting_fold_number,
      threshold_multiplier = threshold_multiplier,
      mean_difference_threshold_multiplier = (
        mean_difference_threshold_multiplier
      )
    )
    projected_null_anchored_result <- anchored_lasso_mean_comparison(
      group_1_sample = projected_null_sample$group_1_sample,
      group_2_sample = projected_null_sample$group_2_sample,
      fold_number = cross_fitting_fold_number,
      threshold_multiplier = threshold_multiplier
    )
    projected_alternative_pc_results <- cross_fitted_pc_mean_comparison(
      group_1_sample = projected_alternative_sample$group_1_sample,
      group_2_sample = projected_alternative_sample$group_2_sample,
      fold_number = cross_fitting_fold_number,
      threshold_multiplier = threshold_multiplier,
      mean_difference_threshold_multiplier = (
        mean_difference_threshold_multiplier
      )
    )
    projected_null_oracle_result <- oracle_pc_mean_comparison(
      group_1_sample = projected_null_sample$group_1_sample,
      group_2_sample = projected_null_sample$group_2_sample,
      population_leading_eigenvector = population_leading_eigenvector
    )
    projected_alternative_oracle_result <- oracle_pc_mean_comparison(
      group_1_sample = projected_alternative_sample$group_1_sample,
      group_2_sample = projected_alternative_sample$group_2_sample,
      population_leading_eigenvector = population_leading_eigenvector
    )

    testing_result_specifications <- list(
      list("Global null", "Type I error", "Plug-in PC", global_null_pc_results$plug_in),
      list("Global null", "Type I error", "Anchored Lasso", global_null_anchored_result),
      list("Global null", "Power", "Plug-in PC", projected_null_pc_results$plug_in),
      list("Global null", "Power", "Anchored Lasso", projected_null_anchored_result),
      list("Projected null", "Type I error", "Oracle PC", projected_null_oracle_result),
      list("Projected null", "Type I error", "Plug-in PC", projected_null_pc_results$plug_in),
      list("Projected null", "Type I error", "Debiased PC", projected_null_pc_results$debiased),
      list("Projected null", "Power", "Oracle PC", projected_alternative_oracle_result),
      list("Projected null", "Power", "Plug-in PC", projected_alternative_pc_results$plug_in),
      list("Projected null", "Power", "Debiased PC", projected_alternative_pc_results$debiased)
    )

    for (testing_result_specification in testing_result_specifications) {
      mean_comparison_simulation_results[[simulation_result_index]] <- (
        create_one_simulation_result(
          sample_size_per_group = sample_size_per_group,
          repeat_index = repeat_index,
          test_family = testing_result_specification[[1]],
          simulation_target = testing_result_specification[[2]],
          method = testing_result_specification[[3]],
          testing_result = testing_result_specification[[4]]
        )
      )
      simulation_result_index <- simulation_result_index + 1
    }

    if (
      repeat_index %% 10 == 0 ||
        repeat_index == simulation_repeat_number
    ) {
      checkpointed_simulation_results <- do.call(
        rbind,
        mean_comparison_simulation_results
      )
      saveRDS(
        list(
          simulation_results = checkpointed_simulation_results,
          settings = simulation_settings
        ),
        file = mean_comparison_simulation_checkpoint_file
      )
    }
  }

  message(
    "Completed sample size per group: ",
    sample_size_per_group
  )
}

mean_comparison_simulation_results <- do.call(
  rbind,
  mean_comparison_simulation_results
)
simulation_elapsed_seconds <- (
  proc.time()[["elapsed"]] - simulation_start_time
)

mean_comparison_simulation_output <- list(
  simulation_results = mean_comparison_simulation_results,
  settings = simulation_settings,
  current_session_elapsed_seconds = simulation_elapsed_seconds
)

dir.create(
  dirname(mean_comparison_simulation_output_file),
  recursive = TRUE,
  showWarnings = FALSE
)
saveRDS(
  mean_comparison_simulation_output,
  file = mean_comparison_simulation_output_file
)

message(
  "Current session elapsed minutes: ",
  round(simulation_elapsed_seconds / 60, 2)
)
