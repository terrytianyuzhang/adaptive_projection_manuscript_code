generate_gaussian_two_sample_from_cholesky <- function(
    sample_size_per_group,
    population_covariance_cholesky_factor,
    population_mean_difference) {
  feature_number <- ncol(population_covariance_cholesky_factor)

  group_1_sample <- (
    matrix(
      stats::rnorm(sample_size_per_group * feature_number),
      nrow = sample_size_per_group,
      ncol = feature_number
    ) %*%
      population_covariance_cholesky_factor
  )
  group_1_sample <- sweep(
    group_1_sample,
    MARGIN = 2,
    STATS = population_mean_difference / 2,
    FUN = "+"
  )

  group_2_sample <- (
    matrix(
      stats::rnorm(sample_size_per_group * feature_number),
      nrow = sample_size_per_group,
      ncol = feature_number
    ) %*%
      population_covariance_cholesky_factor
  )
  group_2_sample <- sweep(
    group_2_sample,
    MARGIN = 2,
    STATS = population_mean_difference / 2,
    FUN = "-"
  )

  list(
    group_1_sample = group_1_sample,
    group_2_sample = group_2_sample
  )
}


create_one_mean_comparison_result <- function(
    feature_number,
    sample_size_per_group,
    repeat_index,
    simulation_target,
    method,
    testing_result,
    significance_level) {
  data.frame(
    feature_number = feature_number,
    sample_size_per_group = sample_size_per_group,
    repeat_index = repeat_index,
    simulation_target = simulation_target,
    method = method,
    parameter_estimate = testing_result$parameter_estimate,
    standard_error = testing_result$standard_error,
    test_statistic = testing_result$test_statistic,
    p_value = testing_result$p_value,
    rejected = testing_result$p_value <= significance_level,
    stringsAsFactors = FALSE
  )
}


run_one_mean_comparison_simulation <- function(
    sample_size_per_group,
    repeat_index,
    simulation_seed,
    shared_simulation_input,
    method) {
  simulation_settings <- shared_simulation_input$simulation_settings
  feature_number <- shared_simulation_input$feature_number

  if (!method %in% names(simulation_settings$simulation_method_ids)) {
    stop("method is not listed in simulation_method_ids.")
  }

  set.seed(simulation_seed)
  projected_null_sample <- generate_gaussian_two_sample_from_cholesky(
    sample_size_per_group = sample_size_per_group,
    population_covariance_cholesky_factor = (
      shared_simulation_input$population_covariance_cholesky_factor
    ),
    population_mean_difference = (
      shared_simulation_input$projected_null_mean_difference
    )
  )
  projected_alternative_sample <- generate_gaussian_two_sample_from_cholesky(
    sample_size_per_group = sample_size_per_group,
    population_covariance_cholesky_factor = (
      shared_simulation_input$population_covariance_cholesky_factor
    ),
    population_mean_difference = (
      shared_simulation_input$projected_alternative_mean_difference
    )
  )

  if (method == "Oracle PC") {
    projected_null_testing_result <- oracle_pc_mean_comparison(
      group_1_sample = projected_null_sample$group_1_sample,
      group_2_sample = projected_null_sample$group_2_sample,
      population_leading_eigenvector = (
        shared_simulation_input$population_leading_eigenvector
      )
    )
    projected_alternative_testing_result <- oracle_pc_mean_comparison(
      group_1_sample = projected_alternative_sample$group_1_sample,
      group_2_sample = projected_alternative_sample$group_2_sample,
      population_leading_eigenvector = (
        shared_simulation_input$population_leading_eigenvector
      )
    )
  } else {
    projected_null_pc_results <- cross_fitted_pc_mean_comparison(
      group_1_sample = projected_null_sample$group_1_sample,
      group_2_sample = projected_null_sample$group_2_sample,
      fold_number = simulation_settings$cross_fitting_fold_number,
      threshold_multiplier = simulation_settings$threshold_multiplier,
      mean_difference_threshold_multiplier = (
        simulation_settings$mean_difference_threshold_multiplier
      ),
      cross_fitting_seed = simulation_seed + 1L
    )
    projected_alternative_pc_results <- cross_fitted_pc_mean_comparison(
      group_1_sample = projected_alternative_sample$group_1_sample,
      group_2_sample = projected_alternative_sample$group_2_sample,
      fold_number = simulation_settings$cross_fitting_fold_number,
      threshold_multiplier = simulation_settings$threshold_multiplier,
      mean_difference_threshold_multiplier = (
        simulation_settings$mean_difference_threshold_multiplier
      ),
      cross_fitting_seed = simulation_seed + 2L
    )
    result_name <- if (method == "Debiased PC") "debiased" else "plug_in"
    projected_null_testing_result <- projected_null_pc_results[[result_name]]
    projected_alternative_testing_result <- (
      projected_alternative_pc_results[[result_name]]
    )
  }

  testing_result_specifications <- list(
    list("Type I error", method, projected_null_testing_result),
    list("Power", method, projected_alternative_testing_result)
  )
  simulation_results <- vector(
    "list",
    length(testing_result_specifications)
  )

  for (
    result_number in seq_along(testing_result_specifications)
  ) {
    testing_result_specification <- (
      testing_result_specifications[[result_number]]
    )
    simulation_results[[result_number]] <- create_one_mean_comparison_result(
      feature_number = feature_number,
      sample_size_per_group = sample_size_per_group,
      repeat_index = repeat_index,
      simulation_target = testing_result_specification[[1]],
      method = testing_result_specification[[2]],
      testing_result = testing_result_specification[[3]],
      significance_level = simulation_settings$significance_level
    )
  }

  do.call(rbind, simulation_results)
}


create_mean_comparison_task_result_file <- function(
    simulation_task_result_directory,
    task_specification,
    simulation_method_ids) {
  method_id <- unname(
    simulation_method_ids[[as.character(task_specification$method)]]
  )

  file.path(
    simulation_task_result_directory,
    method_id,
    as.character(task_specification$setting),
    paste0(
      "012_mean_comparison_",
      as.character(task_specification$unique_id),
      ".rds"
    )
  )
}


validate_mean_comparison_task_output <- function(
    task_output,
    expected_task_specification,
    expected_setting_specification,
    expected_simulation_settings) {
  if (!is.list(task_output)) {
    return(FALSE)
  }

  if (
    !all(
      c(
        "task_specification",
        "simulation_results",
        "simulation_settings"
      ) %in% names(task_output)
    )
  ) {
    return(FALSE)
  }

  if (!identical(
    task_output$simulation_settings,
    expected_simulation_settings
  )) {
    return(FALSE)
  }

  saved_task_specification <- task_output$task_specification
  if (
    nrow(saved_task_specification) != 1L ||
      !identical(
        names(saved_task_specification),
        names(expected_task_specification)
      ) ||
      !identical(
        as.character(saved_task_specification$unique_id),
        as.character(expected_task_specification$unique_id)
      ) ||
      !identical(
        as.character(saved_task_specification$method),
        as.character(expected_task_specification$method)
      ) ||
      !identical(
        as.character(saved_task_specification$setting),
        as.character(expected_task_specification$setting)
      ) ||
      !identical(
        as.integer(saved_task_specification$batch_number),
        as.integer(expected_task_specification$batch_number)
      )
  ) {
    return(FALSE)
  }

  simulation_results <- task_output$simulation_results
  required_result_names <- c(
    "feature_number",
    "sample_size_per_group",
    "repeat_index",
    "simulation_target",
    "method",
    "parameter_estimate",
    "standard_error",
    "test_statistic",
    "p_value",
    "rejected"
  )
  if (!all(required_result_names %in% names(simulation_results))) {
    return(FALSE)
  }

  expected_repeat_indices <- seq.int(
    (expected_task_specification$batch_number - 1L) *
      expected_simulation_settings$repeat_number_per_batch + 1L,
    expected_task_specification$batch_number *
      expected_simulation_settings$repeat_number_per_batch
  )
  expected_simulation_targets <- c("Type I error", "Power")

  if (
    nrow(simulation_results) != length(expected_repeat_indices) * 2L ||
      !identical(
        sort(unique(simulation_results$repeat_index)),
        as.integer(expected_repeat_indices)
      ) ||
      !identical(
        unique(simulation_results$method),
        as.character(expected_task_specification$method)
      ) ||
      !setequal(
        unique(simulation_results$simulation_target),
        expected_simulation_targets
      ) ||
      any(
        simulation_results$feature_number !=
          expected_setting_specification$feature_number
      ) ||
      any(
        simulation_results$sample_size_per_group !=
          expected_setting_specification$sample_size_per_group
      )
  ) {
    return(FALSE)
  }

  result_combination_counts <- table(
    simulation_results$repeat_index,
    simulation_results$simulation_target,
    simulation_results$method
  )

  all(result_combination_counts == 1L)
}
