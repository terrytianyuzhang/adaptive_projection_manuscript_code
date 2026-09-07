validate_covariance_estimation_sample <- function(sample, sample_name) {
  if (!is.matrix(sample) && !is.data.frame(sample)) {
    stop(paste(sample_name, "must be a numeric matrix or data frame."))
  }

  sample <- as.matrix(sample)

  if (!is.numeric(sample)) {
    stop(paste(sample_name, "must contain only numeric values."))
  }

  if (nrow(sample) < 2) {
    stop(paste(sample_name, "must contain at least two observations."))
  }

  if (ncol(sample) < 2) {
    stop(paste(sample_name, "must contain at least two features."))
  }

  if (any(!is.finite(sample))) {
    stop(paste(sample_name, "must contain only finite values."))
  }

  sample
}


validate_threshold_covariance_numeric_scalar <- function(value, value_name) {
  if (
    !is.numeric(value) ||
      length(value) != 1 ||
      is.na(value) ||
      !is.finite(value)
  ) {
    stop(paste(value_name, "must be one finite numeric value."))
  }
}


calculate_smallest_symmetric_eigenvalue <- function(symmetric_matrix) {
  matrix_dimension <- nrow(symmetric_matrix)

  if (matrix_dimension <= 3) {
    return(
      min(
        eigen(
          as.matrix(symmetric_matrix),
          symmetric = TRUE,
          only.values = TRUE
        )$values
      )
    )
  }

  symmetric_matrix_for_eigensolver <- as(
    symmetric_matrix,
    "generalMatrix"
  )

  smallest_eigenvalue_result <- RSpectra::eigs_sym(
    symmetric_matrix_for_eigensolver,
    k = 1,
    which = "SA",
    opts = list(retvec = FALSE)
  )

  if (is.list(smallest_eigenvalue_result)) {
    return(smallest_eigenvalue_result$values[1])
  }

  smallest_eigenvalue_result[1]
}


calculate_symmetric_operator_norm <- function(symmetric_matrix) {
  matrix_dimension <- nrow(symmetric_matrix)

  if (matrix_dimension <= 3) {
    return(
      max(
        abs(
          eigen(
            as.matrix(symmetric_matrix),
            symmetric = TRUE,
            only.values = TRUE
          )$values
        )
      )
    )
  }

  symmetric_matrix_for_eigensolver <- as(
    symmetric_matrix,
    "generalMatrix"
  )

  largest_magnitude_eigenvalue_result <- RSpectra::eigs_sym(
    symmetric_matrix_for_eigensolver,
    k = 1,
    which = "LM",
    opts = list(retvec = FALSE)
  )

  if (is.list(largest_magnitude_eigenvalue_result)) {
    return(abs(largest_magnitude_eigenvalue_result$values[1]))
  }

  abs(largest_magnitude_eigenvalue_result[1])
}


calculate_leading_symmetric_eigenpairs <- function(
    symmetric_matrix,
    eigenpair_number = 2) {
  matrix_dimension <- nrow(symmetric_matrix)

  if (eigenpair_number < 1 || eigenpair_number >= matrix_dimension) {
    stop(
      paste(
        "eigenpair_number must be at least 1 and smaller than the",
        "matrix dimension."
      )
    )
  }

  if (matrix_dimension <= 3) {
    eigen_result <- eigen(
      as.matrix(symmetric_matrix),
      symmetric = TRUE
    )

    return(list(
      values = eigen_result$values[seq_len(eigenpair_number)],
      vectors = eigen_result$vectors[
        ,
        seq_len(eigenpair_number),
        drop = FALSE
      ]
    ))
  }

  symmetric_matrix_for_eigensolver <- as(
    symmetric_matrix,
    "generalMatrix"
  )
  eigen_result <- RSpectra::eigs_sym(
    symmetric_matrix_for_eigensolver,
    k = eigenpair_number,
    which = "LA"
  )
  decreasing_eigenvalue_order <- order(
    eigen_result$values,
    decreasing = TRUE
  )

  list(
    values = eigen_result$values[decreasing_eigenvalue_order],
    vectors = eigen_result$vectors[
      ,
      decreasing_eigenvalue_order,
      drop = FALSE
    ]
  )
}


estimate_thresholded_common_covariance <- function(
    group_1_sample,
    group_2_sample,
    threshold_multiplier = 1.25,
    minimum_estimated_eigenvalue = 1e-6) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The Matrix package is required.")
  }

  if (!requireNamespace("RSpectra", quietly = TRUE)) {
    stop("The RSpectra package is required.")
  }

  group_1_sample <- validate_covariance_estimation_sample(
    group_1_sample,
    "group_1_sample"
  )
  group_2_sample <- validate_covariance_estimation_sample(
    group_2_sample,
    "group_2_sample"
  )

  if (ncol(group_1_sample) != ncol(group_2_sample)) {
    stop("group_1_sample and group_2_sample must have the same features.")
  }

  validate_threshold_covariance_numeric_scalar(
    threshold_multiplier,
    "threshold_multiplier"
  )
  validate_threshold_covariance_numeric_scalar(
    minimum_estimated_eigenvalue,
    "minimum_estimated_eigenvalue"
  )

  if (threshold_multiplier < 0) {
    stop("threshold_multiplier must be nonnegative.")
  }

  if (minimum_estimated_eigenvalue <= 0) {
    stop("minimum_estimated_eigenvalue must be greater than 0.")
  }

  estimated_group_1_mean <- colMeans(group_1_sample)
  estimated_group_2_mean <- colMeans(group_2_sample)

  centered_group_1_sample <- sweep(
    group_1_sample,
    MARGIN = 2,
    STATS = estimated_group_1_mean,
    FUN = "-"
  )
  centered_group_2_sample <- sweep(
    group_2_sample,
    MARGIN = 2,
    STATS = estimated_group_2_mean,
    FUN = "-"
  )
  pooled_centered_sample <- rbind(
    centered_group_1_sample,
    centered_group_2_sample
  )

  pooled_sample_size <- nrow(pooled_centered_sample)
  covariance_residual_degree_of_freedom <- pooled_sample_size - 2
  feature_number <- ncol(pooled_centered_sample)
  pooled_sample_covariance <- (
    crossprod(pooled_centered_sample) /
      covariance_residual_degree_of_freedom
  )

  estimated_variance_products <- outer(
    diag(pooled_sample_covariance),
    diag(pooled_sample_covariance),
    FUN = "*"
  )
  estimated_covariance_product_variances <- (
    estimated_variance_products + pooled_sample_covariance^2
  )
  entrywise_thresholds <- (
    threshold_multiplier *
      sqrt(
        estimated_covariance_product_variances *
          log(feature_number) /
          covariance_residual_degree_of_freedom
      )
  )

  retained_entry_indicator <- (
    abs(pooled_sample_covariance) >= entrywise_thresholds
  )
  diag(retained_entry_indicator) <- TRUE

  thresholded_covariance_matrix <- pooled_sample_covariance
  thresholded_covariance_matrix[!retained_entry_indicator] <- 0
  thresholded_covariance_matrix <- (
    thresholded_covariance_matrix + t(thresholded_covariance_matrix)
  ) / 2
  thresholded_covariance_matrix <- Matrix::Matrix(
    thresholded_covariance_matrix,
    sparse = TRUE
  )

  smallest_eigenvalue_before_shift <- (
    calculate_smallest_symmetric_eigenvalue(thresholded_covariance_matrix)
  )
  positive_semidefinite_shift <- max(
    0,
    minimum_estimated_eigenvalue - smallest_eigenvalue_before_shift
  )

  estimated_covariance_matrix <- thresholded_covariance_matrix
  if (positive_semidefinite_shift > 0) {
    estimated_covariance_matrix <- (
      estimated_covariance_matrix +
        Matrix::Diagonal(
          feature_number,
          x = positive_semidefinite_shift
        )
    )
  }

  off_diagonal_retained_entry_indicator <- retained_entry_indicator
  diag(off_diagonal_retained_entry_indicator) <- FALSE
  retained_off_diagonal_entry_number <- (
    sum(off_diagonal_retained_entry_indicator) / 2
  )

  diagnostics <- data.frame(
    feature_number = feature_number,
    pooled_sample_size = pooled_sample_size,
    covariance_residual_degree_of_freedom = (
      covariance_residual_degree_of_freedom
    ),
    threshold_multiplier = threshold_multiplier,
    retained_off_diagonal_entry_number = retained_off_diagonal_entry_number,
    retained_off_diagonal_proportion = (
      retained_off_diagonal_entry_number /
        (feature_number * (feature_number - 1) / 2)
    ),
    smallest_eigenvalue_before_shift = smallest_eigenvalue_before_shift,
    positive_semidefinite_shift = positive_semidefinite_shift,
    smallest_eigenvalue_after_shift = (
      smallest_eigenvalue_before_shift + positive_semidefinite_shift
    )
  )

  list(
    estimated_covariance_matrix = estimated_covariance_matrix,
    estimated_group_1_mean = estimated_group_1_mean,
    estimated_group_2_mean = estimated_group_2_mean,
    diagnostics = diagnostics
  )
}


orient_projection_direction <- function(projection_direction) {
  largest_loading_index <- which.max(abs(projection_direction))

  if (projection_direction[largest_loading_index] < 0) {
    projection_direction <- -projection_direction
  }

  as.numeric(projection_direction)
}


create_cross_fitting_fold_indices <- function(sample_size, fold_number) {
  if (fold_number < 2 || fold_number > sample_size) {
    stop("fold_number must be between 2 and the sample size.")
  }

  shuffled_observation_indices <- sample(seq_len(sample_size))
  split(
    shuffled_observation_indices,
    rep(seq_len(fold_number), length.out = sample_size)
  )
}


calculate_leading_eigenvector_pseudoinverse_action <- function(
    estimated_covariance_matrix,
    estimated_mean_difference,
    minimum_leading_eigenvalue_gap = 1e-8) {
  covariance_eigendecomposition <- eigen(
    as.matrix(estimated_covariance_matrix),
    symmetric = TRUE
  )
  leading_eigenvalue_gap <- (
    covariance_eigendecomposition$values[1] -
      covariance_eigendecomposition$values[2]
  )

  if (leading_eigenvalue_gap <= minimum_leading_eigenvalue_gap) {
    stop("The estimated leading eigenvalue is not sufficiently separated.")
  }

  estimated_leading_eigenvector <- orient_projection_direction(
    covariance_eigendecomposition$vectors[, 1]
  )
  remaining_eigenvectors <- covariance_eigendecomposition$vectors[
    ,
    -1,
    drop = FALSE
  ]
  inverse_eigenvalue_differences <- 1 / (
    covariance_eigendecomposition$values[1] -
      covariance_eigendecomposition$values[-1]
  )
  estimated_pseudoinverse_action <- as.numeric(
    remaining_eigenvectors %*%
      (
        inverse_eigenvalue_differences *
          crossprod(
            remaining_eigenvectors,
            estimated_mean_difference
          )
      )
  )

  list(
    estimated_leading_eigenvalue = (
      covariance_eigendecomposition$values[1]
    ),
    estimated_leading_eigenvalue_gap = leading_eigenvalue_gap,
    estimated_leading_eigenvector = estimated_leading_eigenvector,
    estimated_pseudoinverse_action = estimated_pseudoinverse_action
  )
}


estimate_debiased_pc_nuisance_parameters <- function(
    group_1_training_sample,
    group_2_training_sample,
    threshold_multiplier = 1.25,
    mean_difference_threshold_multiplier = 1) {
  covariance_estimation_result <- estimate_thresholded_common_covariance(
    group_1_sample = group_1_training_sample,
    group_2_sample = group_2_training_sample,
    threshold_multiplier = threshold_multiplier
  )
  unthresholded_estimated_mean_difference <- (
    covariance_estimation_result$estimated_group_1_mean -
      covariance_estimation_result$estimated_group_2_mean
  )
  estimated_mean_difference_standard_errors <- sqrt(
    apply(group_1_training_sample, 2, stats::var) /
      nrow(group_1_training_sample) +
      apply(group_2_training_sample, 2, stats::var) /
      nrow(group_2_training_sample)
  )
  mean_difference_thresholds <- (
    mean_difference_threshold_multiplier *
      sqrt(2 * log(ncol(group_1_training_sample))) *
      estimated_mean_difference_standard_errors
  )
  estimated_mean_difference <- unthresholded_estimated_mean_difference
  estimated_mean_difference[
    abs(estimated_mean_difference) < mean_difference_thresholds
  ] <- 0
  eigendecomposition_result <- (
    calculate_leading_eigenvector_pseudoinverse_action(
      estimated_covariance_matrix = (
        covariance_estimation_result$estimated_covariance_matrix
      ),
      estimated_mean_difference = estimated_mean_difference
    )
  )

  c(
    list(
      estimated_group_1_mean = (
        covariance_estimation_result$estimated_group_1_mean
      ),
      estimated_group_2_mean = (
        covariance_estimation_result$estimated_group_2_mean
      ),
      estimated_covariance_matrix = (
        covariance_estimation_result$estimated_covariance_matrix
      ),
      estimated_mean_difference = estimated_mean_difference,
      unthresholded_estimated_mean_difference = (
        unthresholded_estimated_mean_difference
      )
    ),
    eigendecomposition_result
  )
}


evaluate_debiased_pc_fold <- function(
    group_1_evaluation_sample,
    group_2_evaluation_sample,
    nuisance_parameter_estimates) {
  estimated_leading_eigenvector <- (
    nuisance_parameter_estimates$estimated_leading_eigenvector
  )
  estimated_pseudoinverse_action <- (
    nuisance_parameter_estimates$estimated_pseudoinverse_action
  )

  centered_group_1_evaluation_sample <- sweep(
    group_1_evaluation_sample,
    MARGIN = 2,
    STATS = nuisance_parameter_estimates$estimated_group_1_mean,
    FUN = "-"
  )
  centered_group_2_evaluation_sample <- sweep(
    group_2_evaluation_sample,
    MARGIN = 2,
    STATS = nuisance_parameter_estimates$estimated_group_2_mean,
    FUN = "-"
  )
  influence_function_centering_constant <- as.numeric(
    crossprod(
      estimated_pseudoinverse_action,
      nuisance_parameter_estimates$estimated_covariance_matrix %*%
        estimated_leading_eigenvector
    )
  )

  group_1_influence_function <- (
    as.numeric(
      centered_group_1_evaluation_sample %*%
        estimated_pseudoinverse_action
    ) *
      as.numeric(
        centered_group_1_evaluation_sample %*%
          estimated_leading_eigenvector
      ) -
      influence_function_centering_constant
  )
  group_2_influence_function <- (
    as.numeric(
      centered_group_2_evaluation_sample %*%
        estimated_pseudoinverse_action
    ) *
      as.numeric(
        centered_group_2_evaluation_sample %*%
          estimated_leading_eigenvector
      ) -
      influence_function_centering_constant
  )
  group_1_projected_observations <- as.numeric(
    group_1_evaluation_sample %*% estimated_leading_eigenvector
  )
  group_2_projected_observations <- as.numeric(
    group_2_evaluation_sample %*% estimated_leading_eigenvector
  )
  evaluation_sample_size_1 <- nrow(group_1_evaluation_sample)
  evaluation_sample_size_2 <- nrow(group_2_evaluation_sample)
  combined_evaluation_sample_size <- (
    evaluation_sample_size_1 + evaluation_sample_size_2
  )
  group_1_influence_weight <- (
    evaluation_sample_size_1 / combined_evaluation_sample_size
  )
  group_2_influence_weight <- (
    evaluation_sample_size_2 / combined_evaluation_sample_size
  )

  group_1_debiased_scores <- (
    group_1_projected_observations +
      group_1_influence_weight * group_1_influence_function
  )
  group_2_debiased_scores <- (
    group_2_projected_observations -
      group_2_influence_weight * group_2_influence_function
  )

  list(
    plug_in_estimate = (
      mean(group_1_projected_observations) -
        mean(group_2_projected_observations)
    ),
    plug_in_variance = (
      stats::var(group_1_projected_observations) /
        evaluation_sample_size_1 +
      stats::var(group_2_projected_observations) /
        evaluation_sample_size_2
    ),
    debiased_estimate = (
      mean(group_1_debiased_scores) - mean(group_2_debiased_scores)
    ),
    debiased_variance = (
      stats::var(group_1_debiased_scores) / evaluation_sample_size_1 +
      stats::var(group_2_debiased_scores) / evaluation_sample_size_2
    )
  )
}


combine_cross_fitted_fold_results <- function(
    fold_estimates,
    fold_variances) {
  fold_number <- length(fold_estimates)
  parameter_estimate <- mean(fold_estimates)
  standard_error <- sqrt(sum(fold_variances) / fold_number^2)
  test_statistic <- parameter_estimate / standard_error

  list(
    parameter_estimate = parameter_estimate,
    standard_error = standard_error,
    test_statistic = test_statistic,
    p_value = 2 * stats::pnorm(-abs(test_statistic))
  )
}


cross_fitted_pc_mean_comparison <- function(
    group_1_sample,
    group_2_sample,
    fold_number = 5,
    threshold_multiplier = 1.25,
    mean_difference_threshold_multiplier = 1,
    cross_fitting_seed = 1) {
  group_1_sample <- validate_covariance_estimation_sample(
    group_1_sample,
    "group_1_sample"
  )
  group_2_sample <- validate_covariance_estimation_sample(
    group_2_sample,
    "group_2_sample"
  )
  set.seed(cross_fitting_seed)
  group_1_fold_indices <- create_cross_fitting_fold_indices(
    nrow(group_1_sample),
    fold_number
  )
  group_2_fold_indices <- create_cross_fitting_fold_indices(
    nrow(group_2_sample),
    fold_number
  )
  fold_results <- vector("list", fold_number)

  for (fold_index in seq_len(fold_number)) {
    group_1_evaluation_indices <- group_1_fold_indices[[fold_index]]
    group_2_evaluation_indices <- group_2_fold_indices[[fold_index]]
    nuisance_parameter_estimates <- estimate_debiased_pc_nuisance_parameters(
      group_1_training_sample = group_1_sample[
        -group_1_evaluation_indices,
        ,
        drop = FALSE
      ],
      group_2_training_sample = group_2_sample[
        -group_2_evaluation_indices,
        ,
        drop = FALSE
      ],
      threshold_multiplier = threshold_multiplier,
      mean_difference_threshold_multiplier = (
        mean_difference_threshold_multiplier
      )
    )
    fold_results[[fold_index]] <- evaluate_debiased_pc_fold(
      group_1_evaluation_sample = group_1_sample[
        group_1_evaluation_indices,
        ,
        drop = FALSE
      ],
      group_2_evaluation_sample = group_2_sample[
        group_2_evaluation_indices,
        ,
        drop = FALSE
      ],
      nuisance_parameter_estimates = nuisance_parameter_estimates
    )
  }

  plug_in_result <- combine_cross_fitted_fold_results(
    fold_estimates = vapply(
      fold_results,
      function(fold_result) fold_result$plug_in_estimate,
      numeric(1)
    ),
    fold_variances = vapply(
      fold_results,
      function(fold_result) fold_result$plug_in_variance,
      numeric(1)
    )
  )
  debiased_result <- combine_cross_fitted_fold_results(
    fold_estimates = vapply(
      fold_results,
      function(fold_result) fold_result$debiased_estimate,
      numeric(1)
    ),
    fold_variances = vapply(
      fold_results,
      function(fold_result) fold_result$debiased_variance,
      numeric(1)
    )
  )

  list(
    plug_in = plug_in_result,
    debiased = debiased_result,
    fold_results = fold_results
  )
}


oracle_pc_mean_comparison <- function(
    group_1_sample,
    group_2_sample,
    population_leading_eigenvector) {
  population_leading_eigenvector <- as.numeric(
    population_leading_eigenvector
  )
  group_1_projected_observations <- as.numeric(
    group_1_sample %*% population_leading_eigenvector
  )
  group_2_projected_observations <- as.numeric(
    group_2_sample %*% population_leading_eigenvector
  )
  parameter_estimate <- (
    mean(group_1_projected_observations) -
      mean(group_2_projected_observations)
  )
  standard_error <- sqrt(
    stats::var(group_1_projected_observations) / nrow(group_1_sample) +
      stats::var(group_2_projected_observations) / nrow(group_2_sample)
  )
  test_statistic <- parameter_estimate / standard_error

  list(
    parameter_estimate = parameter_estimate,
    standard_error = standard_error,
    test_statistic = test_statistic,
    p_value = 2 * stats::pnorm(-abs(test_statistic))
  )
}


estimate_logistic_lasso_direction <- function(
    group_1_training_sample,
    group_2_training_sample,
    logistic_cross_validation_fold_number = 5) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("The glmnet package is required for the anchored test.")
  }

  pooled_training_sample <- rbind(
    group_1_training_sample,
    group_2_training_sample
  )
  group_labels <- c(
    rep(1, nrow(group_1_training_sample)),
    rep(0, nrow(group_2_training_sample))
  )
  logistic_fold_indices <- c(
    rep(
      seq_len(logistic_cross_validation_fold_number),
      length.out = nrow(group_1_training_sample)
    ),
    rep(
      seq_len(logistic_cross_validation_fold_number),
      length.out = nrow(group_2_training_sample)
    )
  )
  logistic_lasso_result <- glmnet::cv.glmnet(
    x = pooled_training_sample,
    y = group_labels,
    family = "binomial",
    foldid = logistic_fold_indices,
    type.measure = "deviance",
    standardize = TRUE,
    nlambda = 50
  )

  as.numeric(
    stats::coef(
      logistic_lasso_result,
      s = "lambda.1se"
    )[-1, 1]
  )
}


anchored_lasso_mean_comparison <- function(
    group_1_sample,
    group_2_sample,
    fold_number = 5,
    threshold_multiplier = 1.25,
    cross_fitting_seed = 1) {
  set.seed(cross_fitting_seed)
  group_1_fold_indices <- create_cross_fitting_fold_indices(
    nrow(group_1_sample),
    fold_number
  )
  group_2_fold_indices <- create_cross_fitting_fold_indices(
    nrow(group_2_sample),
    fold_number
  )
  fold_estimates <- numeric(fold_number)
  fold_variances <- numeric(fold_number)

  for (fold_index in seq_len(fold_number)) {
    group_1_evaluation_indices <- group_1_fold_indices[[fold_index]]
    group_2_evaluation_indices <- group_2_fold_indices[[fold_index]]
    group_1_training_sample <- group_1_sample[
      -group_1_evaluation_indices,
      ,
      drop = FALSE
    ]
    group_2_training_sample <- group_2_sample[
      -group_2_evaluation_indices,
      ,
      drop = FALSE
    ]
    nuisance_parameter_estimates <- estimate_debiased_pc_nuisance_parameters(
      group_1_training_sample = group_1_training_sample,
      group_2_training_sample = group_2_training_sample,
      threshold_multiplier = threshold_multiplier
    )
    estimated_logistic_lasso_direction <- estimate_logistic_lasso_direction(
      group_1_training_sample,
      group_2_training_sample
    )
    effective_training_sample_size <- (
      nrow(group_1_training_sample) + nrow(group_2_training_sample)
    )
    logistic_direction_threshold <- (
      effective_training_sample_size^(-1 / 3)
    )
    estimated_projection_direction <- (
      nuisance_parameter_estimates$estimated_leading_eigenvector
    )

    if (
      sqrt(sum(estimated_logistic_lasso_direction^2)) >=
        logistic_direction_threshold
    ) {
      estimated_projection_direction <- (
        estimated_projection_direction +
          sqrt(effective_training_sample_size) *
            estimated_logistic_lasso_direction
      )
    }
    estimated_projection_direction <- (
      estimated_projection_direction /
        sqrt(sum(estimated_projection_direction^2))
    )

    group_1_projected_observations <- as.numeric(
      group_1_sample[group_1_evaluation_indices, , drop = FALSE] %*%
        estimated_projection_direction
    )
    group_2_projected_observations <- as.numeric(
      group_2_sample[group_2_evaluation_indices, , drop = FALSE] %*%
        estimated_projection_direction
    )
    fold_estimates[fold_index] <- (
      mean(group_1_projected_observations) -
        mean(group_2_projected_observations)
    )
    fold_variances[fold_index] <- (
      stats::var(group_1_projected_observations) /
        length(group_1_projected_observations) +
      stats::var(group_2_projected_observations) /
        length(group_2_projected_observations)
    )
  }

  combine_cross_fitted_fold_results(
    fold_estimates = fold_estimates,
    fold_variances = fold_variances
  )
}


generate_gaussian_two_sample <- function(
    sample_size_per_group,
    population_covariance_matrix,
    population_mean_difference) {
  feature_number <- nrow(population_covariance_matrix)
  population_covariance_cholesky_factor <- chol(
    population_covariance_matrix
  )
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
