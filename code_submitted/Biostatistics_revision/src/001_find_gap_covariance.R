validate_finite_numeric_scalar <- function(value, value_name) {
  if (
    !is.numeric(value) ||
      length(value) != 1 ||
      is.na(value) ||
      !is.finite(value)
  ) {
    stop(paste(value_name, "must be one finite numeric value."))
  }
}


validate_positive_integer <- function(value, value_name, minimum_value = 1) {
  validate_finite_numeric_scalar(value, value_name)

  if (value < minimum_value || value != floor(value)) {
    stop(
      paste(
        value_name,
        "must be an integer greater than or equal to",
        paste0(minimum_value, ".")
      )
    )
  }
}


validate_sparse_covariance_settings <- function(
    feature_number,
    average_nonzero_off_diagonal_number_per_row,
    covariance_perturbation_size,
    minimum_leading_eigenvalue_gap,
    maximum_search_attempt_number,
    covariance_seed,
    leading_eigenpair_number,
    minimum_edge_weight,
    maximum_edge_weight) {
  validate_positive_integer(feature_number, "feature_number", minimum_value = 3)
  validate_finite_numeric_scalar(
    average_nonzero_off_diagonal_number_per_row,
    "average_nonzero_off_diagonal_number_per_row"
  )
  validate_finite_numeric_scalar(
    covariance_perturbation_size,
    "covariance_perturbation_size"
  )
  validate_finite_numeric_scalar(
    minimum_leading_eigenvalue_gap,
    "minimum_leading_eigenvalue_gap"
  )
  validate_positive_integer(
    maximum_search_attempt_number,
    "maximum_search_attempt_number"
  )
  validate_positive_integer(covariance_seed, "covariance_seed", minimum_value = 0)
  validate_positive_integer(
    leading_eigenpair_number,
    "leading_eigenpair_number",
    minimum_value = 2
  )
  validate_finite_numeric_scalar(minimum_edge_weight, "minimum_edge_weight")
  validate_finite_numeric_scalar(maximum_edge_weight, "maximum_edge_weight")

  if (
    average_nonzero_off_diagonal_number_per_row <= 0 ||
      average_nonzero_off_diagonal_number_per_row > feature_number - 1
  ) {
    stop(
      paste(
        "average_nonzero_off_diagonal_number_per_row must be greater than 0",
        "and no greater than feature_number - 1."
      )
    )
  }

  if (covariance_perturbation_size <= 0 || covariance_perturbation_size >= 1) {
    stop("covariance_perturbation_size must be strictly between 0 and 1.")
  }

  if (minimum_leading_eigenvalue_gap <= 0) {
    stop("minimum_leading_eigenvalue_gap must be greater than 0.")
  }

  if (minimum_leading_eigenvalue_gap >= 2 * covariance_perturbation_size) {
    stop(
      paste(
        "minimum_leading_eigenvalue_gap must be smaller than",
        "2 * covariance_perturbation_size."
      )
    )
  }

  if (leading_eigenpair_number >= feature_number) {
    stop("leading_eigenpair_number must be smaller than feature_number.")
  }

  if (covariance_seed > .Machine$integer.max) {
    stop("covariance_seed must not exceed .Machine$integer.max.")
  }

  if (minimum_edge_weight <= 0) {
    stop("minimum_edge_weight must be greater than 0.")
  }

  if (maximum_edge_weight < minimum_edge_weight) {
    stop("maximum_edge_weight must be at least minimum_edge_weight.")
  }
}


generate_random_sparse_symmetric_matrix <- function(
    feature_number,
    edge_number,
    minimum_edge_weight,
    maximum_edge_weight) {
  selected_edge_keys <- numeric(0)

  while (length(selected_edge_keys) < edge_number) {
    remaining_edge_number <- edge_number - length(selected_edge_keys)
    candidate_edge_number <- max(100, 2 * remaining_edge_number)

    candidate_first_features <- sample.int(
      feature_number,
      size = candidate_edge_number,
      replace = TRUE
    )
    candidate_second_features <- sample.int(
      feature_number,
      size = candidate_edge_number,
      replace = TRUE
    )

    distinct_feature_indicator <- (
      candidate_first_features != candidate_second_features
    )
    candidate_first_features <- candidate_first_features[
      distinct_feature_indicator
    ]
    candidate_second_features <- candidate_second_features[
      distinct_feature_indicator
    ]

    candidate_smaller_features <- pmin(
      candidate_first_features,
      candidate_second_features
    )
    candidate_larger_features <- pmax(
      candidate_first_features,
      candidate_second_features
    )

    candidate_edge_keys <- (
      candidate_smaller_features +
        feature_number * (candidate_larger_features - 1)
    )
    selected_edge_keys <- unique(c(selected_edge_keys, candidate_edge_keys))

    if (length(selected_edge_keys) > edge_number) {
      selected_edge_keys <- selected_edge_keys[seq_len(edge_number)]
    }
  }

  smaller_features <- ((selected_edge_keys - 1) %% feature_number) + 1
  larger_features <- (
    (selected_edge_keys - smaller_features) / feature_number
  ) + 1

  edge_weights <- stats::runif(
    edge_number,
    min = minimum_edge_weight,
    max = maximum_edge_weight
  )

  Matrix::sparseMatrix(
    i = c(smaller_features, larger_features),
    j = c(larger_features, smaller_features),
    x = c(edge_weights, edge_weights),
    dims = c(feature_number, feature_number),
    symmetric = FALSE
  )
}


calculate_extreme_eigenvalues <- function(sparse_symmetric_matrix) {
  two_largest_eigenpairs <- tryCatch(
    RSpectra::eigs_sym(
      sparse_symmetric_matrix,
      k = 2,
      which = "LA"
    ),
    error = function(error_condition) NULL
  )

  smallest_eigenpair <- tryCatch(
    RSpectra::eigs_sym(
      sparse_symmetric_matrix,
      k = 1,
      which = "SA"
    ),
    error = function(error_condition) NULL
  )

  if (is.null(two_largest_eigenpairs) || is.null(smallest_eigenpair)) {
    return(NULL)
  }

  decreasing_eigenvalue_order <- order(
    two_largest_eigenpairs$values,
    decreasing = TRUE
  )

  list(
    two_largest_eigenvalues = two_largest_eigenpairs$values[
      decreasing_eigenvalue_order
    ],
    smallest_eigenvalue = smallest_eigenpair$values[1]
  )
}


calculate_leading_eigenpairs <- function(
    sparse_symmetric_matrix,
    leading_eigenpair_number) {
  leading_eigenpairs <- RSpectra::eigs_sym(
    sparse_symmetric_matrix,
    k = leading_eigenpair_number,
    which = "LA"
  )

  decreasing_eigenvalue_order <- order(
    leading_eigenpairs$values,
    decreasing = TRUE
  )

  leading_eigenvalues <- leading_eigenpairs$values[
    decreasing_eigenvalue_order
  ]
  leading_eigenvectors <- leading_eigenpairs$vectors[
    ,
    decreasing_eigenvalue_order,
    drop = FALSE
  ]

  for (eigenvector_number in seq_len(ncol(leading_eigenvectors))) {
    largest_loading_position <- which.max(
      abs(leading_eigenvectors[, eigenvector_number])
    )

    if (leading_eigenvectors[largest_loading_position, eigenvector_number] < 0) {
      leading_eigenvectors[, eigenvector_number] <- (
        -leading_eigenvectors[, eigenvector_number]
      )
    }
  }

  list(
    leading_eigenvalues = leading_eigenvalues,
    leading_eigenvectors = leading_eigenvectors
  )
}


preserve_random_seed <- function() {
  random_seed_existed <- exists(
    ".Random.seed",
    envir = .GlobalEnv,
    inherits = FALSE
  )

  if (random_seed_existed) {
    previous_random_seed <- get(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )
  } else {
    previous_random_seed <- NULL
  }

  list(
    random_seed_existed = random_seed_existed,
    previous_random_seed = previous_random_seed
  )
}


restore_random_seed <- function(random_seed_state) {
  if (random_seed_state$random_seed_existed) {
    assign(
      ".Random.seed",
      random_seed_state$previous_random_seed,
      envir = .GlobalEnv
    )
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    remove(".Random.seed", envir = .GlobalEnv)
  }
}


find_sparse_covariance_matrix <- function(
    feature_number,
    average_nonzero_off_diagonal_number_per_row,
    covariance_perturbation_size,
    minimum_leading_eigenvalue_gap,
    maximum_search_attempt_number,
    covariance_seed,
    leading_eigenpair_number = 5,
    minimum_edge_weight = 0.5,
    maximum_edge_weight = 1.5) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The Matrix package is required.")
  }

  if (!requireNamespace("RSpectra", quietly = TRUE)) {
    stop("The RSpectra package is required.")
  }

  validate_sparse_covariance_settings(
    feature_number = feature_number,
    average_nonzero_off_diagonal_number_per_row = (
      average_nonzero_off_diagonal_number_per_row
    ),
    covariance_perturbation_size = covariance_perturbation_size,
    minimum_leading_eigenvalue_gap = minimum_leading_eigenvalue_gap,
    maximum_search_attempt_number = maximum_search_attempt_number,
    covariance_seed = covariance_seed,
    leading_eigenpair_number = leading_eigenpair_number,
    minimum_edge_weight = minimum_edge_weight,
    maximum_edge_weight = maximum_edge_weight
  )

  edge_number <- round(
    feature_number * average_nonzero_off_diagonal_number_per_row / 2
  )
  maximum_edge_number <- feature_number * (feature_number - 1) / 2
  edge_number <- min(edge_number, maximum_edge_number)

  random_seed_state <- preserve_random_seed()
  on.exit(restore_random_seed(random_seed_state), add = TRUE)
  set.seed(covariance_seed)

  failed_eigenvalue_attempt_number <- 0

  for (search_attempt_number in seq_len(maximum_search_attempt_number)) {
    sparse_symmetric_matrix <- generate_random_sparse_symmetric_matrix(
      feature_number = feature_number,
      edge_number = edge_number,
      minimum_edge_weight = minimum_edge_weight,
      maximum_edge_weight = maximum_edge_weight
    )

    extreme_eigenvalues <- calculate_extreme_eigenvalues(
      sparse_symmetric_matrix
    )

    if (is.null(extreme_eigenvalues)) {
      failed_eigenvalue_attempt_number <- failed_eigenvalue_attempt_number + 1
      next
    }

    spectral_scaling_value <- max(
      abs(extreme_eigenvalues$two_largest_eigenvalues[1]),
      abs(extreme_eigenvalues$smallest_eigenvalue)
    )
    leading_eigenvalue_gap <- (
      covariance_perturbation_size *
        (
          extreme_eigenvalues$two_largest_eigenvalues[1] -
            extreme_eigenvalues$two_largest_eigenvalues[2]
        ) /
        spectral_scaling_value
    )

    if (leading_eigenvalue_gap < minimum_leading_eigenvalue_gap) {
      next
    }

    leading_eigenpairs <- calculate_leading_eigenpairs(
      sparse_symmetric_matrix = sparse_symmetric_matrix,
      leading_eigenpair_number = leading_eigenpair_number
    )

    covariance_matrix <- (
      Matrix::Diagonal(feature_number) +
        covariance_perturbation_size *
        sparse_symmetric_matrix /
        spectral_scaling_value
    )
    covariance_leading_eigenvalues <- (
      1 +
        covariance_perturbation_size *
        leading_eigenpairs$leading_eigenvalues /
        spectral_scaling_value
    )
    covariance_smallest_eigenvalue <- (
      1 +
        covariance_perturbation_size *
        extreme_eigenvalues$smallest_eigenvalue /
        spectral_scaling_value
    )

    nonzero_off_diagonal_number_per_row <- Matrix::rowSums(
      sparse_symmetric_matrix != 0
    )
    realized_average_nonzero_off_diagonal_number_per_row <- mean(
      nonzero_off_diagonal_number_per_row
    )

    if (!isTRUE(Matrix::isSymmetric(covariance_matrix))) {
      stop("The accepted covariance matrix is not symmetric.")
    }

    if (covariance_smallest_eigenvalue <= 0) {
      stop("The accepted covariance matrix is not positive definite.")
    }

    if (
      covariance_leading_eigenvalues[1] -
        covariance_leading_eigenvalues[2] <
        minimum_leading_eigenvalue_gap
    ) {
      stop("The accepted covariance matrix does not have the requested gap.")
    }

    generation_settings <- list(
      feature_number = feature_number,
      average_nonzero_off_diagonal_number_per_row = (
        average_nonzero_off_diagonal_number_per_row
      ),
      covariance_perturbation_size = covariance_perturbation_size,
      minimum_leading_eigenvalue_gap = minimum_leading_eigenvalue_gap,
      maximum_search_attempt_number = maximum_search_attempt_number,
      covariance_seed = covariance_seed,
      leading_eigenpair_number = leading_eigenpair_number,
      minimum_edge_weight = minimum_edge_weight,
      maximum_edge_weight = maximum_edge_weight
    )

    diagnostics <- data.frame(
      feature_number = feature_number,
      edge_number = edge_number,
      average_nonzero_off_diagonal_number_per_row = (
        realized_average_nonzero_off_diagonal_number_per_row
      ),
      minimum_nonzero_off_diagonal_number_per_row = min(
        nonzero_off_diagonal_number_per_row
      ),
      maximum_nonzero_off_diagonal_number_per_row = max(
        nonzero_off_diagonal_number_per_row
      ),
      leading_eigenvalue = covariance_leading_eigenvalues[1],
      second_eigenvalue = covariance_leading_eigenvalues[2],
      leading_eigenvalue_gap = (
        covariance_leading_eigenvalues[1] -
          covariance_leading_eigenvalues[2]
      ),
      smallest_eigenvalue = covariance_smallest_eigenvalue,
      condition_number = (
        covariance_leading_eigenvalues[1] /
          covariance_smallest_eigenvalue
      ),
      search_attempt_number = search_attempt_number,
      failed_eigenvalue_attempt_number = failed_eigenvalue_attempt_number,
      covariance_seed = covariance_seed
    )

    return(list(
      covariance_matrix = covariance_matrix,
      leading_eigenvalues = covariance_leading_eigenvalues,
      leading_eigenvectors = leading_eigenpairs$leading_eigenvectors,
      nonzero_off_diagonal_number_per_row = (
        nonzero_off_diagonal_number_per_row
      ),
      generation_settings = generation_settings,
      diagnostics = diagnostics
    ))
  }

  stop(
    paste(
      "No sparse covariance matrix attained the requested leading eigenvalue",
      "gap after", maximum_search_attempt_number, "attempts.",
      "Increase maximum_search_attempt_number or decrease",
      "minimum_leading_eigenvalue_gap."
    )
  )
}
