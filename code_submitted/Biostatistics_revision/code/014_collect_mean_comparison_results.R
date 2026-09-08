# File input and output -------------------------------------------------------

simulation_task_source_file <- file.path(
  "..",
  "src",
  "012_mean_comparison_simulation_task.R"
)
simulation_task_manifest_file <- file.path(
  "..",
  "data",
  "intermediate",
  "011_mean_comparison_task_manifest.csv"
)
simulation_setting_table_file <- file.path(
  "..",
  "data",
  "intermediate",
  "011_mean_comparison_setting_table.csv"
)
shared_simulation_input_files <- c(
  feature_100 = file.path(
    "..", "data", "intermediate",
    "011_mean_comparison_shared_input_p100.rds"
  ),
  feature_1000 = file.path(
    "..", "data", "intermediate",
    "011_mean_comparison_shared_input_p1000.rds"
  )
)
simulation_task_result_directory <- file.path(
  "..",
  "data",
  "intermediate",
  "012_mean_comparison_task_result"
)
collected_simulation_output_file <- file.path(
  "..",
  "data",
  "intermediate",
  "014_mean_comparison_simulation.rds"
)
rejection_rate_summary_output_file <- file.path(
  "..",
  "data",
  "final",
  "014_mean_comparison_rejection_rate_summary.rds"
)


# Source validation functions and read task definitions ----------------------

source(simulation_task_source_file)

simulation_task_manifest <- read.csv(
  simulation_task_manifest_file,
  stringsAsFactors = FALSE
)
simulation_setting_table <- read.csv(
  simulation_setting_table_file,
  stringsAsFactors = FALSE
)
shared_simulation_inputs <- lapply(
  shared_simulation_input_files,
  readRDS
)
reference_simulation_settings <- (
  shared_simulation_inputs[[1]]$simulation_settings
)

simulation_settings_agree <- vapply(
  shared_simulation_inputs,
  function(shared_simulation_input) {
    identical(
      shared_simulation_input$simulation_settings,
      reference_simulation_settings
    )
  },
  logical(1)
)
if (!all(simulation_settings_agree)) {
  stop("The shared simulation inputs contain different settings.")
}
if (!identical(
  names(simulation_setting_table),
  c(
    "setting", "feature_number", "sample_size_per_group",
    "population_mean_structure"
  )
) || anyDuplicated(simulation_setting_table$setting) > 0L) {
  stop("The simulation setting table is invalid.")
}
if (!identical(
  simulation_setting_table,
  reference_simulation_settings$simulation_setting_table
)) {
  stop("The simulation setting table does not match the shared inputs.")
}

expected_task_number <- (
  nrow(reference_simulation_settings$simulation_setting_table) *
    reference_simulation_settings$simulation_repeat_number /
    reference_simulation_settings$repeat_number_per_batch *
    length(reference_simulation_settings$simulation_method_ids)
)
if (
  nrow(simulation_task_manifest) != expected_task_number ||
    !identical(
      names(simulation_task_manifest),
      c("unique_id", "method", "setting", "batch_number")
    ) ||
    anyDuplicated(simulation_task_manifest$unique_id) > 0L
) {
  stop("The simulation task manifest is incomplete or has duplicate IDs.")
}

simulation_task_result_files <- vapply(
  seq_len(nrow(simulation_task_manifest)),
  function(task_number) {
    create_mean_comparison_task_result_file(
      simulation_task_result_directory = simulation_task_result_directory,
      task_specification = simulation_task_manifest[
        task_number,
        ,
        drop = FALSE
      ],
      simulation_method_ids = (
        reference_simulation_settings$simulation_method_ids
      )
    )
  },
  character(1)
)
missing_task_ids <- simulation_task_manifest$unique_id[
  !file.exists(simulation_task_result_files)
]

if (length(missing_task_ids) > 0L) {
  stop(
    "The following simulation tasks have no result: ",
    paste(missing_task_ids, collapse = ", ")
  )
}


# Validate and collect task results -----------------------------------------

task_outputs <- vector("list", nrow(simulation_task_manifest))
invalid_task_ids <- character()

for (task_number in seq_len(nrow(simulation_task_manifest))) {
  task_specification <- simulation_task_manifest[
    task_number,
    ,
    drop = FALSE
  ]
  task_output <- tryCatch(
    readRDS(simulation_task_result_files[task_number]),
    error = function(read_error) NULL
  )
  setting_specification <- simulation_setting_table[
    simulation_setting_table$setting == task_specification$setting,
    ,
    drop = FALSE
  ]
  if (nrow(setting_specification) != 1L) {
    stop("setting was not found uniquely: ", task_specification$setting)
  }
  feature_name <- paste0(
    "feature_",
    setting_specification$feature_number
  )
  task_output_is_valid <- validate_mean_comparison_task_output(
    task_output = task_output,
    expected_task_specification = task_specification,
    expected_setting_specification = setting_specification,
    expected_simulation_settings = (
      shared_simulation_inputs[[feature_name]]$simulation_settings
    )
  )

  if (!task_output_is_valid) {
    invalid_task_ids <- c(
      invalid_task_ids,
      task_specification$unique_id
    )
  } else {
    task_outputs[[task_number]] <- task_output
  }
}

if (length(invalid_task_ids) > 0L) {
  stop(
    "The following simulation task results are invalid: ",
    paste(invalid_task_ids, collapse = ", ")
  )
}

simulation_results <- do.call(
  rbind,
  lapply(
    task_outputs,
    function(task_output) task_output$simulation_results
  )
)
simulation_results <- simulation_results[
  order(
    simulation_results$feature_number,
    simulation_results$sample_size_per_group,
    simulation_results$repeat_index,
    simulation_results$simulation_target,
    simulation_results$method
  ),
  ,
  drop = FALSE
]

result_combination_counts <- aggregate(
  list(result_number = simulation_results$rejected),
  by = list(
    feature_number = simulation_results$feature_number,
    sample_size_per_group = simulation_results$sample_size_per_group,
    repeat_index = simulation_results$repeat_index,
    simulation_target = simulation_results$simulation_target,
    method = simulation_results$method
  ),
  FUN = length
)

if (
  nrow(result_combination_counts) != nrow(simulation_results) ||
    any(result_combination_counts$result_number != 1L)
) {
  stop("Collected results contain duplicated simulation result rows.")
}

simulation_coverage <- aggregate(
  repeat_index ~ feature_number + sample_size_per_group +
    simulation_target + method,
  data = simulation_results,
  FUN = function(repeat_index) length(unique(repeat_index))
)
expected_coverage_row_number <- (
  nrow(reference_simulation_settings$simulation_setting_table) *
    2L * length(reference_simulation_settings$simulation_method_ids)
)
if (
  nrow(simulation_coverage) != expected_coverage_row_number ||
    any(
      simulation_coverage$repeat_index !=
        reference_simulation_settings$simulation_repeat_number
    )
) {
  stop(
    "The collected results do not contain every requested repetition ",
    "for every simulation setting and method."
  )
}


# Summarize rejection rates --------------------------------------------------

rejection_rate_summary <- aggregate(
  rejected ~ feature_number + sample_size_per_group +
    simulation_target + method,
  data = simulation_results,
  FUN = function(rejected) {
    c(
      rejection_rate = mean(rejected),
      monte_carlo_standard_error = sqrt(
        mean(rejected) * (1 - mean(rejected)) / length(rejected)
      ),
      repeat_number = length(rejected)
    )
  }
)
rejection_rate_summary <- data.frame(
  feature_number = rejection_rate_summary$feature_number,
  sample_size_per_group = rejection_rate_summary$sample_size_per_group,
  simulation_target = rejection_rate_summary$simulation_target,
  method = rejection_rate_summary$method,
  rejection_rate = rejection_rate_summary$rejected[, "rejection_rate"],
  monte_carlo_standard_error = (
    rejection_rate_summary$rejected[, "monte_carlo_standard_error"]
  ),
  repeat_number = rejection_rate_summary$rejected[, "repeat_number"],
  significance_level = reference_simulation_settings$significance_level,
  stringsAsFactors = FALSE
)
rejection_rate_summary$lower_plot_limit <- pmax(
  0,
  rejection_rate_summary$rejection_rate -
    1.96 * rejection_rate_summary$monte_carlo_standard_error
)
rejection_rate_summary$upper_plot_limit <- pmin(
  1,
  rejection_rate_summary$rejection_rate +
    1.96 * rejection_rate_summary$monte_carlo_standard_error
)

collected_simulation_output <- list(
  simulation_results = simulation_results,
  simulation_settings = reference_simulation_settings,
  task_manifest = simulation_task_manifest
)

dir.create(
  dirname(collected_simulation_output_file),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  dirname(rejection_rate_summary_output_file),
  recursive = TRUE,
  showWarnings = FALSE
)
saveRDS(
  collected_simulation_output,
  file = collected_simulation_output_file
)
saveRDS(
  rejection_rate_summary,
  file = rejection_rate_summary_output_file
)

message(
  "Collected ", nrow(simulation_results),
  " result rows from ", nrow(simulation_task_manifest),
  " simulation tasks."
)
print(rejection_rate_summary)
