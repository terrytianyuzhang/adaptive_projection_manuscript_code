# Run as: Rscript 012_run_mean_comparison_task.R simulation_task_id

# File input and output -------------------------------------------------------

mean_comparison_source_file <- file.path(
  "..",
  "src",
  "002_threshold_sparse_covariance.R"
)
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


# Source computational functions --------------------------------------------

source(mean_comparison_source_file)
source(simulation_task_source_file)


# Identify one task ----------------------------------------------------------

command_arguments <- commandArgs(trailingOnly = TRUE)
if (length(command_arguments) != 1L) {
  stop("Provide exactly one unique_id.")
}
unique_id <- command_arguments[1]

simulation_task_manifest <- read.csv(
  simulation_task_manifest_file,
  stringsAsFactors = FALSE
)
task_match_indicator <- (
  simulation_task_manifest$unique_id == unique_id
)
if (sum(task_match_indicator) != 1L) {
  stop("unique_id was not found uniquely: ", unique_id)
}
task_specification <- simulation_task_manifest[
  task_match_indicator,
  ,
  drop = FALSE
]
simulation_setting_table <- read.csv(
  simulation_setting_table_file,
  stringsAsFactors = FALSE
)
setting_match_indicator <- (
  simulation_setting_table$setting == task_specification$setting
)
if (sum(setting_match_indicator) != 1L) {
  stop("setting was not found uniquely: ", task_specification$setting)
}
setting_specification <- simulation_setting_table[
  setting_match_indicator,
  ,
  drop = FALSE
]

feature_name <- paste0(
  "feature_",
  setting_specification$feature_number
)
shared_simulation_input <- readRDS(
  shared_simulation_input_files[[feature_name]]
)
current_result_file <- create_mean_comparison_task_result_file(
  simulation_task_result_directory,
  task_specification,
  shared_simulation_input$simulation_settings$simulation_method_ids
)


# Resume a completed task ----------------------------------------------------

if (file.exists(current_result_file)) {
  existing_task_output <- tryCatch(
    readRDS(current_result_file),
    error = function(read_error) NULL
  )
  existing_result_is_valid <- validate_mean_comparison_task_output(
    task_output = existing_task_output,
    expected_task_specification = task_specification,
    expected_setting_specification = setting_specification,
    expected_simulation_settings = (
      shared_simulation_input$simulation_settings
    )
  )

  if (existing_result_is_valid) {
    message("Task ", unique_id, " is already complete.")
    quit(save = "no")
  }

  stop(
    "The existing result for ", unique_id,
    " is incomplete or incompatible. Remove it before rerunning this task."
  )
}


# Run one batch --------------------------------------------------------------

message(
  "Running ", unique_id,
  ": ", task_specification$method,
  ", p = ", setting_specification$feature_number,
  ", n = ", setting_specification$sample_size_per_group,
  ", batch ", task_specification$batch_number, "."
)

repeat_indices <- seq.int(
  (task_specification$batch_number - 1L) *
    shared_simulation_input$simulation_settings$
      repeat_number_per_batch + 1L,
  task_specification$batch_number *
    shared_simulation_input$simulation_settings$
      repeat_number_per_batch
)
batch_simulation_results <- vector(
  "list",
  length(repeat_indices)
)

for (repeat_number in seq_along(repeat_indices)) {
  repeat_index <- repeat_indices[repeat_number]
  simulation_seed <- (
    shared_simulation_input$simulation_settings$
      mean_comparison_simulation_seed +
      setting_specification$feature_number * 1000000L +
      setting_specification$sample_size_per_group * 1000L +
      repeat_index
  )

  batch_simulation_results[[repeat_number]] <- (
    run_one_mean_comparison_simulation(
      sample_size_per_group = setting_specification$sample_size_per_group,
      repeat_index = repeat_index,
      simulation_seed = simulation_seed,
      shared_simulation_input = shared_simulation_input,
      method = task_specification$method
    )
  )
}

task_output <- list(
  task_specification = task_specification,
  simulation_results = do.call(rbind, batch_simulation_results),
  simulation_settings = shared_simulation_input$simulation_settings
)

if (!validate_mean_comparison_task_output(
  task_output = task_output,
  expected_task_specification = task_specification,
  expected_setting_specification = setting_specification,
  expected_simulation_settings = shared_simulation_input$simulation_settings
)) {
  stop("Internal validation failed for task ", unique_id, ".")
}

temporary_result_file <- paste0(
  current_result_file,
  ".temporary_",
  Sys.getpid()
)
on.exit(unlink(temporary_result_file), add = TRUE)
saveRDS(task_output, file = temporary_result_file)

if (!file.rename(temporary_result_file, current_result_file)) {
  stop("Could not move the completed task result into place.")
}

message("Completed task ", unique_id, ".")
