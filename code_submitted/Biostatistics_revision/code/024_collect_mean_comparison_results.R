# File input and output ------------------------------------------------------

simulation_task_table_input_file <- file.path(
  "..",
  "data",
  "intermediate",
  "011_mean_comparison_task_table.csv"
)
simulation_config_source_file <- file.path(
  "..",
  "src",
  "011_mean_comparison_simulation_config.R"
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


# Testing constant ----------------------------------------------------------

significance_level <- 0.05


# Read task results ---------------------------------------------------------

simulation_task_table <- read.csv(
  simulation_task_table_input_file,
  stringsAsFactors = FALSE
)
simulation_task_result_files <- file.path(
  simulation_task_result_directory,
  paste0(
    "012_mean_comparison_task_",
    simulation_task_table$id,
    ".rds"
  )
)

missing_task_ids <- simulation_task_table$id[
  !file.exists(simulation_task_result_files)
]
if (length(missing_task_ids) > 0L) {
  stop(
    "Missing results for task IDs: ",
    paste(missing_task_ids, collapse = ", ")
  )
}

task_result_tables <- vector(
  "list",
  nrow(simulation_task_table)
)

for (task_number in 1:nrow(simulation_task_table)) {
  repetition_results <- readRDS(
    simulation_task_result_files[task_number]
  )
  repetition_result_table <- do.call(
    rbind,
    lapply(
      repetition_results,
      as.data.frame,
      stringsAsFactors = FALSE
    )
  )
  repetition_result_table$task_id <- (
    simulation_task_table$id[task_number]
  )
  repetition_result_table$batch_number <- (
    simulation_task_table$batch_number[task_number]
  )
  task_result_tables[[task_number]] <- repetition_result_table
}

simulation_results <- do.call(
  rbind,
  task_result_tables
)
rownames(simulation_results) <- NULL
names(simulation_results)[
  names(simulation_results) == "method"
] <- "method_id"


# Add setting information ---------------------------------------------------

setting_information_tables <- vector(
  "list",
  length(unique(simulation_task_table$setting_id))
)
unique_setting_ids <- sort(unique(simulation_task_table$setting_id))

for (setting_number in 1:length(unique_setting_ids)) {
  setting_id <- unique_setting_ids[setting_number]
  source(simulation_config_source_file)

  if (setting_id %in% c(
    1L, 2L, 3L, 4L, 5L,
    6L, 7L, 8L, 9L, 10L
  )) {
    simulation_target <- "Type I error"
  } else {
    simulation_target <- "Power"
  }

  setting_information_tables[[setting_number]] <- data.frame(
    setting_id = setting_id,
    feature_number = simulation_config$dimension,
    sample_size_per_group = simulation_config$control_sample_size,
    simulation_target = simulation_target,
    stringsAsFactors = FALSE
  )
}

setting_information_table <- do.call(
  rbind,
  setting_information_tables
)
simulation_results <- merge(
  simulation_results,
  setting_information_table,
  by = "setting_id"
)

simulation_results$method <- NA_character_
simulation_results$method[
  simulation_results$method_id == "debiased_pc"
] <- "Debiased PC"
simulation_results$method[
  simulation_results$method_id == "oracle_pc"
] <- "Oracle PC"
simulation_results$method[
  simulation_results$method_id == "plug_in_pc"
] <- "Plug-in PC"
simulation_results$rejected <- (
  simulation_results$p_value <= significance_level
)


# Summarize rejection rates -------------------------------------------------

summary_group <- interaction(
  simulation_results$feature_number,
  simulation_results$sample_size_per_group,
  simulation_results$simulation_target,
  simulation_results$method,
  drop = TRUE
)
simulation_result_groups <- split(
  simulation_results,
  summary_group
)

rejection_rate_summary_tables <- lapply(
  simulation_result_groups,
  function(simulation_result_group) {
    rejection_rate <- mean(simulation_result_group$rejected)
    monte_carlo_standard_error <- sqrt(
      rejection_rate * (1 - rejection_rate) /
        nrow(simulation_result_group)
    )

    data.frame(
      feature_number = simulation_result_group$feature_number[1L],
      sample_size_per_group = (
        simulation_result_group$sample_size_per_group[1L]
      ),
      simulation_target = simulation_result_group$simulation_target[1L],
      method = simulation_result_group$method[1L],
      rejection_rate = rejection_rate,
      monte_carlo_standard_error = monte_carlo_standard_error,
      repeat_number = nrow(simulation_result_group),
      significance_level = significance_level,
      lower_plot_limit = max(
        0,
        rejection_rate - 1.96 * monte_carlo_standard_error
      ),
      upper_plot_limit = min(
        1,
        rejection_rate + 1.96 * monte_carlo_standard_error
      ),
      stringsAsFactors = FALSE
    )
  }
)
rejection_rate_summary <- do.call(
  rbind,
  rejection_rate_summary_tables
)
rownames(rejection_rate_summary) <- NULL


# Save collected results ----------------------------------------------------

collected_simulation_output <- list(
  simulation_results = simulation_results,
  simulation_task_table = simulation_task_table
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
