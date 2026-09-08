# File output ---------------------------------------------------------------

simulation_task_table_output_file <- file.path(
  "..",
  "data",
  "intermediate",
  "011_mean_comparison_task_table.csv"
)


# Task grid -----------------------------------------------------------------

method_names <- c(
  "debiased_pc",
  "oracle_pc",
  "plug_in_pc"
)
setting_ids <- seq_len(20L)
batch_numbers <- seq_len(100L)

task_grid <- expand.grid(
  batch_number = batch_numbers,
  method = method_names,
  setting_id = setting_ids,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

simulation_task_table <- data.frame(
  id = seq_len(nrow(task_grid)),
  method = task_grid$method,
  setting_id = task_grid$setting_id,
  batch_number = task_grid$batch_number,
  stringsAsFactors = FALSE
)


# Save task table -----------------------------------------------------------

write.csv(
  simulation_task_table,
  file = simulation_task_table_output_file,
  row.names = FALSE,
  quote = FALSE
)
