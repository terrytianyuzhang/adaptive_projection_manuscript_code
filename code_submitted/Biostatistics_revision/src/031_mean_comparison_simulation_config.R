# File input ----------------------------------------------------------------

original_simulation_config_source_file <- file.path(
  "..",
  "src",
  "011_mean_comparison_simulation_config.R"
)


source(
  original_simulation_config_source_file,
  local = TRUE
)


# Mean allocation -----------------------------------------------------------

original_projected_null_mean_difference <- projected_null_mean_difference
enlarged_projected_null_mean_difference <- (
  5 * original_projected_null_mean_difference
)

if (setting_id %in% c(
  1L, 2L, 3L, 4L, 5L,
  6L, 7L, 8L, 9L, 10L
)) {
  treatment_minus_control_mean <- enlarged_projected_null_mean_difference
} else if (setting_id %in% c(
  11L, 12L, 13L, 14L, 15L
)) {
  treatment_minus_control_mean <- (
    original_projected_null_mean_difference +
      projected_alternative_addition
  )
} else {
  treatment_minus_control_mean <- (
    original_projected_null_mean_difference +
      8 * projected_alternative_addition
  )
}

control_mean <- -treatment_minus_control_mean / 2
treatment_mean <- treatment_minus_control_mean / 2

simulation_config$control_mean <- control_mean
simulation_config$treatment_mean <- treatment_mean
