# Simulation grid ------------------------------------------------------------

simulation_setting_ids <- c(
  "p0100_n0050", "p0100_n0200", "p0100_n0400",
  "p0100_n0600", "p0100_n0800", "p1000_n0050",
  "p1000_n0200", "p1000_n0400", "p1000_n0600",
  "p1000_n0800"
)

# Assign every setting attribute from its setting ID.  These vectors are the
# only place where the simulation grid is specified.
feature_100_setting_ids <- c(
  "p0100_n0050", "p0100_n0200", "p0100_n0400",
  "p0100_n0600", "p0100_n0800"
)
feature_1000_setting_ids <- c(
  "p1000_n0050", "p1000_n0200", "p1000_n0400",
  "p1000_n0600", "p1000_n0800"
)
sample_size_50_setting_ids <- c("p0100_n0050", "p1000_n0050")
sample_size_200_setting_ids <- c("p0100_n0200", "p1000_n0200")
sample_size_400_setting_ids <- c("p0100_n0400", "p1000_n0400")
sample_size_600_setting_ids <- c("p0100_n0600", "p1000_n0600")
sample_size_800_setting_ids <- c("p0100_n0800", "p1000_n0800")
projected_mean_setting_ids <- simulation_setting_ids

simulation_setting_table <- data.frame(
  setting = simulation_setting_ids,
  feature_number = NA_integer_,
  sample_size_per_group = NA_integer_,
  population_mean_structure = NA_character_,
  stringsAsFactors = FALSE
)
simulation_setting_table$feature_number[
  simulation_setting_table$setting %in% feature_100_setting_ids
] <- 100L
simulation_setting_table$feature_number[
  simulation_setting_table$setting %in% feature_1000_setting_ids
] <- 1000L
simulation_setting_table$sample_size_per_group[
  simulation_setting_table$setting %in% sample_size_50_setting_ids
] <- 50L
simulation_setting_table$sample_size_per_group[
  simulation_setting_table$setting %in% sample_size_200_setting_ids
] <- 200L
simulation_setting_table$sample_size_per_group[
  simulation_setting_table$setting %in% sample_size_400_setting_ids
] <- 400L
simulation_setting_table$sample_size_per_group[
  simulation_setting_table$setting %in% sample_size_600_setting_ids
] <- 600L
simulation_setting_table$sample_size_per_group[
  simulation_setting_table$setting %in% sample_size_800_setting_ids
] <- 800L
simulation_setting_table$population_mean_structure[
  simulation_setting_table$setting %in% projected_mean_setting_ids
] <- "projected_null_and_projected_alternative"

if (anyNA(simulation_setting_table)) {
  stop("Every simulation setting ID must specify dimension, sample size, and population means.")
}
feature_numbers <- sort(unique(simulation_setting_table$feature_number))
simulation_repeat_number <- 1000L
repeat_number_per_batch <- 10L


# Mean-comparison settings ---------------------------------------------------

cross_fitting_fold_number <- 5L
threshold_multiplier <- 1.25
mean_difference_threshold_multiplier <- 1.25
projected_null_mean_difference_size <- 1.00
projected_null_support_size <- 10L
projected_alternative_mean_difference_size <- 0.25
projected_alternative_support_size <- 10L
significance_level <- 0.05
mean_comparison_simulation_seed <- 20260829L

simulation_method_ids <- c(
  "Debiased PC" = "debiased_pc",
  "Oracle PC" = "oracle_pc",
  "Plug-in PC" = "plug_in_pc"
)
