# File input and output -------------------------------------------------------

mean_comparison_simulation_input_file <- file.path(
  data_directory,
  "intermediate",
  "003_mean_comparison_simulation_p100.rds"
)

mean_comparison_summary_output_file <- file.path(
  data_directory,
  "final",
  "004_mean_comparison_rejection_rate_summary_p100.rds"
)

type_1_error_plot_output_file <- file.path(
  data_directory,
  "final",
  "004_mean_comparison_type_1_error_p100.pdf"
)

power_plot_output_file <- file.path(
  data_directory,
  "final",
  "004_mean_comparison_power_p100.pdf"
)


# Read and summarize ---------------------------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The ggplot2 package is required.")
}

mean_comparison_simulation_output <- readRDS(
  mean_comparison_simulation_input_file
)
mean_comparison_simulation_results <- (
  mean_comparison_simulation_output$simulation_results
)

rejection_rate_summary <- aggregate(
  rejected ~ sample_size_per_group + test_family +
    simulation_target + method,
  data = mean_comparison_simulation_results,
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
  sample_size_per_group = rejection_rate_summary$sample_size_per_group,
  test_family = rejection_rate_summary$test_family,
  simulation_target = rejection_rate_summary$simulation_target,
  method = rejection_rate_summary$method,
  rejection_rate = rejection_rate_summary$rejected[, "rejection_rate"],
  monte_carlo_standard_error = (
    rejection_rate_summary$rejected[, "monte_carlo_standard_error"]
  ),
  repeat_number = rejection_rate_summary$rejected[, "repeat_number"]
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

dir.create(
  dirname(mean_comparison_summary_output_file),
  recursive = TRUE,
  showWarnings = FALSE
)
saveRDS(
  rejection_rate_summary,
  file = mean_comparison_summary_output_file
)


# Plot -----------------------------------------------------------------------

create_rejection_rate_plot <- function(
    plot_data,
    vertical_axis_label,
    vertical_axis_maximum,
    vertical_axis_breaks,
    include_significance_reference = FALSE) {
  rejection_rate_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = sample_size_per_group,
      y = rejection_rate,
      color = method,
      group = method
    )
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = lower_plot_limit,
        ymax = upper_plot_limit
      ),
      width = 20,
      linewidth = 0.45
    ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~test_family) +
    ggplot2::scale_x_continuous(
      breaks = sort(unique(plot_data$sample_size_per_group))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, vertical_axis_maximum),
      breaks = vertical_axis_breaks
    ) +
    ggplot2::labs(
      x = "Sample size per group",
      y = vertical_axis_label,
      color = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(
        hjust = 0.5,
        margin = ggplot2::margin(t = 8)
      )
    )

  if (include_significance_reference) {
    rejection_rate_plot <- rejection_rate_plot +
      ggplot2::geom_hline(
        yintercept = mean_comparison_simulation_output$settings$
          significance_level,
        linetype = "dashed",
        color = "grey35"
      )
  }

  rejection_rate_plot
}

type_1_error_plot <- create_rejection_rate_plot(
  plot_data = rejection_rate_summary[
    rejection_rate_summary$simulation_target == "Type I error",
    ,
    drop = FALSE
  ],
  vertical_axis_label = "Type I error",
  vertical_axis_maximum = 0.20,
  vertical_axis_breaks = seq(0, 0.20, by = 0.05),
  include_significance_reference = TRUE
)
power_plot <- create_rejection_rate_plot(
  plot_data = rejection_rate_summary[
    rejection_rate_summary$simulation_target == "Power",
    ,
    drop = FALSE
  ],
  vertical_axis_label = "Power",
  vertical_axis_maximum = 1,
  vertical_axis_breaks = seq(0, 1, by = 0.2)
)

ggplot2::ggsave(
  filename = type_1_error_plot_output_file,
  plot = type_1_error_plot,
  device = grDevices::pdf,
  width = 7.5,
  height = 4.2
)
ggplot2::ggsave(
  filename = power_plot_output_file,
  plot = power_plot,
  device = grDevices::pdf,
  width = 7.5,
  height = 4.2
)

print(rejection_rate_summary)
