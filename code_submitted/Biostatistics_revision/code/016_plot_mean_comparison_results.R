# File input and output -------------------------------------------------------

rejection_rate_summary_input_file <- file.path(
  "..",
  "data",
  "final",
  "014_mean_comparison_rejection_rate_summary.rds"
)
type_1_error_and_power_plot_output_file <- file.path(
  "..",
  "data",
  "final",
  "016_mean_comparison_type_1_error_and_power.pdf"
)


# Read results ---------------------------------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The ggplot2 package is required.")
}
if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("The patchwork package is required.")
}

rejection_rate_summary <- readRDS(
  rejection_rate_summary_input_file
)
rejection_rate_summary$method <- factor(
  rejection_rate_summary$method,
  levels = c("Debiased PC", "Oracle PC", "Plug-in PC")
)
rejection_rate_summary$feature_number_label <- factor(
  paste0("p = ", rejection_rate_summary$feature_number),
  levels = c("p = 100", "p = 1000")
)
significance_levels <- unique(rejection_rate_summary$significance_level)
if (length(significance_levels) != 1L) {
  stop("The rejection-rate summary must contain one significance level.")
}


# Plot helpers ---------------------------------------------------------------

create_rejection_rate_plot <- function(
    plot_data,
    plot_title,
    vertical_axis_maximum,
    vertical_axis_breaks,
    significance_reference = NULL) {
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
    ggplot2::facet_wrap(~feature_number_label, ncol = 1) +
    ggplot2::scale_x_continuous(
      breaks = sort(unique(plot_data$sample_size_per_group))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, vertical_axis_maximum),
      breaks = vertical_axis_breaks
    ) +
    ggplot2::labs(
      title = plot_title,
      x = "Sample size per group",
      y = "Rejection rate",
      color = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )

  if (!is.null(significance_reference)) {
    rejection_rate_plot <- rejection_rate_plot +
      ggplot2::geom_hline(
        yintercept = significance_reference,
        linetype = "dashed",
        color = "grey35"
      )
  }

  rejection_rate_plot
}


# Create combined figure -----------------------------------------------------

type_1_error_plot_data <- rejection_rate_summary[
  rejection_rate_summary$simulation_target == "Type I error",
  ,
  drop = FALSE
]
power_plot_data <- rejection_rate_summary[
  rejection_rate_summary$simulation_target == "Power",
  ,
  drop = FALSE
]
type_1_error_vertical_axis_maximum <- max(
  0.20,
  ceiling(max(type_1_error_plot_data$upper_plot_limit) / 0.05) * 0.05
)

type_1_error_plot <- create_rejection_rate_plot(
  plot_data = type_1_error_plot_data,
  plot_title = "Type I error",
  vertical_axis_maximum = type_1_error_vertical_axis_maximum,
  vertical_axis_breaks = seq(
    0,
    type_1_error_vertical_axis_maximum,
    by = 0.05
  ),
  significance_reference = significance_levels
)
power_plot <- create_rejection_rate_plot(
  plot_data = power_plot_data,
  plot_title = "Power",
  vertical_axis_maximum = 1,
  vertical_axis_breaks = seq(0, 1, by = 0.2)
)
type_1_error_and_power_plot <- (
  type_1_error_plot | power_plot
) +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(
  filename = type_1_error_and_power_plot_output_file,
  plot = type_1_error_and_power_plot,
  device = grDevices::pdf,
  width = 9,
  height = 7
)
