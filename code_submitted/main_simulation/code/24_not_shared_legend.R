library(ggplot2)
library(ggpubr)
library(tidyverse)
library(ggtext)

rm(list = ls())

work_directory <- '/Users/tianyuzhang/Documents/adaptive_projection_manuscript_code/code_submitted/main_simulation/'
settings <- c('InfProject_hard', 'InfAlternative_hard')

# --- legend spec (shared aesthetics, but legends remain per-panel) ---
legend_levels <- c(
  "T<sub>anc</sub>(v<sub>1</sub>)",
  "T<sub>anc</sub>(v<sub>2</sub>)",
  "T<sub>1s</sub>(v<sub>1</sub>)",
  "T<sub>1s</sub>(v<sub>2</sub>)",
  "T<sub>pi</sub>(v<sub>1</sub>)"
)

col_map <- c(
  "T<sub>anc</sub>(v<sub>1</sub>)" = "#E69F00",
  "T<sub>anc</sub>(v<sub>2</sub>)"  = "#D55E00",
  "T<sub>1s</sub>(v<sub>1</sub>)"               = "#56B4E9",
  "T<sub>1s</sub>(v<sub>2</sub>)"               = "#009E73",
  "T<sub>pi</sub>(v<sub>1</sub>)"                = "#CC79A7"
)

lt_map <- c(
  "T<sub>anc</sub>(v<sub>1</sub>)" = 6,
  "T<sub>anc</sub>(v<sub>2</sub>)"                 = 6,
  "T<sub>1s</sub>(v<sub>1</sub>)"               = 1,
  "T<sub>1s</sub>(v<sub>2</sub>)"               = 1,
  "T<sub>pi</sub>(v<sub>1</sub>)"                = 3
)

shape_map <- c(
  "T<sub>anc</sub>(v<sub>1</sub>)" = 16,
  "T<sub>anc</sub>(v<sub>2</sub>)"                 = 17,
  "T<sub>1s</sub>(v<sub>1</sub>)"               = 15,
  "T<sub>1s</sub>(v<sub>2</sub>)"               = 3,
  "T<sub>pi</sub>(v<sub>1</sub>)"                = 1
)

# =======================
# Panel 1: Projected Null
# =======================
SIMULATION_BATCH <- settings[1]
num_features <- 1e3
result_file_name <- file.path(
  work_directory, "result",
  paste0("SETTING_", SIMULATION_BATCH, "_my_methods_", "dimension_", num_features, ".rds")
)

result_avg_p1 <- readRDS(result_file_name) |>
  tibble() |>
  mutate(rejection = abs(test_statistics_pc1) > 1.96) |>
  group_by(method, sample_size_1) |>
  summarize(mean_rejection = mean(rejection), .groups = "drop") |>
  filter(method %in% c('simple', 'debiased', 'lasso_pc1')) |>
  mutate(method = case_when(
    method == "lasso_pc1" ~ "T<sub>anc</sub>(v<sub>1</sub>)",
    method == "simple"   ~ "T<sub>pi</sub>(v<sub>1</sub>)",
    method == "debiased" ~ "T<sub>1s</sub>(v<sub>1</sub>)",
    TRUE ~ method
  )) |>
  mutate(method = factor(method, levels = legend_levels))

p1 <- ggplot(
  result_avg_p1,
  aes(sample_size_1, mean_rejection, color = method, linetype = method, shape = method)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "#AAAAAA") +
  scale_color_manual(values = col_map, breaks = legend_levels) +
  scale_linetype_manual(values = lt_map, breaks = legend_levels) +
  scale_shape_manual(values = shape_map, breaks = legend_levels) +
  theme_bw() +
  labs(y = "rejection %", x = "group 1 sample size", title = "Projected Null") +
  theme(plot.title = element_text(size = 15, face = "bold"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 12))

# =======================
# Panel 2: Alternative
# =======================
SIMULATION_BATCH <- settings[2]
result_file_name <- file.path(
  work_directory, "result",
  paste0("SETTING_", SIMULATION_BATCH, "_my_methods_", "dimension_", num_features, ".rds")
)

result <- readRDS(result_file_name) |>
  tibble() |>
  select(repeat_index, test_statistics_pc1, test_statistics_pc2, method, sample_size_1)

result_debiased <- result |>
  filter(method %in% c("debiased")) |>
  pivot_longer(cols = c(test_statistics_pc1, test_statistics_pc2),
               names_to = "pc_index", values_to = "test_stat") |>
  mutate(pc_index = gsub("^.*_", "", pc_index),
         method = paste0(method, "_", pc_index)) |>
  select(-pc_index) |>
  mutate(method = case_when(
    method == "debiased_pc1" ~ "T<sub>1s</sub>(v<sub>1</sub>)",
    method == "debiased_pc2" ~ "T<sub>1s</sub>(v<sub>2</sub>)",
    TRUE ~ method
  ))

result_anchor <- result |>
  filter(method %in% c("lasso_pc1", "lasso_pc2")) |>
  transmute(
    repeat_index,
    sample_size_1,
    test_stat = test_statistics_pc1,
    method = case_when(
      method == "lasso_pc1" ~ "T<sub>anc</sub>(v<sub>1</sub>)",
      method == "lasso_pc2" ~ "T<sub>anc</sub>(v<sub>2</sub>)"
    )
  )

result_comb <- bind_rows(result_debiased, result_anchor)

result_avg_p2 <- result_comb |>
  mutate(rejection = abs(test_stat) > 1.96) |>
  group_by(method, sample_size_1) |>
  summarize(mean_rejection = mean(rejection), .groups = "drop") |>
  mutate(method = factor(method, levels = legend_levels))

p2 <- ggplot(
  result_avg_p2,
  aes(sample_size_1, mean_rejection, color = method, linetype = method, shape = method)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = col_map, breaks = legend_levels) +
  scale_linetype_manual(values = lt_map, breaks = legend_levels) +
  scale_shape_manual(values = shape_map, breaks = legend_levels) +
  theme_bw() +
  labs(y = "rejection %", x = "group 1 sample size", title = "Alternative") +
  theme(plot.title = element_text(size = 15, face = "bold"),
        legend.position = "bottom",
        legend.text = element_markdown(size = 12))

# =======================
# Arrange: each subplot keeps its own legend
# =======================
inf_plot <- ggarrange(
  p1, p2,
  ncol = 2,
  labels = c("A", "B"),
  common.legend = FALSE
)

inf_plot

main_simulation_plot_file <- file.path(work_directory, "result", "24_main_simulation.pdf")
ggsave(main_simulation_plot_file, plot = inf_plot,
       width = 250, height = 125, units = "mm")


