library(ggplot2)
library(ggpubr)
library(tidyverse)
rm(list = ls())

work_directory <- '/Users/tianyuzhang/Documents/adaptive_projection_manuscript_code/code_submitted/main_simulation/'
settings <- c('InfProject_hard', 'InfAlternative_hard')

shape_vector <- c(0, 1, 8, 9)
linetype_vector <- c(6,1,3,4)
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7",
               "#0072B2", "#D55E00", "#F0E442")
# label_vector <- c("chen2010", "debiased", "anchor", 'plug-in')

SIMULATION_BATCH <- settings[1]
num_features <- 1e3
result_file_name <- file.path(work_directory,
                              "result",
                              paste0("SETTING_", SIMULATION_BATCH,
                                     "_my_methods_",
                                     "dimension_", num_features,
                                     ".rds" ))

result_avg_p1 <- readRDS(result_file_name) |> tibble() |>
              mutate(rejection = abs(test_statistics_pc1) > 1.96) |>
              group_by(method, sample_size_1) |>
              summarize(mean_rejection = mean(rejection), .groups = "drop") |>
              filter(method %in% c('simple', 'debiased', 'lasso_pc1')) |>
              mutate(method = case_when(method == "lasso_pc1" ~ "anchor (PC1)",
                                        method == "simple"   ~ "plug-in (PC1)",
                                        method == "debiased" ~ "debiased (PC1)",
                                        TRUE ~ method))
glimpse(result_avg_p1)



SIMULATION_BATCH <- settings[2]
num_features <- 1e3
result_file_name <- file.path(work_directory,
                              "result",
                              paste0("SETTING_", SIMULATION_BATCH,
                                     "_my_methods_",
                                     "dimension_", num_features,
                                     ".rds"))
result <- readRDS(result_file_name) |> tibble() |>
          select(repeat_index, test_statistics_pc1, 
                 test_statistics_pc2, method, sample_size_1) 

result_debiase_simple <-  result |>
          filter(method %in% c("debiased")) |>
          pivot_longer(cols = c(test_statistics_pc1, test_statistics_pc2),
                       names_to = "pc_index",
                       values_to = "test_stat") |>
          mutate(pc_index = gsub("^.*_", "", pc_index)) |>
          mutate(method = paste0(method, "_", pc_index)) |>
          select(-pc_index) |>
          mutate(method = case_when(method == "debiased_pc1" ~ "debiased (PC1)",
                                    method == "debiased_pc2" ~ "debiased (PC2)"))
  

result_anchor <-  result |>
  filter(method %in% c("lasso_pc1", "lasso_pc2")) |>
  select(!any_of("test_statistics_pc2")) |>
  rename(test_stat = test_statistics_pc1) |>
  mutate(method = case_when(method == "lasso_pc1" ~ "anchor (PC1)",
                            method == "lasso_pc2" ~ "anchor (PC2)"))

result_comb <- rbind(result_debiase_simple, result_anchor) 

result_avg_p2 <- result_comb |>
  mutate(rejection = abs(test_stat) > 1.96) |>
  group_by(method, sample_size_1) |>
  summarize(mean_rejection = mean(rejection), .groups = "drop")
result_avg_p2


# --- define a single shared legend spec (order + aesthetics) ---
legend_levels <- c(
  "anchor (PC1)", "anchor (PC2)",
  "debiased (PC1)", "debiased (PC2)",
  "plug-in (PC1)"
)

col_map <- c(
  "anchor (PC1)"   = "#E69F00",
  "anchor (PC2)"   = "#D55E00",
  "debiased (PC1)" = "#56B4E9",
  "debiased (PC2)" = "#009E73",
  "plug-in (PC1)"  = "#CC79A7"
)

lt_map <- c(
  "anchor (PC1)"   = 6,
  "anchor (PC2)"   = 6,
  "debiased (PC1)" = 1,
  "debiased (PC2)" = 1,
  "plug-in (PC1)"  = 3
)

shape_map <- c(
  "anchor (PC1)"   = 16,
  "anchor (PC2)"   = 17,
  "debiased (PC1)" = 15,
  "debiased (PC2)" = 3,
  "plug-in (PC1)"  = 1
)

# ---------------- p1 ----------------
p1 <- result_avg_p1 %>%
  mutate(method = factor(method, levels = legend_levels)) %>%
  ggplot(aes(sample_size_1, mean_rejection,
             color = method, linetype = method, shape = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "#AAAAAA") +
  scale_color_manual(values = col_map, limits = legend_levels, drop = FALSE) +
  scale_linetype_manual(values = lt_map, limits = legend_levels, drop = FALSE) +
  scale_shape_manual(values = shape_map, limits = legend_levels, drop = FALSE) +
  theme_bw() +
  labs(y = "rejection %", x = "group 1 sample size", title = "Projected Null") +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.position = "bottom")

p2 <- result_avg_p2 %>%
  mutate(method = factor(method, levels = legend_levels)) %>%
  ggplot(aes(sample_size_1, mean_rejection,
             color = method, linetype = method, shape = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = col_map, limits = legend_levels, drop = FALSE) +
  scale_linetype_manual(values = lt_map, limits = legend_levels, drop = FALSE) +
  scale_shape_manual(values = shape_map, limits = legend_levels, drop = FALSE) +
  theme_bw() +
  labs(y = "rejection %", x = "group 1 sample size", title = "Alternative") +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.position = "bottom")

# 1) Make BOTH plots map linetype in the main aes (not inside geom_line)
p1 <- p1 + aes(linetype = method) + theme(legend.position = "bottom")
p2 <- p2 + aes(linetype = method) + theme(legend.position = "bottom")

# 2) Create ONE legend from a dummy plot with the FULL set of levels
legend_levels <- c(
  "anchor (PC1)", "anchor (PC2)",
  "debiased (PC1)", "debiased (PC2)",
  "plug-in (PC1)"
)

col_map <- c(
  "anchor (PC1)"   = "#E69F00",
  "anchor (PC2)"   = "#D55E00",
  "debiased (PC1)" = "#56B4E9",
  "debiased (PC2)" = "#009E73",
  "plug-in (PC1)"  = "#CC79A7"
)

lt_map <- c(
  "anchor (PC1)"   = 6,
  "anchor (PC2)"   = 6,
  "debiased (PC1)" = 1,
  "debiased (PC2)" = 1,
  "plug-in (PC1)"  = 3
)

legend_plot <- tibble(
  sample_size_1  = 1,
  mean_rejection = 0,
  method = factor(legend_levels, levels = legend_levels)
) %>%
  ggplot(aes(sample_size_1, mean_rejection,
             color = method, linetype = method, shape = method)) +
  geom_line() +
  geom_point(size = 2) +
  scale_color_manual(values = col_map, breaks = legend_levels, drop = FALSE) +
  scale_linetype_manual(values = lt_map, breaks = legend_levels, drop = FALSE) +
  scale_shape_manual(values = shape_map, breaks = legend_levels, drop = FALSE) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(
    color    = guide_legend(nrow = 1, byrow = TRUE),
    linetype = guide_legend(nrow = 1, byrow = TRUE),
    shape    = guide_legend(nrow = 1, byrow = TRUE)
  )

shared_legend <- ggpubr::get_legend(legend_plot)

# 3) Remove legends from the individual plots and stack the shared legend below
p1_noleg <- p1 + theme(legend.position = "none")
p2_noleg <- p2 + theme(legend.position = "none")

inf_plot <- ggarrange(
  ggarrange(p1_noleg, p2_noleg, ncol = 2, labels = c("A", "B")),
  shared_legend,
  ncol = 1,
  heights = c(1, 0.18)
)

inf_plot

main_simulation_plot_file <- file.path(work_directory,
                                    "result", '23_main_simulation.pdf')
ggsave(file = main_simulation_plot_file,
       width = 200, height = 100, units = "mm")
