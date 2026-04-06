#!/usr/bin/env Rscript
# Rerun brms Bayesian analysis using bayescomp package
# Run: pixi run -e brms Rscript scripts/brms_reanalysis.R
#
# Reads baker_bench.csv, fits negbinomial model, generates 6 supplementary
# figures + conditional effects plot for the main paper.

suppressPackageStartupMessages({
  library(brms)
  library(tidybayes)
  library(bayesplot)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(posterior)
  library(showtext)
  library(sysfonts)
  library(ggthemes)
  library(archive)
})

# Load bayescomp from local install
devtools::load_all("~/Git/Github/Rlang/bayescomp")

# --- Configuration ---
DATA_PATH <- "data/baker_bench.csv"
MODEL_DIR <- "data/models"
PAPER_IMGS <- "imgs/gen/R"
SUPPL_DIR <- file.path(PAPER_IMGS, "suppl")

dir.create(MODEL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPPL_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Data ---
cat("Reading benchmark data...\n")
bench_raw <- bc_read_benchmark(
  DATA_PATH,
  format = "wide",
  method_suffixes = c("CINEB", "MMF"),
  system_col = "System",
  count_col = "Calls",
  time_col = "Time",
  success_col = "Term"
)

bench_long <- bc_pivot_long(
  bench_raw,
  method_pattern = "_(CINEB|MMF)$",
  method_levels = c("CINEB", "MMF"),
  cens_value = "BAD_MAX_ITERATIONS"
)

# Rename MMF -> OCINEB for paper consistency
bench_long$method <- factor(
  ifelse(bench_long$method == "MMF", "OCINEB", as.character(bench_long$method)),
  levels = c("CINEB", "OCINEB")
)

# Verify data
cat(sprintf("Data: %d observations, %d systems, %d methods\n",
    nrow(bench_long), nlevels(bench_long$system_id), nlevels(bench_long$method)))

# Summary stats
bench_long %>%
  group_by(method) %>%
  summarize(n = n(), mean = mean(count), median = median(count)) %>%
  print()

# --- Model ---
cat("\nFitting model...\n")
model <- bc_fit(
  bench_long,
  response = "count",
  spline_by_method = "RMSD_Init_Final",
  model_shape = TRUE,
  chains = 8, iter = 5000, warmup = 2000,
  cores = 8, seed = 1995,
  adapt_delta = 0.999, max_treedepth = 15,
  file = file.path(MODEL_DIR, "bayescomp_ocineb_v2")
)

# --- Effects ---
cat("\nExtracting effects...\n")
effects <- bc_summarize_effects(model)
effect_tbl <- bc_effect_table(model)
cat("\nEffect summary:\n")
print(effect_tbl)

# Save effect table for paper reference
writeLines(
  capture.output(print(effect_tbl)),
  file.path(MODEL_DIR, "effect_summary.txt")
)

# Save full model summary for supplementary verbatim block
writeLines(
  capture.output({
    cat("=== brms Model Summary ===\n\n")
    print(summary(model))
    cat("\n=== LOO Cross-Validation ===\n\n")
    print(loo(model))
    cat("\n=== Effect Table ===\n\n")
    print(effect_tbl)
  }),
  file.path(MODEL_DIR, "brms_model_output.txt")
)

# --- Figures ---
cat("\nGenerating figures...\n")

# Colors matching the paper
method_colors <- c("CINEB" = "#FF655D", "OCINEB" = "#004D40")

# 1. Distance spread histogram
cat("  [1/7] nolog_dist.png\n")
p1 <- ggplot(bench_long, aes(x = RMSD_Init_Final)) +
  geom_histogram(bins = 20, fill = "gray50", color = "white") +
  labs(x = "RMSD(Init, TS) (linear scale)", y = "Count") +
  theme_bayescomp(base_size = 36)
ggsave(file.path(SUPPL_DIR, "nolog_dist.png"), p1, width = 12, height = 8, dpi = 300)

# 2. Posterior predictive density
cat("  [2/7] pp_density.png\n")
p2 <- bc_plot_pp(model, ndraws = 50) +
  labs(title = "Posterior Predictive Check: Density")
ggsave(file.path(SUPPL_DIR, "pp_density.png"), p2, width = 12, height = 8, dpi = 300)

# 3. Posterior predictive group intervals
cat("  [3/7] pp_group.png\n")
pp_plots <- bc_pp_check(model, group_col = "method")
p3 <- pp_plots$grouped +
  labs(title = "Posterior Predictive Check: Intervals by Method") +
  theme_bayescomp(base_size = 36)
ggsave(file.path(SUPPL_DIR, "pp_group.png"), p3, width = 12, height = 8, dpi = 300)

# 4. Shape posterior
cat("  [4/7] brms_shape_posterior.png\n")
p4 <- bc_plot_shape(model, method_col = "method", colors = method_colors)
ggsave(file.path(SUPPL_DIR, "brms_shape_posterior.png"), p4,
       width = 12, height = 8, dpi = 300)

# 5. LOO-PIT
cat("  [5/7] ppc_loo.png\n")
loo_result <- bc_loo(model)
p5 <- bc_plot_loo_pit(model, loo_result, bench_long$count)
ggsave(file.path(SUPPL_DIR, "ppc_loo.png"), p5, width = 12, height = 8, dpi = 300)

# 6. Conditional effects (main paper figure)
cat("  [6/7] brms_pes.png\n")
cond_effects <- conditional_effects(model,
  effects = "RMSD_Init_Final:method", points = TRUE)
plot_data <- cond_effects$`RMSD_Init_Final:method`

p6 <- ggplot(plot_data, aes(x = RMSD_Init_Final, y = estimate__,
                             color = method, fill = method)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.5) +
  geom_jitter(data = bench_long, aes(y = count),
              width = 0.05, height = 0, size = 3, alpha = 0.5, shape = 16) +
  scale_y_log10() +
  scale_color_manual(values = method_colors, labels = c("CI-NEB", "OCI-NEB")) +
  scale_fill_manual(values = method_colors, labels = c("CI-NEB", "OCI-NEB")) +
  labs(
    x = expression("Distance (" * RMSD[Init - TS] * " in " * ring(A) * ")"),
    y = "Gradient Calls (Log Scale)",
    color = "Method", fill = "Method"
  ) +
  theme_bayescomp(base_size = 36) +
  theme(legend.position = "inside", legend.position.inside = c(0.15, 0.85),
        legend.title = element_text(face = "bold", size = rel(1.2)),
        legend.text = element_text(size = rel(1.0)))
ggsave(file.path(PAPER_IMGS, "brms_pes.png"), p6, width = 12, height = 8, dpi = 300)

# 7. Cactus + Violin (v1 style, 2-panel)
cat("  [7/7] cactus_violin.png\n")
library(patchwork)
cactus_data <- bench_long %>%
  filter(success == TRUE) %>%
  group_by(method) %>%
  arrange(time) %>%
  mutate(solved_count = row_number()) %>%
  ungroup()

panel_a <- ggplot(cactus_data, aes(x = time, y = solved_count, color = method)) +
  geom_step(linewidth = 1.5) +
  scale_x_log10(breaks = c(1, 10, 100, 500), labels = scales::comma) +
  scale_color_manual(values = method_colors, labels = c("CI-NEB", "OCI-NEB")) +
  labs(title = "A", x = "Time (s, log10)", y = "Cumulative Problems Solved") +
  theme_bayescomp(base_size = 36) +
  theme(legend.position = c(0.8, 0.2))

panel_b <- bench_long %>%
  filter(success == TRUE) %>%
  ggplot(aes(x = method, y = count, fill = method)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.15, alpha = 0.5, fill = "white") +
  scale_y_log10(breaks = c(10, 30, 100, 300, 1000, 3000), labels = scales::comma) +
  scale_fill_manual(values = method_colors, labels = c("CI-NEB", "OCI-NEB")) +
  labs(title = "B", x = "Method", y = "Grad. Evals (log10)") +
  theme_bayescomp(base_size = 36) +
  theme(legend.position = "none")

combined <- panel_a + panel_b + plot_layout(ncol = 2, widths = c(1.2, 0.8))
ggsave(file.path(SUPPL_DIR, "cactus_violin.png"), combined,
       width = 16, height = 10, dpi = 300)

cat("\nAll figures generated.\n")
cat(sprintf("Supplementary: %s\n", SUPPL_DIR))
cat(sprintf("Main paper: %s\n", PAPER_IMGS))

# Print key numbers for paper text update
cat("\n=== Numbers for paper text ===\n")
print(effect_tbl)
