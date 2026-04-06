#!/usr/bin/env Rscript
# Generate dataset_characterization.png (3-panel: barrier dist, RMSD density, cost vs barrier)
# Usage: pixi run -e brms Rscript scripts/gen_dataset_char.R

suppressPackageStartupMessages({
  library(bayescomp)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

DATA_PATH <- "data/baker_bench.csv"
SUPPL <- "imgs/gen/R/suppl"
dir.create(SUPPL, recursive = TRUE, showWarnings = FALSE)

method_colors <- c(CINEB = "#FF655D", OCINEB = "#004D40")

bench_raw <- bc_read_benchmark(
  DATA_PATH, format = "wide",
  method_suffixes = c("CINEB", "MMF"),
  system_col = "System", count_col = "Calls",
  time_col = "Time", success_col = "Term"
)

bench_long <- bc_pivot_long(
  bench_raw, method_pattern = "_(CINEB|MMF)$",
  method_levels = c("CINEB", "MMF"),
  cens_value = "BAD_MAX_ITERATIONS"
)

bench_long$Method <- factor(
  ifelse(bench_long$method == "MMF", "OCINEB", as.character(bench_long$method)),
  levels = c("CINEB", "OCINEB")
)

dp <- bench_long %>% filter(success == TRUE)

# Panel A: Barrier height distribution
p_barrier <- ggplot(dp, aes(x = Barrier, fill = Method)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 15, color = "white") +
  scale_fill_manual(values = method_colors) +
  labs(title = "A. Reaction Barrier Distribution",
       x = "Barrier Height (eV)", y = "Count") +
  theme_bayescomp(base_size = 36) +
  theme(legend.position = "none")

# Panel B: RMSD saddle accuracy
p_accuracy <- ggplot(dp, aes(x = RMSD_Saddle, fill = Method)) +
  geom_density(alpha = 0.5, color = NA) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = method_colors) +
  labs(title = "B. Saddle Point Similarity",
       x = expression("RMSD to Benchmark Saddle (" * ring(A) * ")"),
       y = "Density") +
  theme_bayescomp(base_size = 36) +
  theme(legend.position = "none")

# Panel C: Cost vs Barrier
p_correlation <- ggplot(dp, aes(x = Barrier, y = count, color = Method)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", alpha = 0.5) +
  scale_y_log10() +
  scale_color_manual(values = method_colors) +
  labs(title = "C. Scaling with Barrier Height",
       x = "Barrier Height (eV)",
       y = "Gradient Calls (log10)") +
  theme_bayescomp(base_size = 36) +
  theme(legend.position = c(0.8, 0.2))

p_combined <- (p_barrier | p_accuracy) / p_correlation
ggsave(file.path(SUPPL, "dataset_characterization.png"),
       p_combined, width = 14, height = 12, dpi = 300)
cat("OK: dataset_characterization.png\n")
