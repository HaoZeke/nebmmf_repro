#!/usr/bin/env Rscript
# Regenerate the v1-style dumbbell plot with updated data
# Run: pixi run -e brms Rscript scripts/regen_dumbbell.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(showtext)
  library(sysfonts)
  library(ggthemes)
})

# Colors
ruhi_colors <- c(
  teal = "#004D40", coral = "#FF655D", sky = "#1E88E5",
  sunshine = "#F1DB4B", green = "#009E73"
)

# Theme
font_add_google("Atkinson Hyperlegible", "Atkinson")
showtext_auto()

theme_Publication <- function(base_size = 22, base_family = "Atkinson") {
  theme_foundation(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
      panel.background = element_rect(colour = NA, fill = "#FFFFFF"),
      plot.background = element_rect(colour = NA, fill = "#FFFFFF"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.title = element_text(face = "bold", size = rel(1.0)),
      axis.text = element_text(size = rel(0.9)),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_line(colour = "#e6e3dd"),
      panel.grid.minor = element_blank(),
      legend.background = element_rect(fill = "#FFFFFF", colour = NA),
      legend.key = element_rect(colour = NA, fill = "#FFFFFF"),
      legend.position = "right",
      strip.background = element_rect(colour = "#FFFFFF", fill = "#FFFFFF"),
      strip.text = element_text(face = "bold", size = rel(1.0))
    )
}

# Read data
data_raw <- read.csv("data/baker_bench.csv", stringsAsFactors = FALSE)
colnames(data_raw) <- gsub("_MMF$", "_OCINEB", colnames(data_raw))

PAPER_IMGS <- "imgs/gen/R"
dir.create(PAPER_IMGS, recursive = TRUE, showWarnings = FALSE)

# Prepare dumbbell data
df_dumbbell <- data_raw %>%
  mutate(
    Calls_CINEB_Plot = ifelse(
      Term_CINEB == "BAD_MAX_ITERATIONS", 8000, Calls_CINEB),
    Calls_OCINEB_Plot = Calls_OCINEB,
    Time_CINEB_Plot = Time_CINEB,
    Time_OCINEB_Plot = Time_OCINEB
  ) %>%
  pivot_longer(
    cols = c(Calls_CINEB_Plot, Calls_OCINEB_Plot,
             Time_CINEB_Plot, Time_OCINEB_Plot),
    names_to = c("Metric", "Method"),
    names_pattern = "(.*)_(.*)_Plot",
    values_to = "Value"
  ) %>%
  mutate(
    Method = factor(Method, levels = c("CINEB", "OCINEB")),
    Metric_Label = case_when(
      Metric == "Calls" ~ "Gradient Evaluations",
      Metric == "Time" ~ "Wall Time (s)"
    ),
    Status_Point = case_when(
      Method == "CINEB" & Term_CINEB == "BAD_MAX_ITERATIONS" ~ "Failed",
      TRUE ~ "Converged"
    )
  )

# Average and median rows
df_avg <- df_dumbbell %>%
  group_by(Metric_Label, Method) %>%
  summarise(
    Value = mean(Value, na.rm = TRUE),
    Status_Point = "Converged",
    Reaction = "bold('AVERAGE CHANGE')",
    .groups = "drop"
  )

df_median <- df_dumbbell %>%
  group_by(Metric_Label, Method) %>%
  summarise(
    Value = median(Value, na.rm = TRUE),
    Status_Point = "Converged",
    Reaction = "bold('MEDIAN CHANGE')",
    .groups = "drop"
  )

df_combined <- bind_rows(df_dumbbell, df_avg, df_median)

# Segments and differences
df_segments <- df_combined %>%
  pivot_wider(names_from = Method, values_from = c(Value, Status_Point)) %>%
  filter(!is.na(Value_CINEB) & !is.na(Value_OCINEB)) %>%
  group_by(Metric_Label) %>%
  mutate(
    Improvement = Value_CINEB - Value_OCINEB,
    Diff_Value = Value_OCINEB - Value_CINEB,
    Diff_Label = case_when(
      Metric_Label == "Gradient Evaluations" ~
        sprintf("%+d", as.integer(Diff_Value)),
      Metric_Label == "Wall Time (s)" ~ sprintf("%+.1fs", Diff_Value)
    ),
    Right_Pos = max(c(Value_CINEB, Value_OCINEB), na.rm = TRUE) * 1.8
  ) %>%
  ungroup()

# Sort by improvement in gradient evaluations
rank_order <- df_segments %>%
  filter(Metric_Label == "Gradient Evaluations") %>%
  arrange(Improvement) %>%
  pull(Reaction)

df_combined$Reaction <- factor(df_combined$Reaction, levels = rank_order)
df_segments$Reaction <- factor(df_segments$Reaction, levels = rank_order)

# Plot
p_dumbbell <- ggplot() +
  geom_segment(
    data = df_segments,
    aes(y = Reaction, yend = Reaction, x = Value_OCINEB, xend = Value_CINEB),
    color = "grey70", linewidth = 0.8
  ) +
  geom_text(
    data = df_segments,
    aes(y = Reaction, x = Right_Pos, label = Diff_Label),
    vjust = 0, size = 14, color = "grey40", family = "Atkinson"
  ) +
  geom_point(
    data = df_combined,
    aes(y = Reaction, x = Value, color = Method, shape = Status_Point),
    size = 7, alpha = 0.9
  ) +
  facet_wrap(~Metric_Label, scales = "free_x", ncol = 2) +
  scale_y_discrete(
    labels = scales::label_parse(),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  scale_color_manual(
    values = c("CINEB" = ruhi_colors[["coral"]],
               "OCINEB" = ruhi_colors[["teal"]])
  ) +
  scale_shape_manual(values = c("Converged" = 19, "Failed" = 4)) +
  scale_x_log10(
    labels = label_number(scale_cut = cut_short_scale()),
    expand = expansion(mult = c(0.1, 0.2))
  ) +
  labs(x = "Log Scale Value", y = NULL,
       color = "Method", shape = "Status") +
  theme_Publication(base_size = 48) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "top",
    legend.box = "horizontal",
    legend.direction = "horizontal"
  )

ggsave(
  file.path(PAPER_IMGS, "dumbbell_plot.png"),
  p_dumbbell, width = 20, height = 14, dpi = 300
)
cat("OK: dumbbell_plot.png\n")
