library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

dir.create("results/qc/figures", recursive = TRUE, showWarnings = FALSE)

variant_qc <- read_csv("results/qc/datasetA_variant_qc.csv", show_col_types = FALSE)

print(variant_qc)

plot_df <- variant_qc %>%
  select(patient, mutect2_candidates, pass_variants) %>%
  pivot_longer(
    cols = c(mutect2_candidates, pass_variants),
    names_to = "metric",
    values_to = "count"
  ) %>%
  mutate(
    metric = recode(
      metric,
      mutect2_candidates = "Mutect2 candidates",
      pass_variants = "PASS variants"
    ),
    patient = factor(patient, levels = c(1094, 755, 999))
  )

p <- ggplot(plot_df, aes(x = patient, y = count, fill = metric)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = count), position = position_dodge(width = 0.9), vjust = -0.3, size = 4) +
  labs(
    title = "Dataset A: variantes somáticas por paciente",
    subtitle = "Comparación entre variantes candidatas de Mutect2 y variantes PASS",
    x = "Paciente",
    y = "Número de variantes",
    fill = "Métrica"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = "results/qc/figures/datasetA_variant_qc.png",
  plot = p,
  width = 8,
  height = 5,
  dpi = 300
)
