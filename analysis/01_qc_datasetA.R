library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================
# BIOINFORMÁTICA I - ENTREGA 3
# QC Dataset A: E-MTAB-11407 / PRJEB50605
# ============================================================

# 1. Cargar datos -----------------------------------------------------------

qc <- read_csv(
  "results/qc/datasetA_alignment_qc.csv",
  show_col_types = FALSE
)

qc <- qc %>%
  mutate(
    sample_type = recode(sample,
                         "N" = "Normal",
                         "T" = "Tumor"),
    sample_id = paste0(patient, "_", sample)
  )

print(qc)


# 2. Preparar porcentajes de alineamiento ----------------------------------

qc_alignment <- qc %>%
  select(
    sample_id,
    pct_mapped,
    pct_properly_paired
  ) %>%
  pivot_longer(
    cols = c(pct_mapped, pct_properly_paired),
    names_to = "metric",
    values_to = "percentage"
  ) %>%
  mutate(
    metric = recode(
      metric,
      "pct_mapped" = "Lecturas mapeadas",
      "pct_properly_paired" = "Correctamente pareadas"
    )
  )


# 3. Figura QC: alineamiento ------------------------------------------------

p_alignment <- ggplot(
  qc_alignment,
  aes(
    x = sample_id,
    y = percentage,
    fill = metric
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_text(
    aes(label = sprintf("%.2f%%", percentage)),
    position = position_dodge(width = 0.8),
    vjust = -0.4,
    size = 3.5
  ) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, 20)
  ) +
  labs(
    title = "Control de calidad del alineamiento",
    subtitle = "Dataset A - subconjuntos de 1 millón de pares de lecturas",
    x = "Muestra",
    y = "Porcentaje (%)",
    fill = "Métrica"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

print(p_alignment)

ggsave(
  "results/qc/figures/datasetA_alignment_qc.png",
  p_alignment,
  width = 9,
  height = 6,
  dpi = 300
)


# 4. Figura QC: duplicación -------------------------------------------------

p_duplication <- ggplot(
  qc,
  aes(
    x = sample_id,
    y = pct_duplication,
    fill = sample_type
  )
) +
  geom_col(width = 0.65) +
  geom_text(
    aes(label = sprintf("%.2f%%", pct_duplication)),
    vjust = -0.4,
    size = 3.8
  ) +
  scale_y_continuous(
    limits = c(0, 60),
    breaks = seq(0, 60, 10)
  ) +
  labs(
    title = "Porcentaje de duplicación por muestra",
    subtitle = "Métricas obtenidas mediante MarkDuplicates",
    x = "Muestra",
    y = "Duplicación (%)",
    fill = "Tipo de muestra"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

print(p_duplication)

ggsave(
  "results/qc/figures/datasetA_duplication_qc.png",
  p_duplication,
  width = 9,
  height = 6,
  dpi = 300
)
