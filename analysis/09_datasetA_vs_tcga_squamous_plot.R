library(readr)
library(dplyr)
library(ggplot2)

dir.create("results/tcga/figures", recursive = TRUE, showWarnings = FALSE)

df <- read_csv(
  "results/tcga/datasetA_functional_vs_tcga_squamous.csv",
  show_col_types = FALSE
)

plot_df <- df %>%
  arrange(TCGA_Squamous_frequency_pct) %>%
  mutate(
    gene = factor(gene, levels = gene),
    datasetA_patients = factor(datasetA_patients)
  )

p <- ggplot(
  plot_df,
  aes(
    x = gene,
    y = TCGA_Squamous_frequency_pct,
    fill = datasetA_patients
  )
) +
  geom_col(width = 0.75) +
  geom_text(
    aes(
      label = ifelse(
        TCGA_Squamous_frequency_pct > 0,
        paste0(
          TCGA_Squamous_mutated,
          "/",
          TCGA_Squamous_total,
          " (",
          round(TCGA_Squamous_frequency_pct, 2),
          "%)"
        ),
        "0/236"
      )
    ),
    hjust = -0.08,
    size = 3
  ) +
  coord_flip() +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.25))
  ) +
  labs(
    title = "Concordancia entre Dataset A y TCGA-CESC escamoso",
    subtitle = "Genes HIGH/MODERATE encontrados en Dataset A y su frecuencia en 236 tumores escamosos TCGA",
    x = "Gen identificado en Dataset A",
    y = "Pacientes TCGA escamosos con mutación (%)",
    fill = "Paciente Dataset A"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  "results/tcga/figures/datasetA_vs_tcga_squamous_concordance.png",
  p,
  width = 10,
  height = 8,
  dpi = 300
)

cat("\n===== RESUMEN =====\n")

cat(
  "Genes HIGH/MODERATE Dataset A:",
  nrow(df),
  "\n"
)

cat(
  "También presentes en TCGA escamoso:",
  sum(df$TCGA_Squamous_mutated > 0),
  "\n"
)

cat(
  "Porcentaje de concordancia:",
  round(
    100 * sum(df$TCGA_Squamous_mutated > 0) / nrow(df),
    1
  ),
  "%\n"
)

cat("\nFigura creada:\n")
cat(
  "results/tcga/figures/datasetA_vs_tcga_squamous_concordance.png\n"
)
