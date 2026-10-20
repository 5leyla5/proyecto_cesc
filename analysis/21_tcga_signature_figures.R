library(readr)
library(dplyr)
library(ggplot2)

# ==========================================================
# Datos
# ==========================================================

dat <- read_csv(
  "results/tcga/COSMIC_SBS_with_metadata.csv",
  show_col_types = FALSE
)

dir.create(
  "results/tcga/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

# Cohorte principal de alta confianza
hc <- dat %>%
  filter(
    total_snvs >= 50,
    cosine >= 0.90
  )

# ==========================================================
# FIGURA 1
# APOBEC G2 vs G3
# ==========================================================

grade <- hc %>%
  filter(
    tumor_grade %in% c("G2", "G3")
  ) %>%
  mutate(
    APOBEC_pct = APOBEC * 100,
    tumor_grade = factor(
      tumor_grade,
      levels = c("G2", "G3")
    )
  )

p1 <- ggplot(
  grade,
  aes(
    x = tumor_grade,
    y = APOBEC_pct
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.45,
    size = 1.5
  ) +
  labs(
    title = "Actividad APOBEC según grado tumoral en TCGA-CESC",
    subtitle = paste0(
      "SBS2 + SBS13 · alta confianza (≥50 SNV; cosine ≥0,90)\n",
      "Wilcoxon p = 0,0358 · FDR = 0,358"
    ),
    x = "Grado tumoral",
    y = "Actividad APOBEC (%)",
    caption = paste0(
      "Mediana: G2 = 47,0% · G3 = 57,6%\n",
      "Diferencia descriptiva; no significativa después de corrección FDR."
    )
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.caption = element_text(
      hjust = 0
    )
  )

ggsave(
  "results/tcga/figures/tcga_APOBEC_G2_vs_G3.png",
  p1,
  width = 7.5,
  height = 6,
  dpi = 300
)

# ==========================================================
# FIGURA 2
# SBS40 por histología
# ==========================================================

hist <- hc %>%
  filter(
    histology_group %in%
      c("Adenocarcinoma", "Squamous")
  ) %>%
  mutate(
    SBS40_pct = SBS40 * 100,
    histology_group = factor(
      histology_group,
      levels = c(
        "Adenocarcinoma",
        "Squamous"
      )
    )
  )

p2 <- ggplot(
  hist,
  aes(
    x = histology_group,
    y = SBS40_pct
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.45,
    size = 1.5
  ) +
  labs(
    title = "Contribución de SBS40 según subtipo histológico",
    subtitle = paste0(
      "Alta confianza (≥50 SNV; cosine ≥0,90)\n",
      "Wilcoxon p = 0,000764 · FDR = 0,0153"
    ),
    x = "Subtipo histológico",
    y = "Contribución SBS40 (%)",
    caption = paste0(
      "Mediana: Adenocarcinoma = 3,37% · Squamous = 0%\n",
      "Sensibilidad ≥100 SNV: Wilcoxon p = 0,000821."
    )
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.caption = element_text(
      hjust = 0
    )
  )

ggsave(
  "results/tcga/figures/tcga_SBS40_histology.png",
  p2,
  width = 8,
  height = 6,
  dpi = 300
)

# ==========================================================
# FIGURA 3
# Proporción de tumores SBS40 >= 1%
# ==========================================================

prev <- hist %>%
  group_by(histology_group) %>%
  summarise(
    n = n(),
    positive = sum(SBS40 >= 0.01),
    prevalence =
      100 * positive / n,
    .groups = "drop"
  )

p3 <- ggplot(
  prev,
  aes(
    x = histology_group,
    y = prevalence
  )
) +
  geom_col(
    width = 0.6
  ) +
  geom_text(
    aes(
      label = paste0(
        round(prevalence, 1),
        "%\n(",
        positive,
        "/",
        n,
        ")"
      )
    ),
    vjust = -0.4,
    size = 4
  ) +
  labs(
    title = "Tumores con contribución SBS40 ≥1%",
    subtitle = paste0(
      "Alta confianza (≥50 SNV; cosine ≥0,90)\n",
      "Fisher p = 0,00116"
    ),
    x = "Subtipo histológico",
    y = "Tumores con SBS40 ≥1% (%)"
  ) +
  scale_y_continuous(
    limits = c(
      0,
      max(prev$prevalence) + 12
    )
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank()
  )

ggsave(
  "results/tcga/figures/tcga_SBS40_prevalence_histology.png",
  p3,
  width = 7.5,
  height = 5.5,
  dpi = 300
)

# ==========================================================
# Resumen
# ==========================================================

cat("\n===== FIGURAS GENERADAS =====\n")

cat(
  "1. results/tcga/figures/tcga_APOBEC_G2_vs_G3.png\n"
)

cat(
  "2. results/tcga/figures/tcga_SBS40_histology.png\n"
)

cat(
  "3. results/tcga/figures/tcga_SBS40_prevalence_histology.png\n"
)

cat("\n===== PREVALENCIA SBS40 >=1% =====\n")
print(prev)
