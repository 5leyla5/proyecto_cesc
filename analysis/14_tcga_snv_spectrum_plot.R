library(readr)
library(dplyr)
library(ggplot2)

dir.create(
  "results/tcga/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

df <- read_csv(
  "results/tcga/snv_spectrum_overall.csv",
  show_col_types = FALSE
)

df <- df %>%
  mutate(
    substitution_class = factor(
      substitution_class,
      levels = c(
        "C>A",
        "C>G",
        "C>T",
        "T>A",
        "T>C",
        "T>G"
      )
    )
  )

p <- ggplot(
  df,
  aes(
    x = substitution_class,
    y = frequency_pct
  )
) +
  geom_col(width = 0.75) +
  geom_text(
    aes(
      label = paste0(
        format(count, big.mark = "."),
        "\n(",
        frequency_pct,
        "%)"
      )
    ),
    vjust = -0.3,
    size = 3.5
  ) +
  labs(
    title = "Espectro de sustituciones somáticas en TCGA-CESC",
    subtitle = "66.646 SNV · representación SBS6",
    x = "Clase de sustitución",
    y = "Frecuencia (%)"
  ) +
  scale_y_continuous(
    limits = c(0, 68),
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme_minimal(base_size = 13)

ggsave(
  "results/tcga/figures/tcga_SBS6_spectrum.png",
  p,
  width = 8,
  height = 5.5,
  dpi = 300
)

cat("\n===== SBS6 =====\n")
print(df)

cat("\nFigura creada:\n")
cat("results/tcga/figures/tcga_SBS6_spectrum.png\n")
