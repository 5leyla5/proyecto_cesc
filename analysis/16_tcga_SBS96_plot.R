library(readr)
library(dplyr)
library(ggplot2)

dir.create(
  "results/tcga/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

df <- read_csv(
  "results/tcga/SBS96_overall.csv",
  show_col_types = FALSE
)

# Orden estándar por clase SBS6
sbs_order <- c(
  "C>A",
  "C>G",
  "C>T",
  "T>A",
  "T>C",
  "T>G"
)

df <- df %>%
  mutate(
    substitution_class = factor(
      substitution_class,
      levels = sbs_order
    ),
    context = factor(
      context,
      levels = context
    )
  )

p <- ggplot(
  df,
  aes(
    x = context,
    y = frequency_pct,
    fill = substitution_class
  )
) +
  geom_col(width = 0.9) +
  facet_grid(
    . ~ substitution_class,
    scales = "free_x",
    space = "free_x"
  ) +
  labs(
    title = "Espectro mutacional SBS96 en TCGA-CESC",
    subtitle = "66.646 sustituciones somáticas · 286 pacientes",
    x = "Contexto trinucleotídico",
    y = "Frecuencia (%)",
    fill = "Sustitución"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 6
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(
      face = "bold"
    ),
    legend.position = "none"
  )

ggsave(
  "results/tcga/figures/tcga_SBS96_spectrum.png",
  p,
  width = 15,
  height = 6,
  dpi = 300
)

cat("\n===== TOP 20 CONTEXTOS SBS96 =====\n")

print(
  df %>%
    arrange(desc(count)) %>%
    select(
      context,
      substitution_class,
      count,
      frequency_pct
    ) %>%
    slice_head(n = 20)
)

cat("\nFigura creada:\n")
cat(
  "results/tcga/figures/tcga_SBS96_spectrum.png\n"
)
