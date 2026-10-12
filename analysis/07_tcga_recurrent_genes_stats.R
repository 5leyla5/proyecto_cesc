library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

dir.create("results/tcga/figures", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# TOP GENES GLOBAL
# ============================================================

overall <- read_csv(
  "results/tcga/recurrent_genes_overall.csv",
  show_col_types = FALSE
)

top15 <- overall %>%
  slice_max(mutated_patients, n = 15) %>%
  arrange(mutated_patients) %>%
  mutate(gene = factor(gene, levels = gene))

p_overall <- ggplot(
  top15,
  aes(x = gene, y = frequency_pct)
) +
  geom_col() +
  coord_flip() +
  geom_text(
    aes(label = paste0(round(frequency_pct, 1), "%")),
    hjust = -0.1,
    size = 3.5
  ) +
  labs(
    title = "Genes recurrentemente mutados en TCGA-CESC",
    subtitle = "Top 15 · mutaciones no sinónimas",
    x = "Gen",
    y = "Pacientes con mutación (%)"
  ) +
  ylim(0, max(top15$frequency_pct) * 1.15) +
  theme_minimal(base_size = 13)

ggsave(
  "results/tcga/figures/tcga_top15_recurrent_genes.png",
  p_overall,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================================
# G2 vs G3
# ============================================================

grade <- read_csv(
  "results/tcga/recurrent_genes_by_grade.csv",
  show_col_types = FALSE
)

grade_stats <- grade %>%
  rowwise() %>%
  mutate(
    p_value = fisher.test(
      matrix(
        c(
          G2_mutated,
          G2_total - G2_mutated,
          G3_mutated,
          G3_total - G3_mutated
        ),
        nrow = 2,
        byrow = TRUE
      )
    )$p.value
  ) %>%
  ungroup() %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH"),
    difference_pct = G3_frequency_pct - G2_frequency_pct
  ) %>%
  arrange(p_value)

write_csv(
  grade_stats,
  "results/tcga/recurrent_genes_grade_fisher.csv"
)

top_grade_genes <- grade %>%
  mutate(total_mutated = G2_mutated + G3_mutated) %>%
  slice_max(total_mutated, n = 15) %>%
  pull(gene)

grade_plot <- grade %>%
  filter(gene %in% top_grade_genes) %>%
  select(
    gene,
    G2_frequency_pct,
    G3_frequency_pct
  ) %>%
  pivot_longer(
    cols = c(G2_frequency_pct, G3_frequency_pct),
    names_to = "grade",
    values_to = "frequency"
  ) %>%
  mutate(
    grade = recode(
      grade,
      G2_frequency_pct = "G2",
      G3_frequency_pct = "G3"
    ),
    gene = factor(gene, levels = rev(top_grade_genes))
  )

p_grade <- ggplot(
  grade_plot,
  aes(x = gene, y = frequency, fill = grade)
) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(
    title = "Genes recurrentes según grado tumoral",
    subtitle = "TCGA-CESC · G2 vs G3",
    x = "Gen",
    y = "Pacientes con mutación (%)",
    fill = "Grado"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  "results/tcga/figures/tcga_recurrent_genes_G2_vs_G3.png",
  p_grade,
  width = 9,
  height = 7,
  dpi = 300
)

# ============================================================
# SQUAMOUS vs ADENOCARCINOMA
# ============================================================

hist <- read_csv(
  "results/tcga/recurrent_genes_by_histology.csv",
  show_col_types = FALSE
)

hist_stats <- hist %>%
  rowwise() %>%
  mutate(
    p_value = fisher.test(
      matrix(
        c(
          Squamous_mutated,
          Squamous_total - Squamous_mutated,
          Adenocarcinoma_mutated,
          Adenocarcinoma_total - Adenocarcinoma_mutated
        ),
        nrow = 2,
        byrow = TRUE
      )
    )$p.value
  ) %>%
  ungroup() %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH"),
    difference_pct =
      Adenocarcinoma_frequency_pct -
      Squamous_frequency_pct
  ) %>%
  arrange(p_value)

write_csv(
  hist_stats,
  "results/tcga/recurrent_genes_histology_fisher.csv"
)

top_hist_genes <- hist %>%
  mutate(total_mutated = Squamous_mutated + Adenocarcinoma_mutated) %>%
  slice_max(total_mutated, n = 15) %>%
  pull(gene)

hist_plot <- hist %>%
  filter(gene %in% top_hist_genes) %>%
  select(
    gene,
    Squamous_frequency_pct,
    Adenocarcinoma_frequency_pct
  ) %>%
  pivot_longer(
    cols = c(
      Squamous_frequency_pct,
      Adenocarcinoma_frequency_pct
    ),
    names_to = "histology",
    values_to = "frequency"
  ) %>%
  mutate(
    histology = recode(
      histology,
      Squamous_frequency_pct = "Squamous",
      Adenocarcinoma_frequency_pct = "Adenocarcinoma"
    ),
    gene = factor(gene, levels = rev(top_hist_genes))
  )

p_hist <- ggplot(
  hist_plot,
  aes(x = gene, y = frequency, fill = histology)
) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(
    title = "Genes recurrentes según subtipo histológico",
    subtitle = "TCGA-CESC · escamoso vs adenocarcinoma",
    x = "Gen",
    y = "Pacientes con mutación (%)",
    fill = "Histología"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  "results/tcga/figures/tcga_recurrent_genes_histology.png",
  p_hist,
  width = 9,
  height = 7,
  dpi = 300
)

cat("\n===== G2 vs G3: genes con menor p-value =====\n")
print(
  grade_stats %>%
    select(
      gene,
      G2_frequency_pct,
      G3_frequency_pct,
      p_value,
      FDR
    ) %>%
    slice_head(n = 15)
)

cat("\n===== SQUAMOUS vs ADENOCARCINOMA: genes con menor p-value =====\n")
print(
  hist_stats %>%
    select(
      gene,
      Squamous_frequency_pct,
      Adenocarcinoma_frequency_pct,
      p_value,
      FDR
    ) %>%
    slice_head(n = 15)
)

cat("\nFiguras y tablas creadas correctamente.\n")
