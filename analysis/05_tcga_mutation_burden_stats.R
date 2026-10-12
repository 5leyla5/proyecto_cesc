library(readr)
library(dplyr)
library(ggplot2)

dir.create("results/tcga/figures", recursive = TRUE, showWarnings = FALSE)

df <- read_csv(
  "results/tcga/tcga_cesc_mutation_burden.csv",
  show_col_types = FALSE
)

# ============================================================
# G2 vs G3
# ============================================================

grade_df <- df %>%
  filter(tumor_grade %in% c("G2", "G3")) %>%
  mutate(tumor_grade = factor(tumor_grade, levels = c("G2", "G3")))

grade_summary <- grade_df %>%
  group_by(tumor_grade) %>%
  summarise(
    n = n(),
    median = median(mutation_burden),
    Q1 = quantile(mutation_burden, 0.25),
    Q3 = quantile(mutation_burden, 0.75),
    mean = mean(mutation_burden),
    min = min(mutation_burden),
    max = max(mutation_burden),
    .groups = "drop"
  )

grade_test <- wilcox.test(
  mutation_burden ~ tumor_grade,
  data = grade_df,
  exact = FALSE
)

cat("\n===== G2 vs G3 =====\n")
print(grade_summary)
cat("\nWilcoxon p-value:", grade_test$p.value, "\n")

write_csv(
  grade_summary,
  "results/tcga/mutation_burden_grade_summary.csv"
)

p_grade <- ggplot(
  grade_df,
  aes(x = tumor_grade, y = mutation_burden)
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.55
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.45,
    size = 1.5
  ) +
  scale_y_log10() +
  labs(
    title = "Carga mutacional según grado tumoral",
    subtitle = "TCGA-CESC: G2 vs G3",
    x = "Grado tumoral",
    y = "Mutaciones no sinónimas por tumor (escala log10)"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  "results/tcga/figures/tcga_mutation_burden_G2_vs_G3.png",
  p_grade,
  width = 7,
  height = 6,
  dpi = 300
)

# ============================================================
# Squamous vs Adenocarcinoma
# ============================================================

hist_df <- df %>%
  filter(histology_group %in% c("Squamous", "Adenocarcinoma")) %>%
  mutate(
    histology_group = factor(
      histology_group,
      levels = c("Squamous", "Adenocarcinoma")
    )
  )

hist_summary <- hist_df %>%
  group_by(histology_group) %>%
  summarise(
    n = n(),
    median = median(mutation_burden),
    Q1 = quantile(mutation_burden, 0.25),
    Q3 = quantile(mutation_burden, 0.75),
    mean = mean(mutation_burden),
    min = min(mutation_burden),
    max = max(mutation_burden),
    .groups = "drop"
  )

hist_test <- wilcox.test(
  mutation_burden ~ histology_group,
  data = hist_df,
  exact = FALSE
)

cat("\n===== SQUAMOUS vs ADENOCARCINOMA =====\n")
print(hist_summary)
cat("\nWilcoxon p-value:", hist_test$p.value, "\n")

write_csv(
  hist_summary,
  "results/tcga/mutation_burden_histology_summary.csv"
)

p_hist <- ggplot(
  hist_df,
  aes(x = histology_group, y = mutation_burden)
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.55
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.45,
    size = 1.5
  ) +
  scale_y_log10() +
  labs(
    title = "Carga mutacional según subtipo histológico",
    subtitle = "TCGA-CESC: carcinoma escamoso vs adenocarcinoma",
    x = "Subtipo histológico",
    y = "Mutaciones no sinónimas por tumor (escala log10)"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  "results/tcga/figures/tcga_mutation_burden_histology.png",
  p_hist,
  width = 7,
  height = 6,
  dpi = 300
)

cat("\nFiguras creadas en results/tcga/figures/\n")
