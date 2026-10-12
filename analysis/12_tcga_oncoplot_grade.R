library(readr)
library(dplyr)
library(ggplot2)

df <- read_csv(
  "results/tcga/tcga_oncoplot_top15_grade.csv",
  show_col_types = FALSE
)

dir.create(
  "results/tcga/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\n===== COHORTE =====\n")
cat("Total:", n_distinct(df$case_id), "\n")
cat(
  "G2:",
  n_distinct(df$case_id[df$tumor_grade == "G2"]),
  "\n"
)
cat(
  "G3:",
  n_distinct(df$case_id[df$tumor_grade == "G3"]),
  "\n"
)

# Orden de genes según frecuencia dentro de G2/G3
gene_order <- df %>%
  filter(variant_class != "No_mutation") %>%
  distinct(case_id, gene) %>%
  count(gene, sort = TRUE) %>%
  pull(gene)

# Incluye TODOS los pacientes, incluso si no tienen
# mutaciones en ninguno de los top 15
sample_order <- df %>%
  mutate(
    mutated = variant_class != "No_mutation"
  ) %>%
  group_by(case_id, tumor_grade) %>%
  summarise(
    n_mutated_genes = sum(mutated),
    .groups = "drop"
  ) %>%
  arrange(
    tumor_grade,
    desc(n_mutated_genes),
    case_id
  ) %>%
  pull(case_id)

df <- df %>%
  mutate(
    gene = factor(
      gene,
      levels = rev(gene_order)
    ),
    case_id = factor(
      case_id,
      levels = sample_order
    ),
    variant_class = recode(
      variant_class,
      "Frame_Shift_Del" = "Frameshift",
      "Frame_Shift_Ins" = "Frameshift",
      "Nonsense_Mutation" = "Nonsense",
      "Nonstop_Mutation" = "Nonsense",
      "Translation_Start_Site" = "Start/stop",
      "Splice_Site" = "Splice site",
      "In_Frame_Del" = "In-frame",
      "In_Frame_Ins" = "In-frame",
      "Missense_Mutation" = "Missense",
      "No_mutation" = "No mutation"
    )
  )

# Comprobar que no perdimos pacientes
cat("\n===== DESPUÉS DEL ORDENAMIENTO =====\n")
cat(
  "Pacientes sin NA:",
  n_distinct(df$case_id[!is.na(df$case_id)]),
  "\n"
)

p <- ggplot(
  df,
  aes(
    x = case_id,
    y = gene,
    fill = variant_class
  )
) +
  geom_tile() +
  facet_grid(
    . ~ tumor_grade,
    scales = "free_x",
    space = "free_x"
  ) +
  labs(
    title = "Paisaje mutacional de TCGA-CESC según grado tumoral",
    subtitle = "Top 15 genes recurrentes · G2 (n=126) vs G3 (n=111)",
    x = "Pacientes",
    y = "Gen",
    fill = "Tipo de variante"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(
  "results/tcga/figures/tcga_oncoplot_top15_G2_vs_G3.png",
  p,
  width = 15,
  height = 7,
  dpi = 300
)

cat("\nOncoplot corregido creado:\n")
cat(
  "results/tcga/figures/tcga_oncoplot_top15_G2_vs_G3.png\n"
)
