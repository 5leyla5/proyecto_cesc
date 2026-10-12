library(readr)
library(dplyr)
library(ggplot2)

df <- read_csv(
  "results/tcga/tcga_oncoplot_top15.csv",
  show_col_types = FALSE
)

dir.create(
  "results/tcga/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

# Orden de genes por frecuencia
gene_order <- df %>%
  filter(variant_class != "No_mutation") %>%
  distinct(case_id, gene) %>%
  count(gene, sort = TRUE) %>%
  pull(gene)

# Ordenar pacientes por histología y número de genes mutados
sample_order <- df %>%
  filter(variant_class != "No_mutation") %>%
  distinct(case_id, histology_group, gene) %>%
  count(case_id, histology_group, name = "n_mutated_genes") %>%
  arrange(
    histology_group,
    desc(n_mutated_genes),
    case_id
  ) %>%
  pull(case_id)

df <- df %>%
  mutate(
    gene = factor(gene, levels = rev(gene_order)),
    case_id = factor(case_id, levels = sample_order),
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
    . ~ histology_group,
    scales = "free_x",
    space = "free_x"
  ) +
  labs(
    title = "Paisaje mutacional de TCGA-CESC",
    subtitle = "Top 15 genes recurrentes · carcinoma escamoso vs adenocarcinoma",
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
  "results/tcga/figures/tcga_oncoplot_top15_histology.png",
  p,
  width = 15,
  height = 7,
  dpi = 300
)

cat("\nOncoplot creado:\n")
cat(
  "results/tcga/figures/tcga_oncoplot_top15_histology.png\n"
)
