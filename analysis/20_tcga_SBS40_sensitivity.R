library(readr)
library(dplyr)

# ==========================================================
# 1. Leer datos
# ==========================================================

dat <- read_csv(
  "results/tcga/COSMIC_SBS_with_metadata.csv",
  show_col_types = FALSE
)

# Solo histologías principales
dat_hist <- dat %>%
  filter(
    histology_group %in%
      c("Adenocarcinoma", "Squamous")
  )

# ==========================================================
# 2. Cohorte original de alta confianza
#    >=50 SNV + cosine >=0.90
# ==========================================================

hc50 <- dat_hist %>%
  filter(
    total_snvs >= 50,
    cosine >= 0.90
  )

# ==========================================================
# 3. Cohorte estricta
#    >=100 SNV + cosine >=0.90
# ==========================================================

hc100 <- dat_hist %>%
  filter(
    total_snvs >= 100,
    cosine >= 0.90
  )

# ==========================================================
# 4. Función resumen SBS40
# ==========================================================

summarise_sbs40 <- function(df, label) {

  cat("\n====================================\n")
  cat(label, "\n")
  cat("====================================\n")

  summary_tbl <- df %>%
    group_by(histology_group) %>%
    summarise(
      n = n(),

      median_SBS40 =
        median(SBS40, na.rm = TRUE) * 100,

      mean_SBS40 =
        mean(SBS40, na.rm = TRUE) * 100,

      Q1 =
        quantile(SBS40, 0.25, na.rm = TRUE) * 100,

      Q3 =
        quantile(SBS40, 0.75, na.rm = TRUE) * 100,

      SBS40_positive_n =
        sum(SBS40 > 0, na.rm = TRUE),

      SBS40_positive_pct =
        mean(SBS40 > 0, na.rm = TRUE) * 100,

      SBS40_ge1pct_n =
        sum(SBS40 >= 0.01, na.rm = TRUE),

      SBS40_ge1pct_pct =
        mean(SBS40 >= 0.01, na.rm = TRUE) * 100,

      .groups = "drop"
    )

  print(summary_tbl)

  # --------------------------------------------------------
  # Wilcoxon: exposición continua
  # --------------------------------------------------------

  adeno <- df %>%
    filter(histology_group == "Adenocarcinoma") %>%
    pull(SBS40)

  squamous <- df %>%
    filter(histology_group == "Squamous") %>%
    pull(SBS40)

  w <- wilcox.test(
    adeno,
    squamous,
    exact = FALSE
  )

  cat("\nWilcoxon SBS40:\n")
  cat("p-value =", format(w$p.value, digits = 6), "\n")

  # --------------------------------------------------------
  # Fisher: SBS40 > 0
  # --------------------------------------------------------

  tab_positive <- table(
    df$histology_group,
    df$SBS40 > 0
  )

  cat("\nTabla SBS40 > 0:\n")
  print(tab_positive)

  fisher_positive <- fisher.test(tab_positive)

  cat(
    "Fisher SBS40 > 0: p =",
    format(fisher_positive$p.value, digits = 6),
    "\n"
  )

  # --------------------------------------------------------
  # Fisher: SBS40 >= 1%
  # --------------------------------------------------------

  tab_1pct <- table(
    df$histology_group,
    df$SBS40 >= 0.01
  )

  cat("\nTabla SBS40 >= 1%:\n")
  print(tab_1pct)

  fisher_1pct <- fisher.test(tab_1pct)

  cat(
    "Fisher SBS40 >= 1%: p =",
    format(fisher_1pct$p.value, digits = 6),
    "\n"
  )

  return(
    list(
      summary = summary_tbl,
      wilcoxon_p = w$p.value,
      fisher_positive_p =
        fisher_positive$p.value,
      fisher_1pct_p =
        fisher_1pct$p.value
    )
  )
}

# ==========================================================
# 5. Ejecutar análisis
# ==========================================================

res50 <- summarise_sbs40(
  hc50,
  "ALTA CONFIANZA: >=50 SNV + cosine >=0.90"
)

res100 <- summarise_sbs40(
  hc100,
  "SENSIBILIDAD: >=100 SNV + cosine >=0.90"
)

# ==========================================================
# 6. Guardar tablas
# ==========================================================

write_csv(
  res50$summary,
  "results/tcga/SBS40_histology_summary_ge50.csv"
)

write_csv(
  res100$summary,
  "results/tcga/SBS40_histology_summary_ge100.csv"
)

# Tabla de tests
tests <- data.frame(
  analysis = c(
    "ge50_cosine90",
    "ge100_cosine90"
  ),

  wilcoxon_p = c(
    res50$wilcoxon_p,
    res100$wilcoxon_p
  ),

  fisher_SBS40_gt0_p = c(
    res50$fisher_positive_p,
    res100$fisher_positive_p
  ),

  fisher_SBS40_ge1pct_p = c(
    res50$fisher_1pct_p,
    res100$fisher_1pct_p
  )
)

write_csv(
  tests,
  "results/tcga/SBS40_histology_sensitivity_tests.csv"
)

cat("\n===== RESUMEN TESTS =====\n")
print(tests)

cat("\nArchivos creados:\n")
cat(
  "results/tcga/SBS40_histology_summary_ge50.csv\n"
)
cat(
  "results/tcga/SBS40_histology_summary_ge100.csv\n"
)
cat(
  "results/tcga/SBS40_histology_sensitivity_tests.csv\n"
)
