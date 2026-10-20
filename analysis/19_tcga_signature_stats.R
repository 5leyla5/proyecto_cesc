library(readr)
library(dplyr)

# ==========================================================
# 1. Leer tabla integrada
# ==========================================================

dat <- read_csv(
  "results/tcga/COSMIC_SBS_with_metadata.csv",
  show_col_types = FALSE
)

# Solo ajustes de alta confianza
hc <- dat %>%
  filter(high_confidence)

cat("\n===== COHORTE ALTA CONFIANZA =====\n")
cat("Total:", nrow(hc), "\n")

# ==========================================================
# 2. Firmas a evaluar
# ==========================================================

signatures <- c(
  "SBS1",
  "SBS2",
  "SBS3",
  "SBS5",
  "SBS6",
  "SBS10a",
  "SBS10b",
  "SBS13",
  "SBS15",
  "SBS17a",
  "SBS17b",
  "SBS26",
  "SBS28",
  "SBS33",
  "SBS40",
  "SBS44",
  "SBS45",
  "SBS50",
  "SBS59",
  "APOBEC"
)

# ==========================================================
# 3. Función para comparar dos grupos
# ==========================================================

compare_signatures <- function(
  data,
  group_var,
  group1,
  group2
) {

  results <- list()

  for (sig in signatures) {

    x <- data[
      data[[group_var]] == group1,
      sig,
      drop = TRUE
    ]

    y <- data[
      data[[group_var]] == group2,
      sig,
      drop = TRUE
    ]

    x <- x[!is.na(x)]
    y <- y[!is.na(y)]

    test <- wilcox.test(
      x,
      y,
      exact = FALSE
    )

    results[[sig]] <- data.frame(
      signature = sig,

      group1 = group1,
      n_group1 = length(x),
      median_group1 = median(x) * 100,
      mean_group1 = mean(x) * 100,

      group2 = group2,
      n_group2 = length(y),
      median_group2 = median(y) * 100,
      mean_group2 = mean(y) * 100,

      median_difference =
        (median(y) - median(x)) * 100,

      p_value = test$p.value
    )
  }

  out <- bind_rows(results)

  out <- out %>%
    mutate(
      FDR = p.adjust(
        p_value,
        method = "BH"
      )
    ) %>%
    arrange(FDR, p_value)

  return(out)
}

# ==========================================================
# 4. G2 vs G3
# ==========================================================

grade_data <- hc %>%
  filter(
    tumor_grade %in% c("G2", "G3")
  )

grade_results <- compare_signatures(
  grade_data,
  "tumor_grade",
  "G2",
  "G3"
)

write_csv(
  grade_results,
  "results/tcga/COSMIC_SBS_G2_vs_G3_stats.csv"
)

# ==========================================================
# 5. Squamous vs Adenocarcinoma
# ==========================================================

histology_data <- hc %>%
  filter(
    histology_group %in%
      c("Squamous", "Adenocarcinoma")
  )

histology_results <- compare_signatures(
  histology_data,
  "histology_group",
  "Adenocarcinoma",
  "Squamous"
)

write_csv(
  histology_results,
  "results/tcga/COSMIC_SBS_histology_stats.csv"
)

# ==========================================================
# 6. Mostrar resultados
# ==========================================================

cat("\n===== G2 vs G3 =====\n")
cat(
  "G2:",
  sum(grade_data$tumor_grade == "G2"),
  "\n"
)
cat(
  "G3:",
  sum(grade_data$tumor_grade == "G3"),
  "\n\n"
)

print(
  grade_results %>%
    select(
      signature,
      median_group1,
      median_group2,
      median_difference,
      p_value,
      FDR
    )
)

cat("\n===== SQUAMOUS vs ADENOCARCINOMA =====\n")
cat(
  "Adenocarcinoma:",
  sum(
    histology_data$histology_group ==
      "Adenocarcinoma"
  ),
  "\n"
)

cat(
  "Squamous:",
  sum(
    histology_data$histology_group ==
      "Squamous"
  ),
  "\n\n"
)

print(
  histology_results %>%
    select(
      signature,
      median_group1,
      median_group2,
      median_difference,
      p_value,
      FDR
    )
)

# ==========================================================
# 7. Resumen APOBEC específico
# ==========================================================

cat("\n===== APOBEC =====\n")

cat("\nG2 vs G3:\n")

print(
  grade_results %>%
    filter(signature == "APOBEC") %>%
    select(
      median_group1,
      median_group2,
      median_difference,
      p_value,
      FDR
    )
)

cat("\nHistología:\n")

print(
  histology_results %>%
    filter(signature == "APOBEC") %>%
    select(
      median_group1,
      median_group2,
      median_difference,
      p_value,
      FDR
    )
)

cat("\nArchivos creados:\n")
cat(
  "results/tcga/COSMIC_SBS_G2_vs_G3_stats.csv\n"
)
cat(
  "results/tcga/COSMIC_SBS_histology_stats.csv\n"
)
