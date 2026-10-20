library(readr)
library(dplyr)

# ==========================================================
# 1. Leer datos
# ==========================================================

rel <- read_csv(
  "results/tcga/COSMIC_SBS_relative.csv",
  show_col_types = FALSE
)

qc <- read_csv(
  "results/tcga/COSMIC_SBS_fit_quality.csv",
  show_col_types = FALSE
)

clinical <- read_csv(
  "samplesheet/tcga_cesc_samplesheet.csv",
  show_col_types = FALSE
)

# ==========================================================
# 2. Revisar nombres de firmas
# ==========================================================

cat("\n===== COLUMNAS DE FIRMAS =====\n")

signature_cols <- setdiff(
  names(rel),
  "case_id"
)

print(signature_cols)

# ==========================================================
# 3. Integrar metadata
# ==========================================================

dat <- rel %>%
  left_join(
    qc,
    by = "case_id"
  ) %>%
  left_join(
    clinical %>%
      select(
        case_id,
        tumor_grade,
        histology_group,
        grade_analysis,
        histology_analysis
      ),
    by = "case_id"
  )

# ==========================================================
# 4. Actividad APOBEC combinada
# ==========================================================

if (all(c("SBS2", "SBS13") %in% names(dat))) {

  dat <- dat %>%
    mutate(
      APOBEC = SBS2 + SBS13
    )

} else {

  stop(
    "No se encontraron columnas SBS2 y SBS13."
  )
}

# ==========================================================
# 5. Definir subconjunto de alta confianza
# ==========================================================

dat <- dat %>%
  mutate(
    high_confidence =
      total_snvs >= 50 &
      cosine >= 0.90
  )

# ==========================================================
# 6. Guardar tabla maestra
# ==========================================================

write_csv(
  dat,
  "results/tcga/COSMIC_SBS_with_metadata.csv"
)

# ==========================================================
# 7. Resumen QC
# ==========================================================

cat("\n===== COHORTE =====\n")
cat("Total:", nrow(dat), "\n")

cat(
  ">=50 SNV:",
  sum(dat$total_snvs >= 50, na.rm = TRUE),
  "\n"
)

cat(
  "Cosine >=0.90:",
  sum(dat$cosine >= 0.90, na.rm = TRUE),
  "\n"
)

cat(
  "Alta confianza (>=50 SNV + cosine >=0.90):",
  sum(dat$high_confidence, na.rm = TRUE),
  "\n"
)

# ==========================================================
# 8. APOBEC global
# ==========================================================

cat("\n===== ACTIVIDAD APOBEC (SBS2 + SBS13) =====\n")

cat(
  "Mediana cohorte completa:",
  round(
    median(dat$APOBEC, na.rm = TRUE) * 100,
    2
  ),
  "%\n"
)

cat(
  "Media cohorte completa:",
  round(
    mean(dat$APOBEC, na.rm = TRUE) * 100,
    2
  ),
  "%\n"
)

hc <- dat %>%
  filter(high_confidence)

cat(
  "Mediana alta confianza:",
  round(
    median(hc$APOBEC, na.rm = TRUE) * 100,
    2
  ),
  "%\n"
)

cat(
  "Media alta confianza:",
  round(
    mean(hc$APOBEC, na.rm = TRUE) * 100,
    2
  ),
  "%\n"
)

# ==========================================================
# 9. Por grado
# ==========================================================

cat("\n===== APOBEC POR GRADO =====\n")

print(
  hc %>%
    filter(
      tumor_grade %in% c("G2", "G3")
    ) %>%
    group_by(tumor_grade) %>%
    summarise(
      n = n(),
      median_APOBEC = median(APOBEC) * 100,
      mean_APOBEC = mean(APOBEC) * 100,
      Q1 = quantile(APOBEC, 0.25) * 100,
      Q3 = quantile(APOBEC, 0.75) * 100,
      .groups = "drop"
    )
)

# ==========================================================
# 10. Por histología
# ==========================================================

cat("\n===== APOBEC POR HISTOLOGÍA =====\n")

print(
  hc %>%
    filter(
      histology_group %in%
        c("Squamous", "Adenocarcinoma")
    ) %>%
    group_by(histology_group) %>%
    summarise(
      n = n(),
      median_APOBEC = median(APOBEC) * 100,
      mean_APOBEC = mean(APOBEC) * 100,
      Q1 = quantile(APOBEC, 0.25) * 100,
      Q3 = quantile(APOBEC, 0.75) * 100,
      .groups = "drop"
    )
)

cat("\nArchivo creado:\n")
cat(
  "results/tcga/COSMIC_SBS_with_metadata.csv\n"
)
