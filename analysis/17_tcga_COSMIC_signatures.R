library(sigminer)

# ==========================================================
# Archivos
# ==========================================================

SBS96_FILE <- "results/tcga/SBS96_by_patient.csv"
COUNTS_FILE <- "results/tcga/titv_by_patient.csv"

OUT_ABS <- "results/tcga/COSMIC_SBS_absolute.csv"
OUT_REL <- "results/tcga/COSMIC_SBS_relative.csv"
OUT_QC  <- "results/tcga/COSMIC_SBS_fit_quality.csv"

# ==========================================================
# 1. Leer datos
# ==========================================================

x <- read.csv(
  SBS96_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

counts <- read.csv(
  COUNTS_FILE,
  stringsAsFactors = FALSE
)

cat("\n===== DATOS =====\n")
cat("Pacientes SBS96:", nrow(x), "\n")
cat("Pacientes conteos:", nrow(counts), "\n")

# ==========================================================
# 2. Construir matriz 96 contextos × pacientes
# ==========================================================

context_cols <- setdiff(
  colnames(x),
  c("case_id", "tumor_sample")
)

mat <- as.matrix(
  x[, context_cols]
)

storage.mode(mat) <- "numeric"

rownames(mat) <- x$case_id

mat <- t(mat)

cat("Dimensiones matriz:", dim(mat)[1], "x", dim(mat)[2], "\n")

if (nrow(mat) != 96) {
  stop("ERROR: la matriz no tiene 96 contextos.")
}

# ==========================================================
# 3. Firmas SBS descritas para Cervix-CA / TCGA / WES
# ==========================================================

cervix_signatures <- c(
  "1",
  "2",
  "3",
  "5",
  "6",
  "10a",
  "10b",
  "13",
  "15",
  "17a",
  "17b",
  "26",
  "28",
  "33",
  "40",
  "44",
  "45",
  "50",
  "59"
)

cat("\n===== FIRMAS CANDIDATAS =====\n")
print(cervix_signatures)

# ==========================================================
# 4. Ajuste NNLS
# ==========================================================

cat("\n===== AJUSTANDO FIRMAS =====\n")

fit <- sig_fit(
  catalogue_matrix = mat,
  sig_index = cervix_signatures,
  sig_db = "SBS",
  db_type = "human-exome",
  method = "NNLS",
  type = "absolute",
  return_error = TRUE,
  show_index = FALSE
)

expo_abs <- fit$expo

cat("\nDimensiones exposiciones:")
cat(nrow(expo_abs), "firmas x", ncol(expo_abs), "pacientes\n")

# ==========================================================
# 5. Exposición relativa
# ==========================================================

expo_rel <- sweep(
  expo_abs,
  2,
  colSums(expo_abs),
  "/"
)

expo_rel[!is.finite(expo_rel)] <- 0

# ==========================================================
# 6. Preparar tablas
# ==========================================================

# Transponer:
# una fila por paciente
abs_df <- as.data.frame(
  t(expo_abs),
  check.names = FALSE
)

abs_df$case_id <- rownames(abs_df)

abs_df <- abs_df[
  c("case_id", setdiff(colnames(abs_df), "case_id"))
]

rel_df <- as.data.frame(
  t(expo_rel),
  check.names = FALSE
)

rel_df$case_id <- rownames(rel_df)

rel_df <- rel_df[
  c("case_id", setdiff(colnames(rel_df), "case_id"))
]

# ==========================================================
# 7. Firma dominante por paciente
# ==========================================================

dominant_index <- apply(
  expo_rel,
  2,
  which.max
)

dominant_signature <- rownames(expo_rel)[
  dominant_index
]

dominant_fraction <- apply(
  expo_rel,
  2,
  max
)

# ==========================================================
# 8. Tabla de QC
# ==========================================================

qc <- data.frame(
  case_id = colnames(expo_abs),
  cosine = as.numeric(
    fit$cosine[colnames(expo_abs)]
  ),
  error = as.numeric(
    fit$errors[colnames(expo_abs)]
  ),
  reconstructed_mutations = as.numeric(
    colSums(expo_abs)
  ),
  dominant_signature = dominant_signature,
  dominant_fraction = dominant_fraction,
  stringsAsFactors = FALSE
)

# Añadir SNV observadas
idx <- match(
  qc$case_id,
  counts$case_id
)

qc$total_snvs <- counts$total_snvs[idx]

qc$reconstructed_fraction <- (
  qc$reconstructed_mutations /
  qc$total_snvs
)

# Flag para análisis de sensibilidad
qc$snv_ge_50 <- qc$total_snvs >= 50

# Orden
qc <- qc[
  c(
    "case_id",
    "total_snvs",
    "cosine",
    "error",
    "reconstructed_mutations",
    "reconstructed_fraction",
    "dominant_signature",
    "dominant_fraction",
    "snv_ge_50"
  )
]

# ==========================================================
# 9. Guardar resultados
# ==========================================================

write.csv(
  abs_df,
  OUT_ABS,
  row.names = FALSE
)

write.csv(
  rel_df,
  OUT_REL,
  row.names = FALSE
)

write.csv(
  qc,
  OUT_QC,
  row.names = FALSE
)

# ==========================================================
# 10. Resumen
# ==========================================================

cat("\n===== CALIDAD DEL AJUSTE =====\n")

cat(
  "Cosine mínimo:",
  round(min(qc$cosine, na.rm = TRUE), 4),
  "\n"
)

cat(
  "Cosine mediano:",
  round(median(qc$cosine, na.rm = TRUE), 4),
  "\n"
)

cat(
  "Cosine medio:",
  round(mean(qc$cosine, na.rm = TRUE), 4),
  "\n"
)

cat(
  "Cosine máximo:",
  round(max(qc$cosine, na.rm = TRUE), 4),
  "\n"
)

cat("\nCosine >= 0.90:",
    sum(qc$cosine >= 0.90, na.rm = TRUE),
    "/", nrow(qc), "\n")

cat("Cosine >= 0.95:",
    sum(qc$cosine >= 0.95, na.rm = TRUE),
    "/", nrow(qc), "\n")

cat("Cosine < 0.90:",
    sum(qc$cosine < 0.90, na.rm = TRUE),
    "\n")

cat("\n===== FIRMA DOMINANTE =====\n")

dominant_table <- sort(
  table(qc$dominant_signature),
  decreasing = TRUE
)

print(dominant_table)

cat("\n===== FIRMA DOMINANTE SOLO >=50 SNV =====\n")

dominant_50 <- sort(
  table(
    qc$dominant_signature[
      qc$snv_ge_50
    ]
  ),
  decreasing = TRUE
)

print(dominant_50)

cat("\n===== 15 MEJORES AJUSTES =====\n")

print(
  head(
    qc[
      order(-qc$cosine),
      c(
        "case_id",
        "total_snvs",
        "cosine",
        "dominant_signature",
        "dominant_fraction"
      )
    ],
    15
  )
)

cat("\n===== 15 AJUSTES MÁS BAJOS =====\n")

print(
  head(
    qc[
      order(qc$cosine),
      c(
        "case_id",
        "total_snvs",
        "cosine",
        "dominant_signature",
        "dominant_fraction"
      )
    ],
    15
  )
)

cat("\nArchivos creados:\n")
cat(OUT_ABS, "\n")
cat(OUT_REL, "\n")
cat(OUT_QC, "\n")
