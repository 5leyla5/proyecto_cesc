library(sigminer)

# ==========================================================
# 1. Leer matriz SBS96
# ==========================================================

x <- read.csv(
  "results/tcga/SBS96_by_patient.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

counts <- read.csv(
  "results/tcga/titv_by_patient.csv",
  stringsAsFactors = FALSE
)

cat("\n===== DATOS =====\n")
cat("Pacientes SBS96:", nrow(x), "\n")
cat("Columnas:", ncol(x), "\n")

# ==========================================================
# 2. Elegir tres pacientes de prueba
#    - mayor número de SNV
#    - cercano a la mediana
#    - cercano a 50 SNV
# ==========================================================

median_snv <- median(counts$total_snvs)

id_max <- counts$case_id[
  which.max(counts$total_snvs)
]

id_med <- counts$case_id[
  which.min(abs(counts$total_snvs - median_snv))
]

id_50 <- counts$case_id[
  which.min(abs(counts$total_snvs - 50))
]

pilot_ids <- unique(
  c(id_max, id_med, id_50)
)

cat("\n===== PACIENTES PILOTO =====\n")

print(
  counts[
    match(pilot_ids, counts$case_id),
    c("case_id", "total_snvs")
  ]
)

# ==========================================================
# 3. Crear matriz 96 contextos × muestras
# ==========================================================

context_cols <- setdiff(
  colnames(x),
  c("case_id", "tumor_sample")
)

pilot <- x[
  x$case_id %in% pilot_ids,
  c("case_id", context_cols)
]

mat <- as.matrix(
  pilot[, context_cols]
)

storage.mode(mat) <- "numeric"

rownames(mat) <- pilot$case_id

# Transponer:
# filas = 96 contextos
# columnas = pacientes
mat <- t(mat)

cat("\n===== MATRIZ PARA sig_fit =====\n")
cat("Dimensiones:", dim(mat)[1], "x", dim(mat)[2], "\n")
cat("Filas esperadas: 96\n")

cat("\nPrimeros contextos:\n")
print(head(rownames(mat)))

cat("\nPacientes:\n")
print(colnames(mat))

# ==========================================================
# 4. Firmas reportadas para Cervix-CA / TCGA / WES
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
# 5. Ajuste piloto
# ==========================================================

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

# ==========================================================
# 6. Inspeccionar estructura exacta del resultado
# ==========================================================

cat("\n===== ESTRUCTURA DEL RESULTADO =====\n")
print(names(fit))

cat("\n===== str(fit) =====\n")
str(
  fit,
  max.level = 2
)

cat("\n===== PILOTO TERMINADO =====\n")
