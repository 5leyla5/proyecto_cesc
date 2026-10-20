library(maftools)

maf_file <- "results/downstream/datasetA_3patients.maf"

cat("\n===== CARGANDO MAF DATASET A =====\n")

maf <- read.maf(
  maf = maf_file,
  isTCGA = FALSE,
  verbose = TRUE
)

cat("\n===== VALIDACIÓN =====\n")

cat(
  "Clase del objeto:",
  class(maf),
  "\n"
)

cat(
  "Muestras detectadas:",
  length(getSampleSummary(maf)$Tumor_Sample_Barcode),
  "\n"
)

cat("\n===== RESUMEN POR MUESTRA =====\n")
print(
  getSampleSummary(maf)
)

cat("\n===== RESUMEN DE GENES =====\n")
print(
  head(
    getGeneSummary(maf),
    20
  )
)

cat("\n===== CAMPOS PRINCIPALES =====\n")

print(
  head(
    maf@data[
      ,
      c(
        "Hugo_Symbol",
        "Chromosome",
        "Start_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele2",
        "Variant_Classification",
        "Variant_Type",
        "Tumor_Sample_Barcode"
      )
    ],
    10
  )
)

cat("\n===== VALIDACIÓN COMPLETADA =====\n")
cat(
  "MAF compatible con maftools.\n"
)
