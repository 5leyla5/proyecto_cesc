import csv
import gzip
from collections import Counter
from pathlib import Path

samplesheet = "samplesheet/tcga_cesc_samplesheet.csv"
outdir = Path("results/tcga")
outdir.mkdir(parents=True, exist_ok=True)

with open(samplesheet, newline="", encoding="utf-8") as f:
    samples = list(csv.DictReader(f))

class_counts = Counter()
total_variants = 0
patients_with_variants = 0
patients_without_variants = []

for s in samples:

    maf = s["maf_file"]
    case_id = s["case_id"]

    opener = gzip.open if maf.endswith(".gz") else open

    n_patient = 0

    with opener(maf, "rt", encoding="utf-8", errors="replace") as fh:

        header = None

        for line in fh:

            if line.startswith("#"):
                continue

            if header is None:
                header = line.rstrip("\n").split("\t")

                try:
                    class_idx = header.index("Variant_Classification")
                except ValueError:
                    raise RuntimeError(
                        f"Variant_Classification no encontrada en {maf}"
                    )

                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) <= class_idx:
                continue

            classification = fields[class_idx]

            class_counts[classification] += 1
            total_variants += 1
            n_patient += 1

    if n_patient > 0:
        patients_with_variants += 1
    else:
        patients_without_variants.append(case_id)

outfile = outdir / "variant_classification_counts.csv"

with open(outfile, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["Variant_Classification", "count"])

    for classification, count in class_counts.most_common():
        writer.writerow([classification, count])

print("===== TCGA-CESC MAF =====")
print("Pacientes analizados:", len(samples))
print("Pacientes con variantes:", patients_with_variants)
print("Pacientes sin variantes:", len(patients_without_variants))
print("Total registros MAF:", total_variants)

print("\n===== VARIANT_CLASSIFICATION =====")
for classification, count in class_counts.most_common():
    print(f"{classification}: {count}")

if patients_without_variants:
    print("\nPacientes sin variantes:")
    for p in patients_without_variants:
        print(p)

print("\nArchivo creado:")
print(outfile)
