import csv
import gzip
from pathlib import Path
from collections import Counter

samplesheet = "samplesheet/tcga_cesc_samplesheet.csv"

outdir = Path("results/tcga")
outdir.mkdir(parents=True, exist_ok=True)

NON_SYNONYMOUS = {
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "Splice_Site",
    "In_Frame_Del",
    "In_Frame_Ins",
    "Nonstop_Mutation",
    "Translation_Start_Site"
}

with open(samplesheet, newline="", encoding="utf-8") as f:
    samples = list(csv.DictReader(f))

results = []

for s in samples:

    maf = s["maf_file"]
    case_id = s["case_id"]

    opener = gzip.open if maf.endswith(".gz") else open

    total_maf = 0
    nonsyn_count = 0
    classes = Counter()

    with opener(maf, "rt", encoding="utf-8", errors="replace") as fh:

        header = None

        for line in fh:

            if line.startswith("#"):
                continue

            if header is None:
                header = line.rstrip("\n").split("\t")

                class_idx = header.index("Variant_Classification")
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) <= class_idx:
                continue

            vc = fields[class_idx]

            total_maf += 1
            classes[vc] += 1

            if vc in NON_SYNONYMOUS:
                nonsyn_count += 1

    results.append({
        "case_id": case_id,
        "tumor_grade": s["tumor_grade"],
        "histology_group": s["histology_group"],
        "grade_analysis": s["grade_analysis"],
        "histology_analysis": s["histology_analysis"],
        "total_maf_records": total_maf,
        "mutation_burden": nonsyn_count
    })

outfile = outdir / "tcga_cesc_mutation_burden.csv"

with open(outfile, "w", newline="", encoding="utf-8") as f:

    fieldnames = [
        "case_id",
        "tumor_grade",
        "histology_group",
        "grade_analysis",
        "histology_analysis",
        "total_maf_records",
        "mutation_burden"
    ]

    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(results)

print("Archivo creado:", outfile)
print("Pacientes:", len(results))
print("Mutaciones no sinónimas totales:",
      sum(r["mutation_burden"] for r in results))

burdens = [r["mutation_burden"] for r in results]

print("Mínimo:", min(burdens))
print("Máximo:", max(burdens))
print("Promedio:", round(sum(burdens)/len(burdens), 2))

print("\nPrimeros 10 pacientes:")
for r in results[:10]:
    print(
        r["case_id"],
        r["tumor_grade"],
        r["histology_group"],
        r["mutation_burden"]
    )
