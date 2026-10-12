import csv
import gzip
from pathlib import Path

samplesheet = "samplesheet/tcga_cesc_samplesheet.csv"
topgenes_file = "results/tcga/recurrent_genes_overall.csv"
outfile = "results/tcga/tcga_oncoplot_top15.csv"

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

# Prioridad cuando un mismo paciente tiene más de una mutación en el mismo gen
priority = {
    "Frame_Shift_Del": 1,
    "Frame_Shift_Ins": 1,
    "Nonsense_Mutation": 2,
    "Nonstop_Mutation": 2,
    "Translation_Start_Site": 2,
    "Splice_Site": 3,
    "In_Frame_Del": 4,
    "In_Frame_Ins": 4,
    "Missense_Mutation": 5
}

# Top 15 genes
with open(topgenes_file, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    top15 = [r["gene"] for r in list(reader)[:15]]

print("Top 15 genes:")
print(", ".join(top15))

with open(samplesheet, newline="", encoding="utf-8") as f:
    samples = list(csv.DictReader(f))

rows = []

for s in samples:

    # Para este oncoplot usamos los dos subtipos principales
    if s["histology_group"] not in {"Squamous", "Adenocarcinoma"}:
        continue

    maf = s["maf_file"]
    opener = gzip.open if maf.endswith(".gz") else open

    gene_events = {}

    with opener(maf, "rt", encoding="utf-8", errors="replace") as fh:

        header = None

        for line in fh:

            if line.startswith("#"):
                continue

            if header is None:
                header = line.rstrip("\n").split("\t")
                gene_idx = header.index("Hugo_Symbol")
                class_idx = header.index("Variant_Classification")
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) <= max(gene_idx, class_idx):
                continue

            gene = fields[gene_idx]
            vc = fields[class_idx]

            if gene not in top15:
                continue

            if vc not in NON_SYNONYMOUS:
                continue

            if gene not in gene_events:
                gene_events[gene] = vc
            else:
                old = gene_events[gene]

                if priority.get(vc, 99) < priority.get(old, 99):
                    gene_events[gene] = vc

    for gene in top15:

        rows.append({
            "case_id": s["case_id"],
            "tumor_grade": s["tumor_grade"],
            "histology_group": s["histology_group"],
            "gene": gene,
            "variant_class": gene_events.get(gene, "No_mutation")
        })

Path("results/tcga").mkdir(parents=True, exist_ok=True)

with open(outfile, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "case_id",
            "tumor_grade",
            "histology_group",
            "gene",
            "variant_class"
        ]
    )
    writer.writeheader()
    writer.writerows(rows)

print("\nPacientes incluidos:",
      len(set(r["case_id"] for r in rows)))

print("Filas creadas:", len(rows))
print("Archivo:", outfile)
