import csv
import gzip
from collections import defaultdict
from pathlib import Path

# ============================================================
# Dataset A: VCF anotados con SnpEff
# ============================================================

vcfs = {
    "1094": "results/test_1094_1M/annotation/1094.snpeff.vcf",
    "755":  "results/test_755_1M/annotation/755.snpeff.vcf",
    "999":  "results/test_999_1M/annotation/999.snpeff.vcf",
}

all_pass_genes = defaultdict(set)
functional_genes = defaultdict(set)
gene_impacts = defaultdict(set)

for patient, vcf in vcfs.items():

    opener = gzip.open if vcf.endswith(".gz") else open

    with opener(vcf, "rt", encoding="utf-8", errors="replace") as fh:

        for line in fh:

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 8:
                continue

            filt = fields[6]
            info = fields[7]

            # Solo variantes PASS
            if filt != "PASS":
                continue

            ann_field = None

            for item in info.split(";"):
                if item.startswith("ANN="):
                    ann_field = item[4:]
                    break

            if not ann_field:
                continue

            for annotation in ann_field.split(","):

                parts = annotation.split("|")

                if len(parts) < 4:
                    continue

                consequence = parts[1]
                impact = parts[2]
                gene = parts[3]

                if not gene or gene == ".":
                    continue

                all_pass_genes[gene].add(patient)
                gene_impacts[gene].add(impact)

                # Comparable de forma aproximada con mutaciones
                # no sinónimas utilizadas en TCGA:
                if impact in {"HIGH", "MODERATE"}:
                    functional_genes[gene].add(patient)


# ============================================================
# Dataset B: frecuencias por histología
# ============================================================

tcga_file = "results/tcga/recurrent_genes_by_histology.csv"

tcga = {}

with open(tcga_file, newline="", encoding="utf-8") as f:

    reader = csv.DictReader(f)

    for row in reader:
        tcga[row["gene"]] = row


outdir = Path("results/tcga")
outdir.mkdir(parents=True, exist_ok=True)


# ============================================================
# Función para generar tablas
# ============================================================

def create_table(gene_dict, outfile):

    rows = []

    for gene, patients in gene_dict.items():

        b = tcga.get(gene)

        if b:
            sq_mut = int(b["Squamous_mutated"])
            sq_total = int(b["Squamous_total"])
            sq_freq = float(b["Squamous_frequency_pct"])

            ad_mut = int(b["Adenocarcinoma_mutated"])
            ad_total = int(b["Adenocarcinoma_total"])
            ad_freq = float(b["Adenocarcinoma_frequency_pct"])

        else:
            sq_mut = 0
            sq_total = 236
            sq_freq = 0

            ad_mut = 0
            ad_total = 46
            ad_freq = 0

        rows.append({
            "gene": gene,
            "datasetA_n_patients": len(patients),
            "datasetA_patients": ";".join(sorted(patients)),
            "datasetA_impacts": ";".join(sorted(gene_impacts[gene])),
            "TCGA_Squamous_mutated": sq_mut,
            "TCGA_Squamous_total": sq_total,
            "TCGA_Squamous_frequency_pct": sq_freq,
            "TCGA_Adenocarcinoma_mutated": ad_mut,
            "TCGA_Adenocarcinoma_total": ad_total,
            "TCGA_Adenocarcinoma_frequency_pct": ad_freq
        })

    rows.sort(
        key=lambda x: (
            x["TCGA_Squamous_frequency_pct"],
            x["datasetA_n_patients"]
        ),
        reverse=True
    )

    with open(outfile, "w", newline="", encoding="utf-8") as f:

        fieldnames = list(rows[0].keys()) if rows else [
            "gene",
            "datasetA_n_patients",
            "datasetA_patients",
            "datasetA_impacts",
            "TCGA_Squamous_mutated",
            "TCGA_Squamous_total",
            "TCGA_Squamous_frequency_pct",
            "TCGA_Adenocarcinoma_mutated",
            "TCGA_Adenocarcinoma_total",
            "TCGA_Adenocarcinoma_frequency_pct"
        ]

        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    return rows


all_rows = create_table(
    all_pass_genes,
    outdir / "datasetA_all_PASS_vs_tcga_squamous.csv"
)

functional_rows = create_table(
    functional_genes,
    outdir / "datasetA_functional_vs_tcga_squamous.csv"
)


# ============================================================
# Resumen
# ============================================================

all_overlap = [
    r for r in all_rows
    if r["TCGA_Squamous_mutated"] > 0
]

functional_overlap = [
    r for r in functional_rows
    if r["TCGA_Squamous_mutated"] > 0
]

print("===== DATASET A vs TCGA-CESC SQUAMOUS =====")

print("\nGenes anotados en variantes PASS del Dataset A:",
      len(all_pass_genes))

print("Genes HIGH/MODERATE del Dataset A:",
      len(functional_genes))

print("\nGenes PASS de A también encontrados en TCGA Squamous:",
      len(all_overlap))

print("Genes HIGH/MODERATE de A también encontrados en TCGA Squamous:",
      len(functional_overlap))

print("\n===== COINCIDENCIAS FUNCIONALES =====")

for r in functional_overlap:

    print(
        f'{r["gene"]:15s}',
        f'A={r["datasetA_patients"]:10s}',
        f'TCGA-SCC={r["TCGA_Squamous_mutated"]}/'
        f'{r["TCGA_Squamous_total"]}',
        f'({r["TCGA_Squamous_frequency_pct"]:.2f}%)'
    )

print("\nArchivos creados:")
print(outdir / "datasetA_all_PASS_vs_tcga_squamous.csv")
print(outdir / "datasetA_functional_vs_tcga_squamous.csv")
