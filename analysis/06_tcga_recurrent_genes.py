import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path

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

overall = Counter()

grade_counts = {
    "G2": Counter(),
    "G3": Counter()
}

hist_counts = {
    "Squamous": Counter(),
    "Adenocarcinoma": Counter()
}

grade_n = Counter()
hist_n = Counter()

for s in samples:

    case_id = s["case_id"]
    maf = s["maf_file"]
    grade = s["tumor_grade"]
    hist = s["histology_group"]

    opener = gzip.open if maf.endswith(".gz") else open

    genes_patient = set()

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
            classification = fields[class_idx]

            if (
                classification in NON_SYNONYMOUS
                and gene
                and gene != "Unknown"
            ):
                genes_patient.add(gene)

    # Una vez por gen y paciente
    for gene in genes_patient:
        overall[gene] += 1

    if grade in grade_counts:
        grade_n[grade] += 1
        for gene in genes_patient:
            grade_counts[grade][gene] += 1

    if hist in hist_counts:
        hist_n[hist] += 1
        for gene in genes_patient:
            hist_counts[hist][gene] += 1


# ------------------------------------------------------------
# Overall
# ------------------------------------------------------------

with open(
    outdir / "recurrent_genes_overall.csv",
    "w",
    newline="",
    encoding="utf-8"
) as f:

    w = csv.writer(f)
    w.writerow(["gene", "mutated_patients", "total_patients", "frequency_pct"])

    for gene, n in overall.most_common():
        w.writerow([
            gene,
            n,
            len(samples),
            round(100 * n / len(samples), 2)
        ])


# ------------------------------------------------------------
# Grade
# ------------------------------------------------------------

all_grade_genes = set(grade_counts["G2"]) | set(grade_counts["G3"])

with open(
    outdir / "recurrent_genes_by_grade.csv",
    "w",
    newline="",
    encoding="utf-8"
) as f:

    w = csv.writer(f)

    w.writerow([
        "gene",
        "G2_mutated",
        "G2_total",
        "G2_frequency_pct",
        "G3_mutated",
        "G3_total",
        "G3_frequency_pct"
    ])

    rows = []

    for gene in all_grade_genes:

        g2 = grade_counts["G2"][gene]
        g3 = grade_counts["G3"][gene]

        rows.append([
            gene,
            g2,
            grade_n["G2"],
            round(100 * g2 / grade_n["G2"], 2),
            g3,
            grade_n["G3"],
            round(100 * g3 / grade_n["G3"], 2)
        ])

    rows.sort(
        key=lambda x: x[1] + x[4],
        reverse=True
    )

    w.writerows(rows)


# ------------------------------------------------------------
# Histology
# ------------------------------------------------------------

all_hist_genes = set(hist_counts["Squamous"]) | set(hist_counts["Adenocarcinoma"])

with open(
    outdir / "recurrent_genes_by_histology.csv",
    "w",
    newline="",
    encoding="utf-8"
) as f:

    w = csv.writer(f)

    w.writerow([
        "gene",
        "Squamous_mutated",
        "Squamous_total",
        "Squamous_frequency_pct",
        "Adenocarcinoma_mutated",
        "Adenocarcinoma_total",
        "Adenocarcinoma_frequency_pct"
    ])

    rows = []

    for gene in all_hist_genes:

        sq = hist_counts["Squamous"][gene]
        ad = hist_counts["Adenocarcinoma"][gene]

        rows.append([
            gene,
            sq,
            hist_n["Squamous"],
            round(100 * sq / hist_n["Squamous"], 2),
            ad,
            hist_n["Adenocarcinoma"],
            round(100 * ad / hist_n["Adenocarcinoma"], 2)
        ])

    rows.sort(
        key=lambda x: x[1] + x[4],
        reverse=True
    )

    w.writerows(rows)


print("===== COHORTE =====")
print("Total:", len(samples))
print("G2:", grade_n["G2"])
print("G3:", grade_n["G3"])
print("Squamous:", hist_n["Squamous"])
print("Adenocarcinoma:", hist_n["Adenocarcinoma"])

print("\n===== TOP 20 GENES RECURRENTES =====")

for gene, n in overall.most_common(20):
    print(
        f"{gene:15s}",
        f"{n:3d}/{len(samples)}",
        f"({100*n/len(samples):.1f}%)"
    )

print("\nArchivos creados:")
print(outdir / "recurrent_genes_overall.csv")
print(outdir / "recurrent_genes_by_grade.csv")
print(outdir / "recurrent_genes_by_histology.csv")
