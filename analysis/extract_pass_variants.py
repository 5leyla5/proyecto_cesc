#!/usr/bin/env python3

import csv
from pathlib import Path

vcfs = {
    "1094": Path("results/test_1094_1M/annotation/1094.snpeff.vcf"),
    "755":  Path("results/test_755_1M/annotation/755.snpeff.vcf"),
    "999":  Path("results/test_999_1M/annotation/999.snpeff.vcf"),
}

output = Path("results/downstream/pass_variants_1M.csv")
output.parent.mkdir(parents=True, exist_ok=True)

impact_rank = {
    "HIGH": 4,
    "MODERATE": 3,
    "LOW": 2,
    "MODIFIER": 1,
    "": 0
}

rows = []

for patient_id, vcf_path in vcfs.items():

    if not vcf_path.exists():
        raise FileNotFoundError(f"No existe: {vcf_path}")

    with open(vcf_path, "r") as vcf:

        for line in vcf:

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            filt = fields[6]
            info = fields[7]

            # Solo variantes que pasaron los filtros
            if filt != "PASS":
                continue

            ann_value = ""

            for item in info.split(";"):
                if item.startswith("ANN="):
                    ann_value = item[4:]
                    break

            annotations = []

            if ann_value:
                for ann in ann_value.split(","):

                    parts = ann.split("|")

                    # ANN de SnpEff posee hasta 16 campos
                    parts += [""] * (16 - len(parts))

                    annotations.append({
                        "allele": parts[0],
                        "consequence": parts[1],
                        "impact": parts[2],
                        "gene": parts[3],
                        "gene_id": parts[4],
                        "feature_type": parts[5],
                        "transcript": parts[6],
                        "biotype": parts[7],
                        "rank": parts[8],
                        "hgvs_c": parts[9],
                        "hgvs_p": parts[10],
                    })

            if annotations:
                best = max(
                    annotations,
                    key=lambda x: impact_rank.get(x["impact"], 0)
                )
            else:
                best = {
                    "gene": "",
                    "consequence": "",
                    "impact": "",
                    "hgvs_c": "",
                    "hgvs_p": "",
                    "transcript": ""
                }

            rows.append({
                "patient_id": patient_id,
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "gene": best["gene"],
                "consequence": best["consequence"],
                "impact": best["impact"],
                "hgvs_c": best["hgvs_c"],
                "hgvs_p": best["hgvs_p"],
                "transcript": best["transcript"],
            })


with open(output, "w", newline="") as csvfile:

    fieldnames = [
        "patient_id",
        "chrom",
        "pos",
        "ref",
        "alt",
        "gene",
        "consequence",
        "impact",
        "hgvs_c",
        "hgvs_p",
        "transcript",
    ]

    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    writer.writerows(rows)


print(f"Tabla creada: {output}")
print(f"Variantes PASS totales: {len(rows)}")

for patient_id in vcfs:
    n = sum(1 for row in rows if row["patient_id"] == patient_id)
    print(f"{patient_id}: {n} PASS")
