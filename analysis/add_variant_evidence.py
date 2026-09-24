#!/usr/bin/env python3

import csv
from pathlib import Path

moderate_file = Path("results/downstream/moderate_variants_1M.csv")
output_file = Path("results/downstream/moderate_variants_evidence_1M.csv")

vcfs = {
    "1094": Path("results/test_1094_1M/annotation/1094.snpeff.vcf"),
    "755":  Path("results/test_755_1M/annotation/755.snpeff.vcf"),
    "999":  Path("results/test_999_1M/annotation/999.snpeff.vcf"),
}

def parse_ad(value):
    if value in ("", ".", None):
        return None, None

    parts = value.split(",")

    try:
        ref = int(parts[0])
        alt = int(parts[1]) if len(parts) > 1 else 0
        return ref, alt
    except ValueError:
        return None, None


# Leer las variantes MODERATE
with open(moderate_file, newline="") as f:
    moderate_rows = list(csv.DictReader(f))

moderate_keys = {
    (
        str(row["patient_id"]),
        row["chrom"],
        row["pos"],
        row["ref"],
        row["alt"],
    ): row
    for row in moderate_rows
}

evidence = {}


for patient_id, vcf_path in vcfs.items():

    with open(vcf_path) as vcf:

        sample_names = []

        for line in vcf:

            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                sample_names = header[9:]
                continue

            fields = line.rstrip("\n").split("\t")

            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            filt = fields[6]
            info = fields[7]

            key = (patient_id, chrom, pos, ref, alt)

            if filt != "PASS" or key not in moderate_keys:
                continue

            # INFO/TLOD
            tlod = ""

            for item in info.split(";"):
                if item.startswith("TLOD="):
                    tlod = item.split("=", 1)[1]
                    break

            format_keys = fields[8].split(":")

            sample_data = {}

            for name, sample_field in zip(sample_names, fields[9:]):

                values = sample_field.split(":")
                values += ["."] * (len(format_keys) - len(values))

                sample_data[name] = dict(zip(format_keys, values))

            normal_name = f"{patient_id}_N"
            tumor_name = f"{patient_id}_T"

            normal = sample_data[normal_name]
            tumor = sample_data[tumor_name]

            normal_ref, normal_alt = parse_ad(normal.get("AD"))
            tumor_ref, tumor_alt = parse_ad(tumor.get("AD"))

            normal_total = (
                normal_ref + normal_alt
                if normal_ref is not None and normal_alt is not None
                else None
            )

            tumor_total = (
                tumor_ref + tumor_alt
                if tumor_ref is not None and tumor_alt is not None
                else None
            )

            normal_vaf = (
                normal_alt / normal_total
                if normal_total not in (None, 0)
                else None
            )

            tumor_vaf = (
                tumor_alt / tumor_total
                if tumor_total not in (None, 0)
                else None
            )

            evidence[key] = {
                "tlod": tlod,
                "normal_dp": normal.get("DP", ""),
                "normal_ref_reads": normal_ref,
                "normal_alt_reads": normal_alt,
                "normal_vaf": normal_vaf,
                "tumor_dp": tumor.get("DP", ""),
                "tumor_ref_reads": tumor_ref,
                "tumor_alt_reads": tumor_alt,
                "tumor_vaf": tumor_vaf,
            }


output_rows = []

for row in moderate_rows:

    key = (
        str(row["patient_id"]),
        row["chrom"],
        row["pos"],
        row["ref"],
        row["alt"],
    )

    ev = evidence.get(key, {})

    new_row = row.copy()
    new_row.update(ev)

    if isinstance(new_row.get("normal_vaf"), float):
        new_row["normal_vaf"] = round(new_row["normal_vaf"], 3)

    if isinstance(new_row.get("tumor_vaf"), float):
        new_row["tumor_vaf"] = round(new_row["tumor_vaf"], 3)

    output_rows.append(new_row)


fieldnames = [
    "patient_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "gene",
    "consequence",
    "hgvs_c",
    "hgvs_p",
    "tlod",
    "normal_dp",
    "normal_ref_reads",
    "normal_alt_reads",
    "normal_vaf",
    "tumor_dp",
    "tumor_ref_reads",
    "tumor_alt_reads",
    "tumor_vaf",
]


with open(output_file, "w", newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=fieldnames,
        extrasaction="ignore"
    )

    writer.writeheader()
    writer.writerows(output_rows)


print(f"Tabla creada: {output_file}")
print(f"Variantes MODERATE: {len(output_rows)}")
print(f"Variantes con evidencia encontrada: {len(evidence)}")
