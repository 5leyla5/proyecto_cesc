import csv
import gzip
from collections import Counter
from pathlib import Path

MAP_FILE = "results/tcga/maf_patient_map_unique.tsv"

OUT_OVERALL = "results/tcga/snv_spectrum_overall.csv"
OUT_PATIENT = "results/tcga/titv_by_patient.csv"

VALID = {"A", "C", "G", "T"}

COMPLEMENT = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C"
}

SBS6_ORDER = [
    "C>A",
    "C>G",
    "C>T",
    "T>A",
    "T>C",
    "T>G"
]


def canonical_sbs(ref, alt):
    """
    Convierte sustituciones a las 6 clases SBS estándar,
    expresadas desde C o T.
    """
    ref = ref.upper()
    alt = alt.upper()

    if ref not in VALID or alt not in VALID:
        return None

    if ref == alt:
        return None

    # Si REF es purina (A/G), usar la hebra complementaria
    if ref in {"A", "G"}:
        ref = COMPLEMENT[ref]
        alt = COMPLEMENT[alt]

    change = f"{ref}>{alt}"

    if change not in SBS6_ORDER:
        return None

    return change


def is_transition(ref, alt):
    pair = f"{ref}>{alt}"

    return pair in {
        "A>G",
        "G>A",
        "C>T",
        "T>C"
    }


# ----------------------------------------------------------
# Leer mapa de pacientes únicos
# ----------------------------------------------------------

samples = []

with open(MAP_FILE, encoding="utf-8") as f:
    reader = csv.reader(f, delimiter="\t")

    for row in reader:
        if len(row) < 3:
            continue

        patient, sample, maf_file = row[:3]

        samples.append({
            "case_id": patient,
            "tumor_sample": sample,
            "maf_file": maf_file
        })


print("Pacientes en mapa:", len(samples))


# ----------------------------------------------------------
# Procesar MAF
# ----------------------------------------------------------

overall_sbs6 = Counter()

total_ti = 0
total_tv = 0
total_snvs = 0

patient_rows = []


for s in samples:

    maf_file = s["maf_file"]

    patient_sbs6 = Counter()

    patient_ti = 0
    patient_tv = 0
    patient_snvs = 0

    with gzip.open(
        maf_file,
        "rt",
        encoding="utf-8",
        errors="replace"
    ) as fh:

        header = None

        for line in fh:

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if header is None:
                header = fields

                ref_idx = header.index("Reference_Allele")
                alt_idx = header.index("Tumor_Seq_Allele2")

                continue

            if len(fields) <= max(ref_idx, alt_idx):
                continue

            ref = fields[ref_idx].upper()
            alt = fields[alt_idx].upper()

            # Solo sustituciones simples de una base
            if (
                len(ref) != 1
                or len(alt) != 1
                or ref not in VALID
                or alt not in VALID
                or ref == alt
            ):
                continue

            patient_snvs += 1
            total_snvs += 1

            if is_transition(ref, alt):
                patient_ti += 1
                total_ti += 1
            else:
                patient_tv += 1
                total_tv += 1

            sbs_class = canonical_sbs(ref, alt)

            if sbs_class:
                patient_sbs6[sbs_class] += 1
                overall_sbs6[sbs_class] += 1


    titv = (
        patient_ti / patient_tv
        if patient_tv > 0
        else None
    )

    row = {
        "case_id": s["case_id"],
        "tumor_sample": s["tumor_sample"],
        "total_snvs": patient_snvs,
        "transitions": patient_ti,
        "transversions": patient_tv,
        "titv_ratio": round(titv, 4) if titv is not None else ""
    }

    for cls in SBS6_ORDER:
        row[cls] = patient_sbs6[cls]

    patient_rows.append(row)


# ----------------------------------------------------------
# Guardar tabla por paciente
# ----------------------------------------------------------

Path("results/tcga").mkdir(parents=True, exist_ok=True)

patient_fields = [
    "case_id",
    "tumor_sample",
    "total_snvs",
    "transitions",
    "transversions",
    "titv_ratio"
] + SBS6_ORDER

with open(
    OUT_PATIENT,
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=patient_fields
    )

    writer.writeheader()
    writer.writerows(patient_rows)


# ----------------------------------------------------------
# Guardar espectro global
# ----------------------------------------------------------

overall_rows = []

for cls in SBS6_ORDER:

    count = overall_sbs6[cls]

    pct = (
        100 * count / total_snvs
        if total_snvs > 0
        else 0
    )

    overall_rows.append({
        "substitution_class": cls,
        "count": count,
        "frequency_pct": round(pct, 2)
    })


with open(
    OUT_OVERALL,
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "substitution_class",
            "count",
            "frequency_pct"
        ]
    )

    writer.writeheader()
    writer.writerows(overall_rows)


# ----------------------------------------------------------
# Resumen
# ----------------------------------------------------------

overall_titv = (
    total_ti / total_tv
    if total_tv > 0
    else None
)

print("\n===== ESPECTRO SNV TCGA-CESC =====")

print("Pacientes procesados:", len(patient_rows))
print("SNV totales:", total_snvs)

print("\n===== Ti / Tv =====")
print("Transitions:", total_ti)
print("Transversions:", total_tv)

if overall_titv is not None:
    print("Ti/Tv:", round(overall_titv, 4))

print("\n===== SBS6 =====")

for cls in SBS6_ORDER:

    count = overall_sbs6[cls]

    pct = (
        100 * count / total_snvs
        if total_snvs > 0
        else 0
    )

    print(
        f"{cls:4s}",
        f"{count:8d}",
        f"{pct:6.2f}%"
    )

print("\nArchivos creados:")
print(OUT_OVERALL)
print(OUT_PATIENT)
