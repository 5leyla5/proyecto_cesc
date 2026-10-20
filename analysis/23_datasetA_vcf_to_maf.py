import csv
from pathlib import Path

# ============================================================
# Entradas
# ============================================================

SAMPLES = {
    "1094": "results/test_1094_1M/annotation/1094.snpeff.vcf",
    "755":  "results/test_755_1M/annotation/755.snpeff.vcf",
    "999":  "results/test_999_1M/annotation/999.snpeff.vcf",
}

OUTDIR = Path("results/downstream/maf")
OUTDIR.mkdir(parents=True, exist_ok=True)

COMBINED_OUT = Path(
    "results/downstream/datasetA_3patients.maf"
)

# ============================================================
# Columnas MAF
# ============================================================

MAF_FIELDS = [
    "Hugo_Symbol",
    "Entrez_Gene_Id",
    "Center",
    "NCBI_Build",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Strand",
    "Variant_Classification",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode",
    "Matched_Norm_Sample_Barcode",
    "Mutation_Status",
    "FILTER",
    "SnpEff_Impact",
    "Transcript_ID",
    "Gene_ID",
    "HGVSc",
    "HGVSp_Short"
]

# ============================================================
# Prioridad de impacto SnpEff
# ============================================================

IMPACT_RANK = {
    "HIGH": 4,
    "MODERATE": 3,
    "LOW": 2,
    "MODIFIER": 1,
    "": 0
}

# ============================================================
# Clasificación MAF desde términos SnpEff
# ============================================================

def maf_classification(annotation, ref, alt):

    terms = annotation.split("&") if annotation else []

    # Prioridad de efectos
    checks = [
        ("stop_gained", "Nonsense_Mutation"),
        ("frameshift_variant",
         "Frame_Shift_Del" if len(ref) > len(alt)
         else "Frame_Shift_Ins"),
        ("splice_acceptor_variant", "Splice_Site"),
        ("splice_donor_variant", "Splice_Site"),
        ("missense_variant", "Missense_Mutation"),
        ("inframe_deletion", "In_Frame_Del"),
        ("conservative_inframe_deletion", "In_Frame_Del"),
        ("disruptive_inframe_deletion", "In_Frame_Del"),
        ("inframe_insertion", "In_Frame_Ins"),
        ("conservative_inframe_insertion", "In_Frame_Ins"),
        ("disruptive_inframe_insertion", "In_Frame_Ins"),
        ("start_lost", "Translation_Start_Site"),
        ("stop_lost", "Nonstop_Mutation"),
        ("synonymous_variant", "Silent"),
        ("splice_region_variant", "Splice_Region"),
        ("5_prime_UTR_variant", "5'UTR"),
        ("3_prime_UTR_variant", "3'UTR"),
        ("intron_variant", "Intron"),
        ("upstream_gene_variant", "5'Flank"),
        ("downstream_gene_variant", "3'Flank"),
        ("intergenic_region", "IGR"),
        ("non_coding_transcript_exon_variant", "RNA"),
    ]

    for effect, maf_class in checks:
        if effect in terms:
            return maf_class

    return "Unknown"


def variant_type(ref, alt):

    if len(ref) == 1 and len(alt) == 1:
        return "SNP"

    if len(ref) > len(alt):
        return "DEL"

    if len(alt) > len(ref):
        return "INS"

    return "ONP"


# ============================================================
# Parsear ANN de SnpEff
# ============================================================

def parse_ann(info):

    ann_value = None

    for item in info.split(";"):
        if item.startswith("ANN="):
            ann_value = item[4:]
            break

    if not ann_value:
        return []

    annotations = []

    for ann in ann_value.split(","):

        f = ann.split("|")

        # ANN estándar SnpEff:
        # 0 Allele
        # 1 Annotation
        # 2 Annotation_Impact
        # 3 Gene_Name
        # 4 Gene_ID
        # 5 Feature_Type
        # 6 Feature_ID
        # 7 Transcript_BioType
        # 8 Rank
        # 9 HGVS.c
        # 10 HGVS.p

        while len(f) < 11:
            f.append("")

        annotations.append({
            "allele": f[0],
            "annotation": f[1],
            "impact": f[2],
            "gene_name": f[3],
            "gene_id": f[4],
            "feature_type": f[5],
            "feature_id": f[6],
            "biotype": f[7],
            "hgvsc": f[9],
            "hgvsp": f[10]
        })

    return annotations


def choose_annotation(annotations):

    if not annotations:
        return None

    # Priorizamos:
    # 1. mayor impacto
    # 2. transcript protein_coding
    ranked = sorted(
        annotations,
        key=lambda x: (
            IMPACT_RANK.get(x["impact"], 0),
            1 if x["biotype"] == "protein_coding" else 0
        ),
        reverse=True
    )

    return ranked[0]


# ============================================================
# Conversión
# ============================================================

all_rows = []

for patient, vcf_file in SAMPLES.items():

    print(f"\nProcesando paciente {patient}")
    print(vcf_file)

    patient_rows = []

    with open(
        vcf_file,
        encoding="utf-8",
        errors="replace"
    ) as fh:

        for line in fh:

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 8:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt_field = fields[4]
            filt = fields[6]
            info = fields[7]

            # Solo variantes PASS
            if filt != "PASS":
                continue

            for alt in alt_field.split(","):

                anns = parse_ann(info)

                # En VCF bialélicos esto normalmente
                # corresponderá directamente al ALT.
                matching = [
                    x for x in anns
                    if x["allele"] == alt
                ]

                if matching:
                    anns = matching

                best = choose_annotation(anns)

                if best is None:
                    best = {
                        "annotation": "",
                        "impact": "",
                        "gene_name": "Unknown",
                        "gene_id": "",
                        "feature_id": "",
                        "hgvsc": "",
                        "hgvsp": ""
                    }

                vtype = variant_type(ref, alt)

                # Coordenadas conservadas desde el VCF.
                # Para SNP: start=end.
                # Para indels: end cubre el alelo REF.
                start = pos
                end = pos + len(ref) - 1

                row = {
                    "Hugo_Symbol":
                        best["gene_name"] or "Unknown",

                    "Entrez_Gene_Id": "",

                    "Center":
                        "Proyecto_CESC",

                    "NCBI_Build":
                        "GRCh38",

                    "Chromosome":
                        chrom,

                    "Start_Position":
                        start,

                    "End_Position":
                        end,

                    "Strand":
                        "+",

                    "Variant_Classification":
                        maf_classification(
                            best["annotation"],
                            ref,
                            alt
                        ),

                    "Variant_Type":
                        vtype,

                    "Reference_Allele":
                        ref,

                    "Tumor_Seq_Allele1":
                        ref,

                    "Tumor_Seq_Allele2":
                        alt,

                    "Tumor_Sample_Barcode":
                        f"{patient}_T",

                    "Matched_Norm_Sample_Barcode":
                        f"{patient}_N",

                    "Mutation_Status":
                        "Somatic",

                    "FILTER":
                        filt,

                    "SnpEff_Impact":
                        best["impact"],

                    "Transcript_ID":
                        best["feature_id"],

                    "Gene_ID":
                        best["gene_id"],

                    "HGVSc":
                        best["hgvsc"],

                    "HGVSp_Short":
                        best["hgvsp"]
                }

                patient_rows.append(row)
                all_rows.append(row)

    outfile = OUTDIR / f"{patient}.maf"

    with outfile.open(
        "w",
        newline="",
        encoding="utf-8"
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=MAF_FIELDS,
            delimiter="\t"
        )

        writer.writeheader()
        writer.writerows(patient_rows)

    print(
        f"PASS convertidas: {len(patient_rows)}"
    )
    print(f"MAF: {outfile}")


# ============================================================
# MAF combinado
# ============================================================

with COMBINED_OUT.open(
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=MAF_FIELDS,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(all_rows)


print("\n======================================")
print("DATASET A → MAF")
print("======================================")

print("Pacientes:", len(SAMPLES))
print("Variantes PASS totales:", len(all_rows))

print("\nPor paciente:")

for patient in SAMPLES:

    n = sum(
        1
        for row in all_rows
        if row["Tumor_Sample_Barcode"] == f"{patient}_T"
    )

    print(f"{patient}: {n}")

print("\nClasificación de variantes:")

from collections import Counter

classes = Counter(
    row["Variant_Classification"]
    for row in all_rows
)

for cls, n in classes.most_common():
    print(f"{cls:25s} {n}")

print("\nTipo de variante:")

types = Counter(
    row["Variant_Type"]
    for row in all_rows
)

for vt, n in types.most_common():
    print(f"{vt:10s} {n}")

print("\nArchivo combinado:")
print(COMBINED_OUT)
