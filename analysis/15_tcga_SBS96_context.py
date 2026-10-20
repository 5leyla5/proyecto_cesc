import csv
import gzip
from collections import Counter
from pathlib import Path

# ============================================================
# Archivos
# ============================================================

MAP_FILE = "results/tcga/maf_patient_map_unique.tsv"

FASTA = "data/reference/grch38/Homo_sapiens_assembly38.fasta"
FAI = FASTA + ".fai"

OUT_OVERALL = "results/tcga/SBS96_overall.csv"
OUT_MATRIX = "results/tcga/SBS96_by_patient.csv"

VALID = {"A", "C", "G", "T"}

COMP = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C"
}

SBS6 = [
    ("C", "A"),
    ("C", "G"),
    ("C", "T"),
    ("T", "A"),
    ("T", "C"),
    ("T", "G")
]


# ============================================================
# Acceso directo al FASTA usando el índice .fai
# No requiere cargar los 3.1 GB en memoria
# ============================================================

class IndexedFasta:

    def __init__(self, fasta_file, fai_file):

        self.fasta = open(fasta_file, "rb")
        self.index = {}

        with open(fai_file, encoding="utf-8") as f:

            for line in f:

                fields = line.rstrip("\n").split("\t")

                chrom = fields[0]

                self.index[chrom] = {
                    "length": int(fields[1]),
                    "offset": int(fields[2]),
                    "line_bases": int(fields[3]),
                    "line_width": int(fields[4])
                }

    def base(self, chrom, pos):

        """
        Devuelve una base usando coordenadas 1-based.
        """

        if chrom not in self.index:
            return None

        info = self.index[chrom]

        if pos < 1 or pos > info["length"]:
            return None

        zero = pos - 1

        byte_offset = (
            info["offset"]
            + (zero // info["line_bases"]) * info["line_width"]
            + (zero % info["line_bases"])
        )

        self.fasta.seek(byte_offset)

        return self.fasta.read(1).decode("ascii").upper()

    def trinucleotide(self, chrom, pos):

        left = self.base(chrom, pos - 1)
        center = self.base(chrom, pos)
        right = self.base(chrom, pos + 1)

        if None in (left, center, right):
            return None

        tri = left + center + right

        if any(b not in VALID for b in tri):
            return None

        return tri

    def close(self):
        self.fasta.close()


# ============================================================
# Reverse complement
# ============================================================

def reverse_complement(seq):

    return "".join(
        COMP[b]
        for b in reversed(seq)
    )


# ============================================================
# Construir las 96 categorías estándar
# ============================================================

SBS96_ORDER = []

for ref, alt in SBS6:

    for left in ["A", "C", "G", "T"]:

        for right in ["A", "C", "G", "T"]:

            SBS96_ORDER.append(
                f"{left}[{ref}>{alt}]{right}"
            )


def canonical_context(tri, ref, alt):

    """
    Convierte cualquier sustitución al formato pirimidínico
    estándar C/T usado en SBS96.

    Ejemplo:
        ACA C>T  -> A[C>T]A

    Si REF es A/G, usa el complemento reverso.
    """

    ref = ref.upper()
    alt = alt.upper()
    tri = tri.upper()

    if ref not in VALID or alt not in VALID:
        return None

    if ref == alt:
        return None

    if len(tri) != 3:
        return None

    # Ya está orientada desde C/T
    if ref in {"C", "T"}:

        left = tri[0]
        right = tri[2]

        return f"{left}[{ref}>{alt}]{right}"

    # Si REF es A/G, invertir orientación
    tri_rc = reverse_complement(tri)

    ref_c = COMP[ref]
    alt_c = COMP[alt]

    left = tri_rc[0]
    right = tri_rc[2]

    return f"{left}[{ref_c}>{alt_c}]{right}"


# ============================================================
# Leer pacientes únicos
# ============================================================

samples = []

with open(MAP_FILE, encoding="utf-8") as f:

    reader = csv.reader(f, delimiter="\t")

    for row in reader:

        if len(row) < 3:
            continue

        samples.append({
            "case_id": row[0],
            "tumor_sample": row[1],
            "maf_file": row[2]
        })


print("Pacientes en mapa:", len(samples))


# ============================================================
# Inicializar FASTA
# ============================================================

reference = IndexedFasta(
    FASTA,
    FAI
)


# ============================================================
# Procesar MAF
# ============================================================

overall = Counter()
patient_rows = []

total_snvs = 0
valid_contexts = 0

ref_mismatch = 0
invalid_context = 0
missing_chrom = 0


for n, sample in enumerate(samples, start=1):

    patient_counts = Counter()

    maf_file = sample["maf_file"]

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

                chr_idx = header.index("Chromosome")
                pos_idx = header.index("Start_Position")
                ref_idx = header.index("Reference_Allele")
                alt_idx = header.index("Tumor_Seq_Allele2")

                continue

            if len(fields) <= max(
                chr_idx,
                pos_idx,
                ref_idx,
                alt_idx
            ):
                continue

            chrom = fields[chr_idx]

            try:
                pos = int(fields[pos_idx])
            except ValueError:
                continue

            ref = fields[ref_idx].upper()
            alt = fields[alt_idx].upper()

            # Solo SNV
            if (
                len(ref) != 1
                or len(alt) != 1
                or ref not in VALID
                or alt not in VALID
                or ref == alt
            ):
                continue

            total_snvs += 1

            if chrom not in reference.index:
                missing_chrom += 1
                continue

            tri = reference.trinucleotide(
                chrom,
                pos
            )

            if tri is None:
                invalid_context += 1
                continue

            # Control crítico:
            # base central del FASTA debe coincidir con REF del MAF
            if tri[1] != ref:
                ref_mismatch += 1
                continue

            context = canonical_context(
                tri,
                ref,
                alt
            )

            if context not in SBS96_ORDER:
                invalid_context += 1
                continue

            valid_contexts += 1

            overall[context] += 1
            patient_counts[context] += 1


    row = {
        "case_id": sample["case_id"],
        "tumor_sample": sample["tumor_sample"]
    }

    for context in SBS96_ORDER:
        row[context] = patient_counts[context]

    patient_rows.append(row)

    if n % 50 == 0:
        print(
            f"Procesados {n}/{len(samples)} pacientes..."
        )


reference.close()


# ============================================================
# Tabla SBS96 global
# ============================================================

overall_rows = []

for context in SBS96_ORDER:

    count = overall[context]

    pct = (
        100 * count / valid_contexts
        if valid_contexts > 0
        else 0
    )

    # Extraer clase central
    substitution = context[2:5]

    overall_rows.append({
        "context": context,
        "substitution_class": substitution,
        "count": count,
        "frequency_pct": round(pct, 4)
    })


Path("results/tcga").mkdir(
    parents=True,
    exist_ok=True
)


with open(
    OUT_OVERALL,
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "context",
            "substitution_class",
            "count",
            "frequency_pct"
        ]
    )

    writer.writeheader()
    writer.writerows(overall_rows)


# ============================================================
# Matriz paciente × SBS96
# ============================================================

matrix_fields = [
    "case_id",
    "tumor_sample"
] + SBS96_ORDER


with open(
    OUT_MATRIX,
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=matrix_fields
    )

    writer.writeheader()
    writer.writerows(patient_rows)


# ============================================================
# Resumen QC
# ============================================================

print("\n===== SBS96 TCGA-CESC =====")

print("Pacientes procesados:", len(patient_rows))
print("SNV totales evaluados:", total_snvs)
print("Contextos SBS96 válidos:", valid_contexts)

print("\n===== CONTROL DE REFERENCIA =====")
print("REF MAF != FASTA:", ref_mismatch)
print("Contextos inválidos:", invalid_context)
print("Cromosomas ausentes en FASTA:", missing_chrom)

if total_snvs > 0:

    pct_valid = 100 * valid_contexts / total_snvs

    print(
        "SNV con contexto válido:",
        round(pct_valid, 2),
        "%"
    )


print("\n===== TOP 15 CONTEXTOS =====")

for context, count in overall.most_common(15):

    pct = (
        100 * count / valid_contexts
        if valid_contexts > 0
        else 0
    )

    print(
        f"{context:10s}",
        f"{count:8d}",
        f"{pct:6.2f}%"
    )


print("\nArchivos creados:")
print(OUT_OVERALL)
print(OUT_MATRIX)
