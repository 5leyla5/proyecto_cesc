import csv
from pathlib import Path

OUT = Path(
    "results/tcga/driver_passenger_interpretation.csv"
)

rows = [
    {
        "gene": "PIK3CA",
        "mutated_patients": 81,
        "frequency_pct": 28.32,
        "category": "Driver conocido",
        "role": "Oncogene",
        "evidence": (
            "Gen significativamente mutado en TCGA cervical; "
            "oncogene reconocido."
        ),
        "interpretation": (
            "Alta recurrencia y evidencia independiente apoyan "
            "relevancia biológica. La oncogenicidad debe evaluarse "
            "a nivel de variante."
        )
    },

    {
        "gene": "EP300",
        "mutated_patients": 36,
        "frequency_pct": 12.59,
        "category": "Driver conocido",
        "role": "Tumor suppressor / chromatin regulation",
        "evidence": (
            "Gen significativamente mutado en cáncer cervical TCGA."
        ),
        "interpretation": (
            "La recurrencia es compatible con selección tumoral "
            "y alteración de regulación transcripcional/cromatina."
        )
    },

    {
        "gene": "FBXW7",
        "mutated_patients": 34,
        "frequency_pct": 11.89,
        "category": "Driver conocido",
        "role": "Tumor suppressor",
        "evidence": (
            "Gen significativamente mutado en cáncer cervical TCGA."
        ),
        "interpretation": (
            "Su recurrencia y función supresora apoyan un papel "
            "driver en una fracción de los tumores."
        )
    },

    {
        "gene": "KMT2C",
        "mutated_patients": 54,
        "frequency_pct": 18.88,
        "category": "Candidato / gen asociado a cáncer",
        "role": "Chromatin regulation",
        "evidence": (
            "Recurrente en nuestra cohorte y descrito entre genes "
            "mutados candidatos en estudios de carcinoma cervical."
        ),
        "interpretation": (
            "Biológicamente plausible, pero la frecuencia por sí sola "
            "no demuestra que cada evento sea driver."
        )
    },

    {
        "gene": "KMT2D",
        "mutated_patients": 39,
        "frequency_pct": 13.64,
        "category": "Candidato / gen asociado a cáncer",
        "role": "Chromatin regulation",
        "evidence": (
            "Recurrente en nuestra cohorte y descrito en estudios "
            "recientes de genes drivers de carcinoma cervical."
        ),
        "interpretation": (
            "Debe interpretarse con apoyo funcional y a nivel "
            "de variante."
        )
    },

    {
        "gene": "TTN",
        "mutated_patients": 85,
        "frequency_pct": 29.72,
        "category": "Posible pasajero / interpretar con cautela",
        "role": "Very large structural gene",
        "evidence": (
            "Gen extremadamente grande; la frecuencia puede aumentar "
            "por longitud y carga mutacional."
        ),
        "interpretation": (
            "Ser el gen más frecuente no permite concluir que sea "
            "el principal driver de CESC."
        )
    },

    {
        "gene": "MUC16",
        "mutated_patients": 49,
        "frequency_pct": 17.13,
        "category": "Posible pasajero / interpretar con cautela",
        "role": "Large mucin gene",
        "evidence": (
            "Gen de gran longitud; su recurrencia puede estar "
            "influida por tamaño génico."
        ),
        "interpretation": (
            "La recurrencia debe interpretarse con cautela y no "
            "equipararse automáticamente a selección positiva."
        )
    },

    {
        "gene": "FLG",
        "mutated_patients": 36,
        "frequency_pct": 12.59,
        "category": "Interpretación cautelosa",
        "role": "Not established as major cervical SMG",
        "evidence": (
            "Frecuente en nuestra cohorte, pero no identificado entre "
            "los principales genes significativamente mutados de TCGA cervical."
        ),
        "interpretation": (
            "Puede reflejar eventos pasajeros o de relevancia incierta; "
            "requiere evidencia adicional."
        )
    },

    {
        "gene": "SYNE1",
        "mutated_patients": 35,
        "frequency_pct": 12.24,
        "category": "Interpretación cautelosa",
        "role": "Large structural gene",
        "evidence": (
            "Alta recurrencia observacional sin evidencia suficiente "
            "en este análisis para asignarlo como driver."
        ),
        "interpretation": (
            "No debe considerarse driver solo por frecuencia."
        )
    },

    {
        "gene": "DMD",
        "mutated_patients": 35,
        "frequency_pct": 12.24,
        "category": "Interpretación cautelosa",
        "role": "Large gene",
        "evidence": (
            "Gen grande y recurrentemente alterado en la cohorte."
        ),
        "interpretation": (
            "La longitud génica puede contribuir a la frecuencia; "
            "no se demuestra selección positiva con este análisis."
        )
    }
]

OUT.parent.mkdir(
    parents=True,
    exist_ok=True
)

with OUT.open(
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "gene",
            "mutated_patients",
            "frequency_pct",
            "category",
            "role",
            "evidence",
            "interpretation"
        ]
    )

    writer.writeheader()
    writer.writerows(rows)


print("\n===== DRIVERS Y POSIBLES PASAJEROS =====\n")

for r in rows:

    print(
        f"{r['gene']:8s} "
        f"{r['frequency_pct']:6.2f}%  "
        f"{r['category']}"
    )

print("\nArchivo creado:")
print(OUT)
