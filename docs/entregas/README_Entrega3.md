## Entrega 3 — Análisis preliminar

### Objetivo
Presentar resultados de control de calidad del Dataset A y primeras figuras estadísticas/biológicas del Dataset B (TCGA-CESC), manteniendo reproducibilidad y trazabilidad entre scripts, tablas y figuras.

### Dataset A — E-MTAB-11407 / PRJEB50605
El Dataset A corresponde a carcinoma escamoso cervicouterino con tres pares tumor-normal. Se utilizaron subconjuntos reproducibles de 1 millón de pares de lecturas por muestra (1M) para validar el pipeline técnico.

Resultados principales:
- Mapeo: 96,09%–99,87%.
- Properly paired: 93,01%–98,00%.
- Duplicación destacada: 1094_N = 55,23% y 1094_T = 17,61%.
- Mutect2: 555 variantes candidatas.
- Variantes PASS: 55 en total.
  - 1094: 6 PASS.
  - 755: 15 PASS.
  - 999: 34 PASS.

Los subconjuntos 1M se utilizan únicamente como validación técnica y no como estimación de carga mutacional real.

### Dataset B — TCGA-CESC
Se descargaron 287 archivos MAF correspondientes a 286 pacientes únicos. El paciente TCGA-UC-A7PG presentaba dos muestras tumorales; se conservó el tumor primario 01A y se excluyó la muestra metastásica 06A.

Cohorte:
- G2: 126 pacientes.
- G3: 111 pacientes.
- Squamous: 236 pacientes.
- Adenocarcinoma: 46 pacientes.
- Other: 4 pacientes.

### Carga mutacional preliminar
Se definió como el número de mutaciones somáticas no sinónimas por tumor.

- Total: 49.458 mutaciones no sinónimas.
- G2: mediana = 69; Q1 = 48,2; Q3 = 118.
- G3: mediana = 75; Q1 = 41; Q3 = 137.
- Wilcoxon G2 vs G3: p = 0,808.

Histología:
- Squamous: mediana = 74; Q1 = 49; Q3 = 134.
- Adenocarcinoma: mediana = 63; Q1 = 37,8; Q3 = 122.
- Wilcoxon: p = 0,0686.

No se observaron diferencias estadísticamente significativas en carga mutacional.

### Genes recurrentes
Top genes globales:
- TTN: 29,72%.
- PIK3CA: 28,32%.
- KMT2C: 18,88%.
- MUC16: 17,13%.
- KMT2D: 13,64%.

Comparaciones por gen:
- G2 vs G3: no se identificaron asociaciones robustas después de corrección FDR.
- Histología: KRAS, TP53 y ERBB2 mostraron diferencias nominales, pero FDR = 1.

### Oncoplots
Se generaron oncoplots de los 15 genes más recurrentes:
- Histología: 282 pacientes (236 Squamous + 46 Adenocarcinoma).
- Grado: 237 pacientes (126 G2 + 111 G3).

Archivos:
- `results/tcga/figures/tcga_oncoplot_top15_histology.png`
- `results/tcga/figures/tcga_oncoplot_top15_G2_vs_G3.png`

### Comparación Dataset A vs TCGA-CESC escamoso
- Genes PASS de A también presentes en TCGA escamoso: 38/52 (73,1%).
- Genes HIGH/MODERATE de A también presentes en TCGA escamoso: 17/19 (89,5%).
- Esta comparación es descriptiva a nivel de gen y no constituye validación de variantes exactas.

### Scripts de la Entrega 3
- `analysis/01_qc_datasetA.R`
- `analysis/02_variant_qc_datasetA.R`
- `analysis/03_check_tcga_maf.py`
- `analysis/04_tcga_mutation_burden.py`
- `analysis/05_tcga_mutation_burden_stats.R`
- `analysis/06_tcga_recurrent_genes.py`
- `analysis/07_tcga_recurrent_genes_stats.R`
- `analysis/08_compare_datasetA_vs_tcga_squamous.py`
- `analysis/09_datasetA_vs_tcga_squamous_plot.R`
- `analysis/10_tcga_oncoplot_data.py`
- `analysis/10b_tcga_oncoplot_grade_data.py`
- `analysis/11_tcga_oncoplot.R`
- `analysis/12_tcga_oncoplot_grade.R`

### Interpretación preliminar
Hasta esta entrega no se observan diferencias estadísticamente significativas en carga mutacional entre G2 y G3 ni entre carcinoma escamoso y adenocarcinoma. El paisaje de genes recurrentes muestra señales nominales distintas entre subtipos, especialmente en KRAS, TP53 y ERBB2, pero estas asociaciones no permanecen significativas después de corregir por comparaciones múltiples.
