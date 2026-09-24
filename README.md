# Paisaje mutacional del cáncer cervicouterino

Proyecto de **Bioinformática I** orientado al análisis reproducible de
variantes somáticas en cáncer cervicouterino.

## Pregunta biomédica

**¿Qué diferencias presenta el paisaje mutacional del cáncer cervicouterino
según el grado tumoral (G2 vs G3) y el subtipo histológico
(carcinoma escamoso vs adenocarcinoma), en términos de carga mutacional
y genes recurrentemente alterados?**

## Conjuntos de datos

Se utilizan dos conjuntos de datos complementarios.

### Dataset A: E-MTAB-11407 / PRJEB50605

Subconjunto de datos WES de carcinoma escamoso cervicouterino con
muestras tumor-normal pareadas.

Se seleccionaron tres pacientes para el desarrollo y validación del pipeline:

| Paciente | Tumor | Normal |
|---|---|---|
| 1094 | ERR8314271 | ERR8314261 |
| 755 | ERR8314276 | ERR8314265 |
| 999 | ERR8314284 | ERR8314269 |

Cada paciente posee una muestra tumoral y una muestra normal, con archivos
FASTQ paired-end R1 y R2. En total corresponden a **6 muestras y
12 archivos FASTQ**.

Para las pruebas locales de la Entrega 2 se utilizaron subconjuntos
reproducibles de **1 millón de pares de lecturas por muestra (1M)**.

### Dataset B: TCGA-CESC

Cohorte pública de cáncer cervicouterino obtenida desde
**NCI Genomic Data Commons (GDC)**.

- 286 pacientes únicos.
- 287 archivos de mutaciones somáticas en formato MAF descargados.
- Para el análisis se utiliza un único MAF por paciente.
- Datos clínicos asociados para grado tumoral e histología.

Composición principal:

- G2: 126 pacientes.
- G3: 111 pacientes.
- Carcinoma escamoso: 236 pacientes.
- Adenocarcinoma: 46 pacientes.
- Otros: 4 pacientes.

Este conjunto será utilizado posteriormente para el análisis de cohorte,
carga mutacional, genes recurrentemente alterados y comparaciones clínicas.

## Archivos principales

- `samplesheet/tcga_cesc_samplesheet.csv`: relación entre paciente,
  archivo MAF, grado tumoral e histología.
- `samplesheet/tcga_cesc_clinical_clean.csv`: información clínica
  consolidada por paciente.
- `samplesheet/gdc_original_samplesheet.tsv`: samplesheet original
  descargado desde GDC.
- `samplesheet/cervical_wes_pipeline_samplesheet.csv`: relación entre
  pacientes, muestras tumor-normal y FASTQ R1/R2.
- `data/tcga_cesc/manifests/`: manifest utilizado para reproducir
  la descarga desde GDC.
- `data/cervical_wes_tumor_normal/metadata/selected_runs.tsv`:
  runs seleccionados para el pipeline.
- `data/cervical_wes_tumor_normal/metadata/selected_fastq_urls.txt`:
  enlaces de descarga de los FASTQ seleccionados.

Para las ejecuciones de prueba de 1M se incluyen además:

- `samplesheet/cervical_wes_1094_1M_samplesheet.csv`
- `samplesheet/cervical_wes_755_1M_samplesheet.csv`
- `samplesheet/cervical_wes_999_1M_samplesheet.csv`

Los archivos genómicos de gran tamaño no se almacenan en GitHub;
el repositorio conserva manifests, metadata, samplesheets y enlaces
necesarios para reproducir su descarga.

## Workflow

![Workflow del proyecto](docs/mermaid-diagram.png)

Para la Entrega 2 se encuentra implementada la rama correspondiente
al procesamiento WES tumor-normal desde archivos FASTQ hasta la
anotación funcional de variantes.

El pipeline está implementado en **Nextflow DSL2** e incluye:

`FastQC → BWA-MEM2 → SAMtools → MarkDuplicates → Mutect2 → FilterMutectCalls → SnpEff`

## Implementación y ambientes

Archivos principales del pipeline:

- `main.nf`
- `nextflow.config`
- `conf/params.yaml`

Ambientes reproducibles:

- `conf/environment_nextflow.yml`
- `conf/environment_gatk4.yml`
- `conf/environment_snpeff.yml`

### Software utilizado

| Herramienta | Versión |
|---|---:|
| Nextflow | 26.04.6 |
| FastQC | 0.12.1 |
| BWA-MEM2 | 2.3 |
| SAMtools | 1.24 |
| seqtk | 1.5 |
| GATK | 4.6.2.0 |
| Picard | 3.4.0 |
| SnpEff | 5.4.0c |

La anotación funcional utiliza la base:

`GRCh38.mane.1.5.ensembl`

## Genoma de referencia

Para el alineamiento se utiliza **GRCh38**.

Archivos requeridos:

- `Homo_sapiens_assembly38.fasta`
- `Homo_sapiens_assembly38.fasta.fai`
- `Homo_sapiens_assembly38.dict`
- índices de BWA-MEM2

Los archivos de referencia e índices no se almacenan en GitHub
debido a su tamaño.

## Resultados generados por el pipeline

Las ejecuciones de prueba con los tres pares tumor-normal generan
VCF filtrados y posteriormente anotados con SnpEff.

Los VCF anotados finales se encuentran en:

- `results/test_1094_1M/annotation/1094.snpeff.vcf`
- `results/test_755_1M/annotation/755.snpeff.vcf`
- `results/test_999_1M/annotation/999.snpeff.vcf`

La descripción y análisis de las variantes obtenidas se presenta
en el informe de la Entrega 2.

## Reportes de ejecución

Para cada ejecución se conservan los reportes generados por Nextflow.

Ubicación:

- `results/test_1094_1M/nextflow_reports/`
- `results/test_755_1M/nextflow_reports/`
- `results/test_999_1M/nextflow_reports/`

Cada carpeta contiene:

- `report.html`
- `timeline.html`
- `trace.txt`
- `dag.html`

Estos archivos documentan la ejecución, tiempos, recursos utilizados
y dependencias entre procesos.

## Procedencia de los datos

### Dataset A

- BioStudies / ArrayExpress E-MTAB-11407:
  https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11407
- ENA PRJEB50605:
  https://www.ebi.ac.uk/ena/browser/view/PRJEB50605
- Publicación asociada:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC5102491/

### Dataset B

- NCI GDC TCGA-CESC:
  https://portal.gdc.cancer.gov/projects/TCGA-CESC
- GDC Data Portal:
  https://portal.gdc.cancer.gov/

## Consideración importante

El **Dataset A** se utiliza para demostrar el procesamiento reproducible
de pares tumor-normal mediante Nextflow.

Los subconjuntos de 1M se utilizan como datos de prueba del pipeline
y no representan la profundidad completa de los WES originales.

Las conclusiones sobre el paisaje mutacional, la carga mutacional,
los genes recurrentemente alterados y las comparaciones por grado
tumoral e histología se realizarán posteriormente utilizando el
**Dataset B: TCGA-CESC**.