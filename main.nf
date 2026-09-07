nextflow.enable.dsl=2

/*
 * Proyecto Bioinformática I
 * Paisaje mutacional del cáncer cervicouterino
 *
 * Dataset de procesamiento:
 * E-MTAB-11407 / PRJEB50605
 *
 * Workflow esperado:
 * FASTQ -> QC -> alineamiento -> BAM -> llamado somático
 * -> filtrado -> VCF anotado
 */

params.samplesheet = "samplesheet/cervical_wes_pipeline_samplesheet.csv"
params.outdir      = "results"


/*
 * PROCESO 1
 * Control de calidad de los FASTQ
 */
process FASTQC {

    tag "${patient_id}_${sample_type}"

    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple val(patient_id),
          val(sample_type),
          path(read1),
          path(read2)

    output:
    path "*_fastqc.html", emit: html
    path "*_fastqc.zip",  emit: zip

    script:
    """
    fastqc \
        --threads 2 \
        ${read1} \
        ${read2}
    """

    stub:
    """
    touch ${read1.simpleName}_fastqc.html
    touch ${read1.simpleName}_fastqc.zip
    touch ${read2.simpleName}_fastqc.html
    touch ${read2.simpleName}_fastqc.zip
    """
}


/*
 * WORKFLOW PRINCIPAL
 */
workflow {

    /*
     * Leer el samplesheet.
     *
     * Cada paciente genera dos entradas:
     * - tumor
     * - normal
     */
    samples_ch = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .flatMap { row ->

            [
                tuple(
                    row.patient_id,
                    "tumor",
                    file(row.tumor_r1),
                    file(row.tumor_r2)
                ),

                tuple(
                    row.patient_id,
                    "normal",
                    file(row.normal_r1),
                    file(row.normal_r2)
                )
            ]
        }


    /*
     * Primera etapa implementada
     */
    FASTQC(samples_ch)


    /*
     * Próximas etapas
     *
     * BWA_MEM2
     *      ↓
     * SAMTOOLS
     *      ↓
     * MUTECT2
     *      ↓
     * FILTER_VARIANTS
     *      ↓
     * ANNOTATION
     */
}
