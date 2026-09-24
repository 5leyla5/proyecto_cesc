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

params.samplesheet     = "samplesheet/cervical_wes_pipeline_samplesheet.csv"
params.outdir          = "results"
params.reference_prefix = "data/reference/grch38/Homo_sapiens_assembly38"
params.reference       = "data/reference/grch38/Homo_sapiens_assembly38.fasta"
params.interval = null
params.snpeff_db   = "GRCh38.mane.1.5.ensembl"
params.snpeff_heap = "8g"

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
      val(run_id),
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
 * PROCESO 2
 * Alineamiento contra GRCh38 con BWA-MEM2
 */
process BWA_MEM2 {

    tag "${patient_id}_${sample_type}"

    publishDir "${params.outdir}/alignment", mode: 'copy'

    maxForks 1

    input:
    tuple val(patient_id),
          val(sample_type),
          val(run_id),
          path(read1),
          path(read2)

    path reference_files

    output:
        tuple val(patient_id),
        val(sample_type),
        val(run_id),
        path("*.sorted.bam"),
        emit: bam

    script:

    def sample_name = sample_type == "tumor" ?
        "${patient_id}_T" :
        "${patient_id}_N"

    def fasta = reference_files.find {
        it.name == "Homo_sapiens_assembly38.fasta"
    }

    """
    bwa-mem2 mem \
        -t 1 \
        -R '@RG\\tID:${run_id}\\tSM:${sample_name}\\tLB:${sample_name}\\tPL:ILLUMINA' \
        ${fasta} \
        ${read1} \
        ${read2} \
        | samtools sort \
            -@ 1 \
            -m 256M \
            -o ${sample_name}.sorted.bam
    """

    stub:

    def sample_name = sample_type == "tumor" ?
        "${patient_id}_T" :
        "${patient_id}_N"

    """
    touch ${sample_name}.sorted.bam
    """
}

/*
 * PROCESO 3
 * Indexar los BAM ordenados
 */
process SAMTOOLS_INDEX {

    tag "${patient_id}_${sample_type}"

    publishDir "${params.outdir}/alignment", mode: 'copy'

    input:
    tuple val(patient_id),
          val(sample_type),
          val(run_id),
          path(bam)

    output:
    tuple val(patient_id),
          val(sample_type),
          val(run_id),
          path(bam),
          path("${bam}.bai"),
          emit: indexed_bam

    script:
    """
    samtools index ${bam}
    """

    stub:
    """
    touch ${bam}.bai
    """
}

/*
 * PROCESO 4
 * Marcado de lecturas duplicadas con GATK
 */
process MARK_DUPLICATES {

    tag "${patient_id}_${sample_type}"

    publishDir "${params.outdir}/alignment/markduplicates", mode: 'copy'

    input:
    tuple val(patient_id),
          val(sample_type),
          val(run_id),
          path(bam),
          path(bai)

    output:
    tuple val(patient_id),
          val(sample_type),
          val(run_id),
          path("*.markdup.bam"),
          path("*.markdup.bai"),
          emit: bam

    path "*_markdup_metrics.txt", emit: metrics

    script:

    def sample_name = sample_type == "tumor" ?
        "${patient_id}_T" :
        "${patient_id}_N"

    """
    micromamba run -n gatk4 gatk MarkDuplicates \
        -I ${bam} \
        -O ${sample_name}.markdup.bam \
        -M ${sample_name}_markdup_metrics.txt \
        --CREATE_INDEX true
    """

    stub:

    def sample_name = sample_type == "tumor" ?
        "${patient_id}_T" :
        "${patient_id}_N"

    """
    touch ${sample_name}.markdup.bam
    touch ${sample_name}.markdup.bai
    touch ${sample_name}_markdup_metrics.txt
    """
}

/*
 * PROCESO 5
 * Llamado de variantes somáticas tumor-normal con Mutect2
 */
process MUTECT2 {

    tag "${patient_id}"

    publishDir "${params.outdir}/mutect2", mode: 'copy'

    maxForks 1

    input:
    tuple val(patient_id),
          val(tumor_run),
          path(tumor_bam),
          path(tumor_bai),
          val(normal_run),
          path(normal_bam),
          path(normal_bai)

    path reference_files

    output:
    tuple val(patient_id),
          path("${patient_id}.unfiltered.vcf.gz"),
          path("${patient_id}.unfiltered.vcf.gz.tbi"),
          path("${patient_id}.unfiltered.vcf.gz.stats"),
          emit: vcf

    script:

    def fasta = reference_files.find {
        it.name == "Homo_sapiens_assembly38.fasta"
    }

    def normal_name = "${patient_id}_N"

    def interval_arg = params.interval ?
        params.interval.split(',').collect { "-L ${it.trim()}" }.join(' ') :
        ""

    """
    micromamba run -n gatk4 gatk \
        --java-options "-Xmx8g" \
        Mutect2 \
        -R ${fasta} \
        -I ${tumor_bam} \
        -I ${normal_bam} \
        -normal ${normal_name} \
        ${interval_arg} \
        -O ${patient_id}.unfiltered.vcf.gz
    """

    stub:
    """
    touch ${patient_id}.unfiltered.vcf.gz
    touch ${patient_id}.unfiltered.vcf.gz.tbi
    touch ${patient_id}.unfiltered.vcf.gz.stats
    """
}

/*
 * PROCESO 6
 * Filtrado de variantes somáticas con FilterMutectCalls
 */
process FILTER_MUTECT_CALLS {

    tag "${patient_id}"

    publishDir "${params.outdir}/mutect2/filtered", mode: 'copy'

    maxForks 1

    input:
    tuple val(patient_id),
          path(vcf),
          path(tbi),
          path(stats)

    path reference_files

    output:
    tuple val(patient_id),
          path("${patient_id}.filtered.vcf.gz"),
          path("${patient_id}.filtered.vcf.gz.tbi"),
          emit: vcf

    path "${patient_id}.filtered.vcf.gz.filteringStats.tsv",
         emit: filtering_stats

    script:

    def fasta = reference_files.find {
        it.name == "Homo_sapiens_assembly38.fasta"
    }

    """
    micromamba run -n gatk4 gatk \
        --java-options "-Xmx4g" \
        FilterMutectCalls \
        -R ${fasta} \
        -V ${vcf} \
        --stats ${stats} \
        -O ${patient_id}.filtered.vcf.gz
    """

    stub:
    """
    touch ${patient_id}.filtered.vcf.gz
    touch ${patient_id}.filtered.vcf.gz.tbi
    touch ${patient_id}.filtered.vcf.gz.filteringStats.tsv
    """
}

/*
 * PROCESO 7
 * Anotación funcional de variantes con SnpEff
 */
process SNPEFF {

    tag "${patient_id}"

    publishDir "${params.outdir}/annotation", mode: 'copy'

    maxForks 1

    input:
    tuple val(patient_id),
          path(vcf),
          path(tbi)

    output:
    tuple val(patient_id),
          path("${patient_id}.snpeff.vcf"),
          emit: vcf

    path "${patient_id}.snpeff_summary.html",
         emit: summary

    path "${patient_id}.snpeff_summary.genes.txt",
         emit: genes

    path "${patient_id}.snpeff.log",
         emit: log

    script:
    """
    micromamba run -n snpeff snpEff \
        -Xmx${params.snpeff_heap} \
        -v \
        -stats ${patient_id}.snpeff_summary.html \
        ${params.snpeff_db} \
        ${vcf} \
        > ${patient_id}.snpeff.vcf \
        2> ${patient_id}.snpeff.log
    """

    stub:
    """
    touch ${patient_id}.snpeff.vcf
    touch ${patient_id}.snpeff_summary.html
    touch ${patient_id}.snpeff_summary.genes.txt
    touch ${patient_id}.snpeff.log
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
                    row.tumor_run,
                    file(row.tumor_r1),
                    file(row.tumor_r2)
                ),

                tuple(
                    row.patient_id,
                    "normal",
                    row.normal_run,
                    file(row.normal_r1),
                    file(row.normal_r2)
                )
            ]
        }

    reference_ch = Channel
        .fromPath("${params.reference_prefix}*", checkIfExists: true)
        .collect()


    /*
    * Control de calidad
    */ 
    FASTQC(samples_ch)
    
    /*
    * Alineamiento contra GRCh38
    */
    BWA_MEM2(samples_ch, reference_ch)

    /*
    * Indexación de los BAM ordenados con samtools
    */
    SAMTOOLS_INDEX(BWA_MEM2.out.bam)

    /*
    * Marcado de duplicados con GATK
    */
    MARK_DUPLICATES(SAMTOOLS_INDEX.out.indexed_bam)

    /*
    * Separar muestras tumor y normal
    * y emparejarlas por paciente
    */
    markdup_split = MARK_DUPLICATES.out.bam.branch {

        tumor:
            it[1] == "tumor"

        normal:
            it[1] == "normal"
    }


    /*
    * Reorganizar los canales dejando patient_id
    * como clave para el emparejamiento
    */
    tumor_ch = markdup_split.tumor.map {
        patient_id, sample_type, run_id, bam, bai ->

        tuple(
            patient_id,
            run_id,
            bam,
            bai
        )
    }

    normal_ch = markdup_split.normal.map {
        patient_id, sample_type, run_id, bam, bai ->

        tuple(
            patient_id,
            run_id,
            bam,
            bai
        )
    }

    /*
    * Unir tumor y normal del mismo paciente
    */
    pairs_ch = tumor_ch.join(normal_ch)

    /*
    * Llamado somático tumor-normal
    */
    MUTECT2(pairs_ch, reference_ch)

    /*
    * Filtrado de variantes somáticas
    */
    FILTER_MUTECT_CALLS(MUTECT2.out.vcf, reference_ch)

    /*
     * Anotación funcional con SnpEff
     */
    SNPEFF(FILTER_MUTECT_CALLS.out.vcf)

}
