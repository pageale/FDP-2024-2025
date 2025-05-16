#!/usr/bin/env nextflow

nextflow.enable.dsl=2
mode = params.publish_dir_mode
ICHORCNA_DIR = params.ichordir

samples_ch = Channel
    .fromPath("${params.reads}*", type: 'dir', checkIfExists: true)
    .ifEmpty { exit 1, "No se encontraron subcarpetas en: ${params.reads}" }
    .filter { dir -> dir.getName() != "unclassified" && dir.listFiles().size() >= 100 }  // this line helps to not considerer barcodees that are noise
    // also, it does not process the unclassified reads
    .map { dir -> tuple(dir.getName(), dir.toString() + "/") }   // add "/" at the end

//samples_ch.view() // to visualize if the tupla is well structured

process dorado_alignment {
    tag "$sampleID"
    publishDir "${params.outdir}/bam", mode: "$mode"
    memory '25GB'
    maxForks 2

    input:
         tuple val(sampleID), path(input_folder)

    output:
         tuple val(sampleID), path("${sampleID}.sorted.aligned.bam"), path("${sampleID}.sorted.aligned.bam.bai")

    script:
    """
    echo "Input folder: ${input_folder}"

    nextflow run epi2me-labs/wf-alignment --bam ${input_folder} --references ${params.references.fasta} \
         --out_dir . --depth_coverage false

    rm combined_refs.fasta
    rm combined_refs.mmi
    rm -rf work/
    """
}

process hmm_readCounter {
    tag "$sampleID"
    publishDir "${params.outdir}/wig", mode: "$mode"

    input:
        tuple val(sampleID), path(bam), path(bai)

    output:
        tuple val(sampleID), path("${sampleID}.wig")

    script:
    """
    readCounter --window 1000000 --quality 20 --chromosome "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y" \
         ${bam} > ${sampleID}.wig
    """
}

process ichorcna {
    tag "$sampleID"
    publishDir "${params.outdir}/ichorcna", mode: "$mode"

    input:
        tuple val(sampleID), path(wig)

    output:
        tuple val(sampleID), path("${sampleID}/**")

    script:
    """
    mkdir -p ${sampleID}
    Rscript ${ICHORCNA_DIR}/scripts/runIchorCNA.R \
      --id ${sampleID} \
      --WIG ${wig} \
      --ploidy "c(2,3)" \
      --normal "c(0.5,0.6,0.7,0.8,0.9)" \
      --maxCN 5 \
      --gcWig ${ICHORCNA_DIR}/inst/extdata/gc_hg38_1000kb.wig \
      --mapWig ${ICHORCNA_DIR}/inst/extdata/map_hg38_1000kb.wig \
      --centromere ${ICHORCNA_DIR}/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
      --normalPanel ${ICHORCNA_DIR}/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds \
      --includeHOMD False \
      --chrs "c(1:22, 'X')" \
      --chrTrain "c(1:22)" \
      --estimateNormal True \
      --estimatePloidy True \
      --estimateScPrevalence True \
      --scStates "c(1,3)" \
      --txnE 0.9999 \
      --txnStrength 10000 \
      --outDir ${sampleID}
    """
}

process tf_table {
    tag "tf_table"
    publishDir "${params.outdir}", mode: "$mode"
    
    input:
        val ichorcna_dirs

    output:
        path "${params.tf_table}"

    script:
    """
    mkdir -p params
    cp ${params.outdir}/ichorcna/*/*.params.txt params/
    python3 ${params.tf_script} params/ ${params.tf_table}
    """
}


workflow {
    def alignment_input

    if ( params.basecalling ) {
        alignment_input = samples_ch | dorado_basecalling
    } else {
        alignment_input = samples_ch
    }

    def aligned = alignment_input 
         | dorado_alignment 
         | hmm_readCounter 
         | ichorcna
         | collect 
         | tf_table
}
