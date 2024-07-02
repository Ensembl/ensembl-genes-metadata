process PARSE_METADATA {
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_work/$gca", mode: 'copy'

    input:
    val gca

    output:
    val gca, emit: gca
    path("${gca}_metadata.json"), emit: metadata
    path("${gca}_metrics.tmp"), emit: metrics_tmp
    path("${gca}_species.json"), emit: species

    script:
    log.info("Executing Python script to get metadata for assembly: $gca")
    """
    python /Users/vianey/Documents/ensembl-genes-metadata/src/python/ensembl/genes/metadata/retrieving_metadata.py --accession $gca
    """
}