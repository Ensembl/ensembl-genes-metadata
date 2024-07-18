process PARSE_METADATA {
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_output/$gca", mode: 'copy'

    input:
    val gca

    output:
    val gca, emit: gca
    path("${gca}_metadata.json"), emit: metadata
    path("${gca}_metrics_bioproject.tmp"), emit: metrics_tmp
    path("${gca}_species.tmp"), emit: species_tmp

    script:
    """
    python /Users/vianey/Documents/ensembl-genes-metadata/src/python/ensembl/genes/metadata/retrieving_metadata.py --accession $gca
    """
}