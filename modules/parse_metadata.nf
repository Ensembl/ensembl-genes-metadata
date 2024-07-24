process PARSE_METADATA {
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_output/$gca", mode: 'copy'

    input:
    val gca

    output:
    val gca, emit: gca
    path("${gca}_assembly.json"), emit: assembly
    path("${gca}_metadata.tmp"), emit: metadata_tmp
    path("${gca}_species.tmp"), emit: species_tmp

    script:
    """
    python /Users/vianey/Documents/ensembl-genes-metadata/src/python/ensembl/genes/metadata/retrieving_metadata.py --accession $gca
    """
}