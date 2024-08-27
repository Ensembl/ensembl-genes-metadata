process CUSTOM_GROUP {
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_output/$gca", mode: 'copy'

    input:
    val gca

    output:
    val gca, emit: gca
    path "${gca}_group.json", emit: group

    script:
    """
    python  /Users/vianey/Documents/ensembl-genes-metadata/src/python/ensembl/genes/metadata/custom_groups.py --accession $gca
    """
}