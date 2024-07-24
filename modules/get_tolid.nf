process GET_TOLID {
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_output/$gca", mode: 'copy'

    input:
    val gca

    output:
    val gca, emit: gca
    path "${gca}_tolid.json", emit: tolid

    script:
    """
    python /Users/vianey/Documents/ensembl-genes-metadata/src/python/ensembl/genes/metadata/get_tolid.py --accession $gca
    """


}