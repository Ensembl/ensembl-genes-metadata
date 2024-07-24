process WRITE2DB_METADATA {
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_output/$gca", mode: 'copy'

    input:
    val gca
    path metadata
    path species_tmp

    output:
    val gca, emit: gca
    path "${metadata.baseName}.last_id", emit: last_id
    path species_tmp, emit: species_tmp 

    script:
    """
    python /Users/vianey/Documents/ensembl-genes-metadata/src/python/ensembl/genes/metadata/write2db.py --file-path $metadata
    """
}