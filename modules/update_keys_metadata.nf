process UPDATE_KEYS_METADATA {
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_output/$gca", mode: 'copy'
    
    input:
    val gca
    path metadata_tmp
    path last_id
    path species_tmp
    
    output:
    val gca, emit: gca
    path "${metadata_tmp.baseName}.json", emit: metadata
    path species_tmp, emit: species_tmp 

    script:
    """
    python /Users/vianey/Documents/ensembl-genes-metadata/src/python/ensembl/genes/metadata/update_keys.py --json-path $metadata_tmp --file-id-path $last_id
    """
}
