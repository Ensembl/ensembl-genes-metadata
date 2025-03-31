process UPDATE_KEYS_METADATA {
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_output/$gca", mode: 'copy'
    
    input:
    tuple val(gca), path(metadata_tmp), path(last_id), path(species_tmp)
    
    output:
    tuple val(gca), path("${metadata_tmp.baseName}.json"), path(species_tmp)
    
    script:
    """
    chmod +x $projectDir/../src/python/ensembl/genes/metadata/update_keys.py
    python $projectDir/../src/python/ensembl/genes/metadata/update_keys.py \
    --json-path $metadata_tmp --file-id-path $last_id \
    --config ${params.db_table_conf}
    """
}
