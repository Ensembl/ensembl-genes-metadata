process UPDATE_KEYS {
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_work/$gca", mode: 'copy'
    
    input:
    val gca
    path json
    path last_id_json
    
    output:
    val gca, emit: gca
    path "${json.baseName}.json", emit: updated_json

    script:
    log.info("Executing Python script to add dependent key to json: $gca")
    """
    python /Users/vianey/Documents/ensembl-genes-metadata/src/python/ensembl/genes/metadata/update_keys.py --json-path $json --file-id-path $last_id_json
    """
}
