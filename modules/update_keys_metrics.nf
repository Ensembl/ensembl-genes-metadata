process UPDATE_KEYS_METRICS {
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_output/$gca", mode: 'copy'
    
    input:
    val gca
    path metrics_tmp
    path last_id
    path species_tmp
    
    output:
    val gca, emit: gca
    path "${metrics_tmp.baseName}.json", emit: metrics
    path species_tmp, emit: species_tmp 

    script:
    """
    python /Users/vianey/Documents/ensembl-genes-metadata/src/python/ensembl/genes/metadata/update_keys.py --json-path $metrics_tmp --file-id-path $last_id
    """
}
