process WRITE2DB_METRICS {
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_output/$gca", mode: 'copy'

    input:
    val gca
    path metrics
    path species_tmp

    output:
    val gca, emit: gca
    path "${metrics.baseName}.last_id", emit: last_id
    path species_tmp, emit: species_tmp 

    script:
    """
    python /Users/vianey/Documents/ensembl-genes-metadata/src/python/ensembl/genes/metadata/write2db.py --file-path $metrics
    """
}