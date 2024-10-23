process WRITE2DB_TOLID {
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_output/$gca", mode: 'copy'

    input:
    val gca
    path tolid

    output:
    val gca, emit: gca
    path "${tolid.baseName}.last_id", emit: last_id

    script:
    """
    chmod +x $projectDir/../src/python/ensembl/genes/metadata/write2db.py
    python $projectDir/../src/python/ensembl/genes/metadata/write2db.py \
    --file-path $tolid --update --metadata ${params.metadata_params} --config ${params.db_table_conf}
    """
}