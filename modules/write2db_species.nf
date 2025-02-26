process WRITE2DB_SPECIES {
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_output/$gca", mode: 'copy'
    maxForks 1

    input:
    tuple val(gca), path(species)
    //val gca
    //path species

    output:
    tuple val(gca), path("${species.baseName}.last_id")
    //val gca, emit: gca
    //path "${species.baseName}.last_id", emit: last_id

    script:
    """
    chmod +x $projectDir/../src/python/ensembl/genes/metadata/write2db.py
    python $projectDir/../src/python/ensembl/genes/metadata/write2db.py \
    --file-path $species --metadata ${params.metadata_params} --config ${params.db_table_conf}
    """
}