process WRITE2DB {
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_work/$gca", mode: 'copy'

    input:
    val gca
    path json

    output:
    val gca, emit: gca
    path "${json.baseName}.last_id", emit: last_id

    script:
    log.info("Executing Python script to populate or update DB: $gca")
    """
    python /Users/vianey/Documents/ensembl-genes-metadata/src/python/ensembl/genes/metadata/write2db.py --file-path $json
    """

}