process WRITE2DB_SPECIES {
    tag "Metrics: $gca"
    publishDir "${params.output_dir}/nextflow_output/$gca", mode: 'copy'

    input:
    val gca
    path species

    output:
    val gca, emit: gca
    path "${species.baseName}.last_id", emit: last_id

    script:
    log.info("Executing Python script to populate Species metrics table: $gca")
    """
    python /Users/vianey/Documents/ensembl-genes-metadata/src/python/ensembl/genes/metadata/write2db.py --file-path $species
    """
}