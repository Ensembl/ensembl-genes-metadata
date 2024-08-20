process WRITE2DB_SPECIES {
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_output/$gca", mode: 'copy'
    maxForks 1

    input:
    val gca
    path species

    output:
    val gca, emit: gca
    path "${species.baseName}.last_id", emit: last_id

    script:
    """
    python /Users/vianey/Documents/ensembl-genes-metadata/src/python/ensembl/genes/metadata/write2db.py --file-path $species
    """
}