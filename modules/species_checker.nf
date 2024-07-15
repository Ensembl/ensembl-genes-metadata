process SPECIES_CHECKER {
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_output/$gca", mode: 'copy'

    input:
    val gca
    path species_tmp

    output:
    val gca, emit: gca
    path "${species_tmp.baseName}.json", emit: species

    script:
    log.info("Executing Python script to fill up missing keys in Species JSON file:")
    """
    python /Users/vianey/Documents/ensembl-genes-metadata/src/python/ensembl/genes/metadata/species_checker.py --json-path $species_tmp
    """
}