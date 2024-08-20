process SPECIES_CHECKER {
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_output/$gca", mode: 'copy'
    maxForks 1

    input:
    val gca
    path species_tmp

    output:
    val gca, emit: gca
    path "${species_tmp.baseName}.json", emit: species

    script:
    """
    python /Users/vianey/Documents/ensembl-genes-metadata/src/python/ensembl/genes/metadata/species_checker.py --json-path $species_tmp
    """
}