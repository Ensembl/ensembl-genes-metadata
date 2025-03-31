process SPECIES_CHECKER {
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_output/$gca", mode: 'copy'

    input:
    tuple val(gca), path(species_tmp), path(last_id)

    output:
    tuple val(gca), path("${species_tmp.baseName}.json")

    script:
    """
    chmod +x $projectDir/../src/python/ensembl/genes/metadata/species_checker.py
    python $projectDir/../src/python/ensembl/genes/metadata/species_checker.py \
    --json-path $species_tmp --ncbi_url ${params.ncbi_url} --enscode ${params.enscode}
    """
}