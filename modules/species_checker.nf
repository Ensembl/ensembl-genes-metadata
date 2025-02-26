process SPECIES_CHECKER {
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_output/$gca", mode: 'copy'
    maxForks 1

    input:
    tuple val(gca), path(species_tmp), path(last_id)
    //val gca
    //path species_tmp

    output:
    tuple val(gca), path("${species_tmp.baseName}.json")
    //val gca, emit: gca
    //path "${species_tmp.baseName}.json", emit: species

    script:
    """
    chmod +x $projectDir/../src/python/ensembl/genes/metadata/species_checker.py
    python $projectDir/../src/python/ensembl/genes/metadata/species_checker.py \
    --json-path $species_tmp --ncbi_url ${params.ncbi_url} --enscode ${params.enscode} \
    --registry ${params.registry_params} --metadata ${params.metadata_params}
    """
}