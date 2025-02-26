process GET_TOLID {
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_output/$gca", mode: 'copy'

    input:
    tuple val(gca), path(last_id)

    output:
    tuple val(gca), path("${gca}_tolid.json")
    //path "${gca}_tolid.json", emit: tolid

    script:
    """
    chmod +x $projectDir/../src/python/ensembl/genes/metadata/get_tolid.py
    python $projectDir/../src/python/ensembl/genes/metadata/get_tolid.py \
    --accession $gca --metadata ${params.metadata_params}
    """


}