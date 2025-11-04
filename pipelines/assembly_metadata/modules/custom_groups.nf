process CUSTOM_GROUP {
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_output/$gca", mode: 'copy'

    input:
    tuple val(gca), path(last_id)

    output:
    tuple val(gca), path("${gca}_group.json")

    script:
    """
    chmod +x $projectDir/../src/python/ensembl/genes/metadata/custom_groups.py
    python $projectDir/../src/python/ensembl/genes/metadata/custom_groups.py \
    --accession $gca --metadata ${params.metadata_params}
    """
}