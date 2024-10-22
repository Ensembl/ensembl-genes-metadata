process CUSTOM_GROUP {

    label 'python'
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_output/$gca", mode: 'copy'

    input:
    val gca

    output:
    val gca, emit: gca
    path "${gca}_group.json", emit: group

    script:
    """
    chmod +x $projectDir/../src/python/ensembl/genes/metadata/custom_groups.py
    python $projectDir/../src/python/ensembl/genes/metadata/custom_groups.py \
    --accession $gca --metadata ${params.metadata_params}
    """
}