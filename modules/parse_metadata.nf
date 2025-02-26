process PARSE_METADATA {
    tag "$gca"
    publishDir "${params.output_dir}/nextflow_output/$gca", mode: 'copy'

    input:
    val gca

    output:
    tuple val(gca), path("${gca}_assembly.json"), path("${gca}_metadata.tmp"), path("${gca}_species.tmp")
    //val gca, emit: gca
    //path("${gca}_assembly.json"), emit: assembly
    //path("${gca}_metadata.tmp"), emit: metadata_tmp
    //path("${gca}_species.tmp"), emit: species_tmp

    script:
    """
    chmod +x $projectDir/../src/python/ensembl/genes/metadata/retrieving_metadata.py
    python $projectDir/../src/python/ensembl/genes/metadata/retrieving_metadata.py \
    --accession $gca --ncbi_url ${params.ncbi_url}
    """
}