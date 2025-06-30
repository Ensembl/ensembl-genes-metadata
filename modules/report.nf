process REPORT {
    publishDir "${params.output_dir}/nextflow_output/", mode: 'copy'

    input:
    path gca_list
    val last_update

    output:
    path "report.txt"
    path "gca_to_run_ncbi.csv"

    script:
    """
    chmod +x $projectDir/../src/python/ensembl/genes/metadata/create_report.py
    python $projectDir/../src/python/ensembl/genes/metadata/create_report.py \
    --file-list $gca_list --metadata ${params.metadata_params} --update-date $last_update
    """
}