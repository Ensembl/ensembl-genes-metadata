process FETCH_GCA {
    tag "taxon:$taxon"

    input:
    val taxon
    val last_update


    output:
    stdout

    script:
    """
    chmod +x $projectDir/../src/python/ensembl/genes/metadata/fetch_new_assemblies.py
    python $projectDir/../src/python/ensembl/genes/metadata/fetch_new_assemblies.py \
    --taxon $taxon --date_update ${params.last_update} --db asm_metadata \
    --registry ${params.registry_params} --metadata ${params.metadata_params} \
    --ncbi ${params.ncbi_params} --ncbi_url ${params.ncbi_url} 
    """

}