process FETCH_GCA {
    tag "taxon:$taxon"

    input:
    val taxon
    path last_update

    output:
    stdout

    script:
    """
    python /Users/vianey/Documents/ensembl-genes-metadata/src/python/ensembl/genes/metadata/fetch_new_assemblies.py --taxon $taxon --file-path $last_update
    """

}