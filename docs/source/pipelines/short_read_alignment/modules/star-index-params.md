# STAR_INDEX_PARAMS

## Process Details

| Property | Value |
|----------|-------|
| Process | `STAR_INDEX_PARAMS` |
| Label | `'python'` |
| Tag | `${meta.taxon_id}:${meta.gca}` |

## Inputs

### Nextflow interface

```nextflow
val(meta)
//when: meta.platform?.toString()?.toLowerCase() == 'illumina'
```

## Outputs

### Nextflow interface

```nextflow
//tuple val(taxon_id), val(genomeDir), val(platform),  val(tissue), val(run_accession), val(pair1), val(pair2)
tuple val(meta), path("stats.json"), emit: genome_stats_output
path "versions.yml", emit: versions_file
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/short_read_alignment/modules/star_index_params.nf`
