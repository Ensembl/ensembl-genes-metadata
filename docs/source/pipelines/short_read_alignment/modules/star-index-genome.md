# STAR_INDEX_GENOME

## Process Details

| Property | Value |
|----------|-------|
| Process | `STAR_INDEX_GENOME` |
| Label | `'star'` |
| Tag | `${meta.taxon_id}:${meta.gca}` |
| Publish directory | `"${meta.fasta_file.parent}", mode: 'copy'` |
| maxForks | `1` |

## Inputs

### Nextflow interface

```nextflow
//tuple val(taxon_id), val(genomeDir), val(gca), val(platform), val(paired), val(tissue), val(run_accession), val(pair1), val(pair2)
tuple val(meta),path(statsJson)
```

## Outputs

### Nextflow interface

```nextflow
//tuple val(taxon_id), val(genomeDir), val(platform),  val(tissue), val(run_accession), val(pair1), val(pair2)
tuple val(meta), path("${meta.fasta_file.parent}/Genome"), emit: genome_index_output
path "versions.yml", emit: versions_file
when: meta.platform?.toString()?.toLowerCase() == 'illumina'
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/short_read_alignment/modules/star_index_genome.nf`
