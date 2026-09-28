# STAR_INDEX_GENOME

STAR_INDEX_GENOME
Generate STAR genome index for RNA-seq alignment.
Input:
- meta: metadata map containing taxon_id, gca, fasta_file, etc.
- statsJson: path to the JSON file containing genome statistics
Output:
- Genome index directory
- Software versions
The module uses STAR to generate the genome index and creates a symbolic link to the output index directory.

## Process Details

| Property | Value |
|----------|-------|
| Process | `STAR_INDEX_GENOME` |
| Label | `'star'` |
| Tag | `${meta.taxon_id}:${meta.gca}` |
| Publish directory | `"${meta.genome_dir}", mode: 'copy'` |
| maxForks | `1` |

## Inputs

### Nextflow interface

```nextflow
tuple val(meta),path(statsJson)
```

## Outputs

### Nextflow interface

```nextflow
val(meta), emit: genome_index_output
path "versions.yml", emit: versions_file
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/short_read_alignment/modules/star_index_genome.nf`
