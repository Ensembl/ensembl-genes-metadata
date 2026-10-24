# MINIMAP2

Align long-read sequencing data to a reference genome using Minimap2.
The process maps Oxford Nanopore (ONT) and PacBio reads against a
pre-built Minimap2 genome index (.mmi). The alignment profile is
selected automatically based on the sequencing platform specified
in the metadata.

## Process Details

| Property | Value |
|----------|-------|
| Process | `MINIMAP2` |
| Label | `'minimap2'` |
| Tag | `$meta.run_accession` |
| storeDir | `"${meta.alignment_dir}"` |

## Inputs

### Nextflow interface

```nextflow
tuple val(meta), path(minimap_index_file)
```

## Outputs

### Nextflow interface

```nextflow
tuple val(meta), path("*.sam") , emit: minimap_alignment
path "versions.yml", emit: versions_file
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/short_read_alignment/modules/minimap2.nf`
