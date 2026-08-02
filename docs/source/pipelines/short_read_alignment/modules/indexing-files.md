# INDEXING_FILES

The process

## Process Details

| Property | Value |
|----------|-------|
| Process | `INDEXING_FILES` |
| Label | `'samtools'` |
| Tag | `${meta.run_accession}` |
| Publish directory | `"${params.outDir}/${meta.taxon_id}/${meta.output_dir}/alignment", mode: 'copy'` |

## Inputs

### Nextflow interface

```nextflow
//tuple val(taxon_id), val(genomeDir),  val(tissue),  path(aligned_file)
//tuple val(taxon_id), val(genomeDir), val(tissue),val(platform),  val(output_dir), path(aligned_file)
tuple val(meta), path(aligned_file)
val extension
```

## Outputs

### Nextflow interface

```nextflow
//tuple val(taxon_id), val(genomeDir), val(tissue), val(platform), val(output_dir), path(aligned_file)
tuple val(meta), path(aligned_file), emit:aligned_output
path "versions.yml", emit: versions_file
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/short_read_alignment/modules/indexing_files.nf`
