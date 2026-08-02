# BAM2CRAM

## Process Details

| Property | Value |
|----------|-------|
| Process | `BAM2CRAM` |
| Label | `'samtools'` |
| Tag | `$aligned_file` |
| Publish directory | `"${meta.output_dir}", mode: "copy"` |

## Inputs

### Nextflow interface

```nextflow
//tuple val(taxon_id), val(genomeDir), val(tissue), val(platform), val(output_dir), path(aligned_file)
tuple val(meta), path(aligned_file)
```

## Outputs

### Nextflow interface

```nextflow
tuple val(meta), path("*.cram"), emit:cram_output
val("versions.yml"), emit: versions_file
```

## Implementation Summary

- Create symbolic links
- Generate software version report

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/short_read_alignment/modules/bam2cram.nf`
