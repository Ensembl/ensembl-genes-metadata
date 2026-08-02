# CHECK_BAM

CHECK_BAM
Validate a BAM file using samtools quickcheck.
Input:
- meta: metadata map containing taxon_id, platform, tissue, etc.
- bamFile: BAM file to validate
Output:
- Validated BAM file (if valid)
- Software versions
The process uses samtools quickcheck to check the integrity of the input BAM file.
If the BAM file is valid, it is passed to the next step; otherwise, an error is raised.

## Process Details

| Property | Value |
|----------|-------|
| Process | `CHECK_BAM` |
| Label | `samtools` |
| Tag | `${meta.tissue}` |

## Inputs

### Nextflow interface

```nextflow
tuple val(meta), path(bamFile)
```

## Outputs

### Nextflow interface

```nextflow
tuple val(meta), path("${meta.tissue}.bam"), emit: good_bam
path "versions.yml", emit: versions_file
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/short_read_alignment/modules/check_bam.nf`
