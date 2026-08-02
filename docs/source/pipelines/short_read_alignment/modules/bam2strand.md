# BAM2STRAND

BAM2STRAND
Split a BAM file into forward and reverse strand BAM files.
Input:
- meta: metadata map containing taxon_id, platform, tissue, etc.
- aligned_file: BAM file to split
Output:
- Forward strand BAM file
- Reverse strand BAM file
- Software versions
The process uses samtools to filter the input BAM file based on the strand information
and generates two separate BAM files for forward and reverse strands.

## Process Details

| Property | Value |
|----------|-------|
| Process | `BAM2STRAND` |
| Label | `'samtools'` |
| Tag | `$aligned_file` |
| storeDir | `"${meta.output_dir}"` |

## Inputs

### Nextflow interface

```nextflow
//tuple val(taxon_id), val(genomeDir), val(tissue), val(platform), val(output_dir), path(aligned_file)
tuple val(meta), path(aligned_file)
```

## Outputs

### Nextflow interface

```nextflow
tuple val(meta), path("*_forward_strand.bam"), path("*_reverse_strand.bam"), emit:aligned_output
path "versions.yml", emit: versions_file
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/short_read_alignment/modules/bam2strand.nf`
