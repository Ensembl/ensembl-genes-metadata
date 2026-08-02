# MERGE_BAM_PER_TISSUE

MERGE_BAM_PER_TISSUE
Merge multiple aligned BAM files belonging to the same tissue into a single BAM.
Input:
- meta: metadata map containing taxon_id, platform, tissue, etc.
- bamFiles: list of BAM files to merge
Output:
- Merged tissue-level BAM file
- BAM index (.bai)
- Software versions
The module validates input BAM files using samtools quickcheck before merging
and validates the merged BAM before indexing.

## Process Details

| Property | Value |
|----------|-------|
| Process | `MERGE_BAM_PER_TISSUE` |
| Label | `samtools` |
| Tag | `${meta.tissue}` |
| maxForks | `2` |
| storeDir | `"${params.outDir}/$meta.taxon_id/$meta.platform/$meta.tissue/alignment"` |

## Inputs

### Nextflow interface

```nextflow
//tuple val(taxon_id), val(genomeDir), val(gca), val(platform), val(paired), val(tissue), val(run_accession), path(aligned_file)
//tuple val(taxon_id), val(genomeDir), val(tissue), val(platform), path(bamFiles)
tuple val(meta), path(bamFiles)
```

## Outputs

### Nextflow interface

```nextflow
//tuple val(taxon_id), val(genomeDir), val(tissue), val(platform), \
//val("${params.outDir}/$taxon_id/$platform/$tissue/alignment"),path("${tissue}.bam")
tuple val(meta), val("${params.outDir}/${meta.taxon_id}/${meta.platform}/${meta.tissue}/alignment"), path("${meta.tissue}.bam"), emit: merged_bam
path "versions.yml", emit: versions_file
```

## Implementation Summary

- Create symbolic links
- Generate software version report

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/short_read_alignment/modules/merge_bam_per_tissue.nf`
