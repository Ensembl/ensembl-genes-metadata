# STAR

The STAR process is used to align short reads to a reference genome.
STAR documentation https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf
STAR \
--runThreadN
--genomeDir
--readFilesIn
--outSAMtype BAM SortedByCoordinate \
--outSAMstrandField intronMotif \  It adds a tag (XS:A:+ or XS:A:-) tag to each alignment based on the strand of the splice junction. Used by StringTie
--twopassMode Basic   Improves splice junction sensitivity/precision. Allows STAR to learn junctions in the first pass and improve mapping in the second.
--outFilterIntronMotifs RemoveNoncanonicalUnannotated   Filters out introns that are non-canonical (non-GT/AG) and not annotated.
--limitSjdbInsertNsj 2000000 \  Limits the number of splice junctions to be inserted into the genome index. This is useful for large genomes or when there are many splice junctions.
Consider using the following parameters for STAR alignment:
--outFilterType BySJout \  Use spliced junctions to filter alignments use known junctions-keep only those reads that contain junctions that passed filtering into SJ.out.tab
--alignIntronMax 100000 \ filter long spurious introns

## Process Details

| Property | Value |
|----------|-------|
| Process | `STAR` |
| Label | `'star'` |
| Tag | `$meta.run_accession` |
| Publish directory | `"${params.outDir}/$meta.taxon_id/$meta.run_accession/alignment/", mode: 'copy'` |

## Inputs

### Nextflow interface

```nextflow
//tuple val(taxon_id), val(genomeDir), val(platform), val(tissue), val(run_accession), val(pair1), val(pair2)
val(meta)
```

## Outputs

### Nextflow interface

```nextflow
//tuple val(taxon_id), val(genomeDir), val(gca), val(platform), val(paired), val(tissue), val(run_accession), path("*_Aligned.sortedByCoord.out.bam")
//tuple val(taxon_id), val(genomeDir), val(tissue), val(platform), val(run_accession), path("*.bam")
tuple val(meta), val("${params.outDir}/${meta.taxon_id}/${meta.run_accession}/alignment"), path("${meta.run_accession}_Aligned.sortedByCoord.out.bam"), emit: star_output
path "versions.yml", emit: versions_file
```

## Implementation Summary

- Copy output files
- Create symbolic links
- Generate software version report

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/short_read_alignment/modules/star.nf`
