# DOWNLOAD_FASTQS

This process downloads the FASTQ files for a given run accession.
It checks if the files already exist in the specified output directory and skips the download if they do.
The downloaded FASTQ files are saved with the names "<run_accession>_1.fastq.gz" and "<run_accession>_2.fastq.gz" (if paired-end).
The process also creates symbolic links to the downloaded files in the current working directory for easy access.
The process takes a tuple containing metadata information as input and outputs a tuple containing the same metadata along with the paths to the downloaded FASTQ files.

## Process Details

| Property | Value |
|----------|-------|
| Process | `DOWNLOAD_FASTQS` |
| Label | `python` |
| Tag | `${meta.taxonId}:${meta.run_accession}` |
| maxForks | `25` |

## Inputs

### Nextflow interface

```nextflow
tuple val(meta)
//tuple val(taxon_id), val(gca), val(platform), val(paired), val(tissue), val(run_accession), val(genomeDir),  val(url1), val(md5_1), val(url2),  val(md5_2)
```

## Outputs

### Nextflow interface

```nextflow
tuple val(meta), path("*_1.fastq.gz"), path("*_2.fastq.gz", optional: true) , emit: fastq_file_output
path "versions.yml", emit: versions_file
```

## Implementation Summary

- Create symbolic links
- Generate software version report

## Source

`/Users/ftricomi/Downloads/ensembl-genes-metadata/pipelines/short_read_alignment/modules/download_fastqs.nf`
