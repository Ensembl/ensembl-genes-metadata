# Short Read Alignment Pipeline

## Quick Links

- **[Input Specification](input.md)**
- **[Output Reference](output.md)**
- **[Parameter Reference](parameters.md)**
- **[Workflow Documentation](workflows/index.md)**
- **[Module Documentation](modules/index.md)**
- **[Troubleshooting](troubleshooting.md)**

---

## Overview

The **Short Read Alignment Pipeline** is a Nextflow DSL2 workflow for aligning RNA sequencing data against reference genomes. It supports both short- and long-read sequencing technologies and produces indexed alignment files with optional downstream processing.

The pipeline automatically:

- Downloads reference genomes from NCBI.
- Downloads sequencing reads from user-provided locations.
- Selects the appropriate aligner based on sequencing platform.
- Validates BAM files after alignment.
- Optionally merges BAM files by tissue.
- Optionally generates stranded BAMs.
- Optionally converts BAM files to CRAM.
- Optionally generates BigWig coverage tracks.
- Collects software versions used during the analysis.

---

## Features

- Nextflow DSL2 modular workflow
- Automatic genome retrieval
- Automatic FASTQ download
- Supports Illumina, PacBio and Oxford Nanopore data
- STAR alignment for Illumina reads
- Minimap2 alignment for long-read sequencing
- BAM integrity validation using `samtools quickcheck`
- Optional tissue-level BAM merging
- Optional stranded BAM generation
- Optional BAM → CRAM conversion
- Optional BigWig generation
- Automatic software version reporting

---

## Workflow

```
                    Samplesheet
                         │
                         ▼
                 FETCH_GENOME
                         │
                         ▼
                DOWNLOAD_FASTQS
                         │
          ┌──────────────┴──────────────┐
          │                             │
          ▼                             ▼
      Illumina                   PacBio / ONT
          │                             │
          ▼                             ▼
 STAR_INDEX_GENOME          MINIMAP2_INDEX_GENOME
          │                             │
          ▼                             ▼
         STAR                      MINIMAP2
          │                             │
          └──────────────┬──────────────┘
                         ▼
                    CHECK_BAM
                         │
                         ▼
                  DELETE_FASTQ
                         │
                         ▼
                    INDEX_BAM
                         │
          Optional tissue merging
                         │
                         ▼
                MERGE_BAM_PER_TISSUE
                         │
                         ▼
                    CHECK_BAM
                         │
                         ▼
                    INDEX_BAM
                         │
         ┌───────────────┼────────────────┐
         │               │                │
         ▼               ▼                ▼
   BAM2STRAND       BAM2CRAM        BAM2BIGWIG
```

---

## Supported sequencing platforms

| Platform | Aligner |
|----------|---------|
| Illumina | STAR |
| PacBio | Minimap2 |
| Oxford Nanopore (ONT) | Minimap2 |

---


## Running the pipeline

Example:

```bash
nextflow run main.nf \
    --csvFile samples.csv \
    --mergeTissue \
    --stranded \
    --bam2cram \
    --bam2bigWig
```

---

## Optional processing

The following optional steps can be enabled independently.

### Merge BAM files by tissue

```
--mergeTissue
```

Reads from the same tissue are merged into a single alignment.

---

### Generate stranded BAMs

```
--stranded
```

Creates forward and reverse strand BAM files.

---

### Generate BigWig files

```
--bam2bigWig
```

Produces genome coverage tracks using **deepTools bamCoverage**.

---

### Convert BAM to CRAM

```
--bam2cram
```

Compresses alignments into CRAM format and generates CRAI indexes.

---

## BAM validation

Every alignment is validated using

```
samtools quickcheck
```

before being used by downstream processes.

Merged BAM files are also validated before indexing.

---

## Software

Major tools used by the pipeline include:

- STAR
- Minimap2
- samtools
- deepTools (bamCoverage)
- Nextflow

---

## Module Documentation

Detailed documentation is available for every module, including:

- Overview
- Inputs
- Outputs
- Parameters
- Implementation
- Dependencies
- Source

See the [Module Documentation](modules/index.md).

## Software versions

Each module records its software version in a `versions.yml` file.

At the end of the workflow all module versions are merged into a single report.

**Module**

- [COLLECT_SOFTWARE_VERSIONS](modules/collect-software-versions.md)
---

## Error handling

The workflow performs several validation steps:

- input parameter validation
- BAM integrity checks
- automatic genome retrieval
- automatic FASTQ download
- indexed BAM generation
- software version tracking

Corrupted BAM files are detected using `samtools quickcheck` before downstream analyses are performed.

---

## License

Licensed under the Apache License 2.0.