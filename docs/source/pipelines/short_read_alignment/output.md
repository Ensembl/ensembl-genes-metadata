# Output

The pipeline organises results by **species**, **sequencing platform**, and **tissue**. Depending on the selected workflow options, the following files are produced.

## Directory structure

```text
<outDir>/
└── <taxon_id>/
    └── <platform>/
        └── <tissue> OR <run_accession>/
            └── alignment/
                ├── <sample>.bam
                ├── <sample>.bam.bai (or .bam.csi)
                ├── <sample>.cram
                ├── <sample>.cram.crai
                ├── <sample>.bw
                ├── <sample>_forward_strand.bam
                ├── <sample>_forward_strand.bam.bai (or .csi)
                ├── <sample>_reverse_strand.bam
                ├── <sample>_reverse_strand.bam.bai (or .csi)
                └── versions.yml
```

The exact contents depend on the options enabled when running the pipeline.

## Output files

| File                      | Description                                                                   | Generated when |
| ------------------------- | ----------------------------------------------------------------------------- | -------------- |
| `*.bam`                   | Coordinate-sorted alignment file.                                             | Always         |
| `*.bam.bai` / `*.bam.csi` | BAM index. CSI is generated automatically when BAI indexing is not supported. | Always         |
| `*.cram`                  | CRAM-compressed alignment.                                                    | `--bam2cram`   |
| `*.cram.crai`             | CRAM index.                                                                   | `--bam2cram`   |
| `*_forward_strand.bam`    | Reads mapped to the forward strand.                                           | `--stranded`   |
| `*_reverse_strand.bam`    | Reads mapped to the reverse strand.                                           | `--stranded`   |
| `*.bw`                    | BigWig genome coverage track generated with deepTools `bamCoverage`.          | `--bam2bigWig` |
| `versions.yml`            | Software versions used by each process.                                       | Always         |

## Tissue-level merging

When `--mergeTissue` is enabled, all BAM files belonging to the same combination of:

* taxon ID
* sequencing platform
* tissue

are merged into a single BAM before downstream analyses. The merged BAM is validated and indexed before optional CRAM conversion, stranded BAM generation, or BigWig creation.

## Software versions

Each process generates a `versions.yml` file containing the versions of the software used during execution.

At the end of the workflow, all version files are collected into a single report, providing a complete record of the software used for the analysis.

## Example output

```text
results/
└── 9606/
    └── illumina/
        └── liver/
            └── alignment/
                ├── liver.bam
                ├── liver.bam.bai
                ├── liver.cram
                ├── liver.cram.crai
                ├── liver_forward_strand.bam
                ├── liver_forward_strand.bam.bai
                ├── liver_reverse_strand.bam
                ├── liver_reverse_strand.bam.bai
                ├── liver.bw
                └── versions.yml
```
