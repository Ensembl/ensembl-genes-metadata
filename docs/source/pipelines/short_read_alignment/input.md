# Input


## Parameters



| Parameter           |                              Default                              | Description                                                              |
| ------------------- | :---------------------------------------------------------------: | ------------------------------------------------------------------------ |
| `--csvFile`         |                            **Required**                           | Input metadata CSV file.                                                 |
| `--outDir`          |                            **Required**                           | Output directory for pipeline results.                                   |
| `--bam2cram`        |                               `true`                              | Convert BAM files to CRAM format.                                        |
| `--mergeTissue`     |                               `true`                              | Merge BAM files belonging to the same tissue.                            |
| `--stranded`        |                               `true`                              | Generate forward and reverse stranded BAM files.                         |
| `--bam2bigWig`      |                               `true`                              | Generate BigWig coverage tracks using deepTools.                         |
| `--cleanOutputDir`  |                               `true`                              | Clean the output directory after completion.                             |
| `--cleanCache`      |                              `false`                              | Clean the cache directory after the workflow completes.                  |
| `--ftpBaseUrl`      |                `ftp://ftp.sra.ebi.ac.uk/vol1/fastq`               | Base FTP URL for FASTQ downloads.                                        |
| `--ncbiBaseUrl`     | `https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/` | Base URL for the NCBI Datasets API.                                      |
| `--cacheDir`        |                              `/cache`                             | Directory used to cache downloaded files.                                |
| `--files_latency`   |                                `60`                               | Delay (seconds) after file operations to accommodate filesystem latency. |
| `--max_intron_size` |                              `100000`                             | Maximum intron size used when building the STAR index.                   |
| `--genome_file`     |                                None                               | Override the genome FASTA with a local file.                             |
| `--enscode`         |                                None                               | Optional Ensembl code used by downstream components.                     |

## CSV File

The pipeline requires a CSV samplesheet.

Required columns include:

| Column | Description |
|---------|-------------|
| taxon_id | NCBI Taxonomy identifier |
| gca | Genome assembly accession |
| platform | Sequencing platform |
| tissue | Tissue name |
| run_accession | Run accession |
| pair1 | FASTQ file or URL |
| pair2 | Second FASTQ (paired-end only) |
| md5_1, md5_2 | MD5 checksum|
| genome_file | Optional local genome |


Example:

```csv
taxon_id,gca,platform,paired,tissue,run_accession,pair1,md5_1,pair2,md5_2,genome_file
9606,GCA_000001405.29,illumina,true,liver,SRR123456,file1.fastq.gz,file1.md5,file2.fastq.gz,file2.md5,fasta path
```
