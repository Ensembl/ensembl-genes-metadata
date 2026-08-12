# Short_Read_Alignment Parameters

Automatically generated from the pipeline `nextflow_schema.json`.

## Input/output options

Define where the pipeline should find input data and save output data.

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `help` | boolean |  | no | Display help text. |
| `bam2cram` | boolean | True | no | Option to convert BAM files to CRAM format. |
| `mergeTissue` | boolean | True | no | Option to merge tissue-specific BAM files into a single file. |
| `stranded` | boolean | True | no | Option to produce stranded BAM files. |
| `bam2bigWig` | boolean | True | no | Option to convert BAM files to BigWig format. |
| `cleanOutputDir` | boolean | True | no | Option to clean the output directory after running the pipeline. |
| `csvFile` | string |  | yes | Path to the input CSV metadata file. |
| `ftpBaseUrl` | string | ftp://ftp.sra.ebi.ac.uk/vol1/fastq | no | Base URL for the FTP server where FASTQ files are stored. |
| `ncbiBaseUrl` | string | https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/ | no | Base URL for the NCBI Datasets API. |
| `enscode` | string |  | no |  |
| `cleanCache` | boolean |  | no | Option to clean the cache directory after running the pipeline. |
| `outDir` | string |  | yes | Output directory for the pipeline results. |
| `cacheDir` | string | /cache | no | Directory for caching files. |
| `files_latency` | integer | 60 | no | Latency in seconds for file operations. |
| `max_intron_size` | integer | 100000 | no | STAR option max_intron_size. |
| `genome_file` | string |  | no | Path for genome file. |
