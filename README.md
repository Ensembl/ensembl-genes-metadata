# Genebuild metadata handling system

This Python module allows users to query a MySQL database to retrieve genome assembly information based on specified BioProject IDs and filter the results based on various metrics. The script fetches relevant genome assembly data, applies user-defined thresholds, and outputs the results in tabular format, both on-screen and as CSV files.

## Features
- Fetches genome assembly data linked to specified BioProject IDs.
- Supports filtering based on metrics such as GC content, sequence length, Contig N50, and others.
- Retrieves taxonomic information and assigns internal clades based on predefined settings.
- Determines whether an assembly is a reference genome using NCBI's Assembly database.
- Generates summary statistics (mean, min, max) for selected metrics.
- Outputs results in a structured table and saves them as CSV files.

## Prerequisites
- Python 3.x
- MySQL database with the necessary tables and data.
- An external JSON file (`db_config.json`) containing MySQL database credentials.
- `clade_settings.json` file for internal clade assignments.

## Installation
1. Clone this repository:
   ```sh
   git clone https://github.com/Ensembl/ensembl-genes-metadata/tree/dev/gb_metadata_handling
   cd filtered-assemblies-query-tool
   ```
2. Install required Python dependencies:
   ```sh
   pip install pandas pymysql requests argparse streamlit
   ```

## Usage

### Command line
Run the script with command-line arguments specifying BioProject IDs and optional filters:
```sh
python script.py --bioproject_id PRJEB40665 PRJEB61747 --gc_percent 40.0 --total_sequence_length 1000000 --asm_level Scaffold --output_dir ./results
```

#### Command-Line Arguments
- `--bioproject_id` (required): One or more BioProject IDs to query.
- `--gc_percent`: Minimum GC content percentage.
- `--total_sequence_length`: Minimum total sequence length (bp).
- `--contig_n50`: Minimum Contig N50 value.
- `--number_of_contigs`: Minimum number of contigs.
- `--number_of_scaffolds`: Minimum number of scaffolds.
- `--scaffold_n50`: Minimum Scaffold N50 value.
- `--genome_coverage`: Minimum genome coverage.
- `--asm_level`: Assembly level(s) to filter by (e.g., Contig, Scaffold, Chromosome, Complete genome).
- `--asm_type`: Assembly type(s) to filter by (e.g., haploid, diploid, etc.).
- `--output_dir`: Directory to save the CSV output files.

#### Output
Upon execution, the script generates the following CSV files in the specified output directory (if no directory is specified it saves the to the working directory):
- `{BioProject}_filtered_assemblies.csv` - Filtered assembly data.
- `{BioProject}_filtered_assemblies_summary_statistics.csv` - Summary statistics of selected metrics.
- `{BioProject}_filtered_assemblies_info_result.csv` - Additional assembly information.
- `{BioProject}_filtered_assemblies_gca_list.csv` - List of GCA IDs.

### GUI in the browser
Start the GUI from the command-line:
```sh
streamlit run app.py
```















