# Genebuild metadata handling system

This Python module allows users to query the MySQL databases to retrieve genome assembly and/or annotation information based on specified BioProject IDs/release date  and filter the results based on various metrics. The script fetches relevant genome assembly and annotation data, applies user-defined thresholds, and outputs the results in tabular format, both on-screen and as CSV files.

## Current modules
- GB_metadata_bioproject_id
- GB_metadata_time
- Fetch_annotations

## Features
### GB_metadata_bioproject_id
- Fetches genome assembly data linked to specified BioProject IDs.
- Supports filtering based on metrics such as GC content, sequence length, Contig N50, and others.
- Retrieves taxonomic information and assigns internal clades based on predefined settings.
- Determines whether an assembly is a reference genome using NCBI's Assembly database.
- Generates summary statistics (mean, min, max) for selected metrics.
- Outputs results in a structured table and saves them as CSV files.

### GB_metadata_time
- Same as GB_metadata_bioproject_id but everything is linked to release date instead of Bioproject ID.

### Fetch_annotations
- Module based on GB_metadata_time.
- Retrieves information on assemblies as well as annotations by genebuild.
- Generates a summary by date.
- Outputs results in a structured table and saves them as CSV files.


## Prerequisites
- Python 3.x
- MySQL database with the necessary tables and data.
- An external JSON file (`db_config.json`) containing MySQL database credentials.
- `clade_settings.json` file for internal clade assignments.

## Installation
1. Clone this repository:
   ```sh
   git clone https://github.com/Ensembl/ensembl-genes-metadata.git
   git checkout dev/gb_metadata_handling
   cd ensembl-genes-metadata
   ```
Optional step: create a virtual environment 

2. Install required Python dependencies:
   ```sh
   pip install -r requirements.txt
   ```

3. Fill out database login credentials in conf/db_config.json and confs/prod_dbs_conf.json



## Usage

### Command line
Run the script with command-line arguments specifying BioProject IDs and optional filters:
```sh
python src/python/ensembl/genes/metadata/GB_metadata_reporting.py --bioproject_id PRJEB40665 PRJEB61747 --asm_level "Complete genome" --output_dir ./results
```

#### Command-Line Arguments
##### GB_metadata_bioproject_id
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
- `--release_date`: Assembly release date to filter by (e.g., 2019-01-01).
- `--output_dir`: Directory to save the CSV output files.

##### GB_metadata_time
- `--release_date`(required): Assembly release date to filter by (e.g., 2019-01-01).
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

##### Fetch_annotations
- `--release_date`(required): Assembly release date to filter by (e.g., 2019-01-01).


#### Output

Upon execution, the script generates the following CSV files in the specified output directory (if no directory is specified it saves the to the working directory):
##### GB_metadata_bioproject_id
- `{BioProject}_filtered_assemblies.csv` - Filtered assembly data.
- `{BioProject}_filtered_assemblies_summary_statistics.csv` - Summary statistics of selected metrics.
- `{BioProject}_filtered_assemblies_info_result.csv` - Additional assembly information.
- `{BioProject}_filtered_assemblies_gca_list.csv` - List of GCA IDs.

##### GB_metadata_time
- `filtered_assemblies_after_{release_date}.csv` - Filtered assembly data.
- `filtered_assemblies_summary_statistics_after_{release_date}.csv` - Summary statistics of selected metrics.
- `filtered_assemblies_info_result_after_{release_date}.csv` - Additional assembly information.
- `filtered_assemblies_gca_list_after_{release_date}.csv`

##### Fetch_annotations
- `assembly_level_summary.csv` - Filtered assembly data binned by assembly level.
- `contig_n50_summary.csv` - Filtered assembly data binned by N50) (bin size = 5000).
- `avg_contig_n50_by_level.csv` - Filtered assembly N50 data binned by assembly level.
- `annotation.csv` - Filtered annotation data binned by annotation method.
- `yearly_summary.csv` - Number of assemblies and annotations by year.


### GUI in the browser
This currently only works with the GB_metadata_biproject_id module.
Start the GUI from the command-line:
```sh
streamlit run src/python/ensembl/genes/metadata/app.py
```

![image](https://github.com/user-attachments/assets/c3aa162d-a616-432a-a0ed-7236ef072c8c)
















