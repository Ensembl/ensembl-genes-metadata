# Genebuild metadata handling system

This Python module allows users to query the MySQL databases to retrieve genome assembly and/or annotation information based on specified BioProject IDs/release date  and filter the results based on various metrics. The script fetches relevant genome assembly and annotation data, applies user-defined thresholds, and outputs the results in tabular format, both on-screen and as CSV files.

## Current modules
- GB_metadata_reporting
- Fetch_annotations

## Features
### GB_metadata_reporting
- Fetches genome assembly data linked to specified BioProject IDs.
- Supports filtering based on metrics such as GC content, sequence length, Contig N50, and others.
- Supports filtering based on taxon ID.
- Retrieves taxonomic information and assigns internal clades based on predefined settings.
- Determines whether an assembly is a reference genome using NCBI's Assembly database.
- Generates summary statistics for selected metrics.
- Outputs results in structured tables and saves them as CSV files.

### Fetch_annotations
- Module based on GB_metadata_reporting.
- Retrieves information on assemblies as well as annotations by genebuild.
- Supports filtering based on taxon ID.
- Generates summaries and check for GCA updates.
- Looks for assemblies that are not annotated.
- Outputs results in structured tables and saves them as CSV files.


## Prerequisites
- Python 3.x
- MySQL database with the necessary tables and data.
- An external JSON files (`conf/db_config.json` and `conf/prod_dbs_conf.json`) containing MySQL database credentials.
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
   Note: list includes GUI app requirements as well. No need to install them id you are only using the command line options.

3. Fill out database login credentials in conf/db_config.json and confs/prod_dbs_conf.json



## Usage

### Command line
Run the script with command-line arguments specifying BioProject IDs and optional filters:
```sh
python src/python/ensembl/genes/metadata/GB_metadata_reporting.py --bioproject_id PRJEB40665 PRJEB61747 --asm_level "Complete genome" --output_dir ./results
```

#### Command-Line Arguments
- `--bioproject_id`: One or more BioProject IDs to query.
- `--gc_percent`: Minimum GC content percentage.
- `--total_sequence_length`: Minimum total sequence length (bp).
- `--contig_n50`: Minimum Contig N50 value.
- `--number_of_contigs`: Minimum number of contigs.
- `--number_of_scaffolds`: Minimum number of scaffolds.
- `--scaffold_n50`: Minimum Scaffold N50 value.
- `--genome_coverage`: Minimum genome coverage.
- `--asm_level`: Assembly level(s) to filter by (e.g., Contig, Scaffold, Chromosome, Complete genome).
- `--asm_type`: Assembly type(s) to filter by (e.g., haploid, diploid, etc.).
- `--taxon_id`: NCBI Taxon ID to filter by (e.g., 40674 for Mammalia)
- `--release_date`: Assembly release date to filter by (e.g., 2019-01-01).
- `--output_dir`: Directory to save the CSV output files.
- `--reference`: Check if GCA is a reference genome (1 for yes, 0 for no). Default 0. Note: NCBI API times out if checking a large number of assemblies.


##### Fetch_annotations
## Parameters

The following parameters can be specified when running the script:

- `--release_date`: Release date in YYYY-MM-DD format (default: `2019-01-01`). Will be applied to both assemblies and annotations.
- `--output_dir`:Directory to save output files.
- `--taxon_id`:  Taxon ID for filtering (e.g., `40674` for Mammalia). Will be applied to both assemblies and annotations.
- `--asm_level`: Assembly level options. Acceptable values include: `Contig`, `Scaffold`, `Chromosome`,`Complete genome`
- `--asm_type`: Assembly type options. Acceptable values include: `haploid`, `alternate-pseudohaplotype`, `unresolved-diploid`, `haploid-with-alt-loci`, `diploid`
- `--contig_n50`:Assembly contig N50 threshold.
- `--gc_percent`: Assembly GC percent threshold.
- `--total_sequence_length`: Assembly sequence length threshold in base pairs (bp).
- `--number_of_contigs`: Assembly contig threshold.
- `--number_of_scaffolds`: Assembly scaffolds threshold.
- `--scaffold_n50`: 
  Scaffold N50 value.
- `--genome_coverage`: 
  Assembly genome coverage threshold .
- `--bioproject_id`: One or more BioProject IDs.

#### Output

Upon execution, the script generates the following CSV files in the specified output directory (or the working directory if none is provided)::
##### GB_metadata_reporting
Each output file follows a structured naming convention based on the BioProject ID and release date provided by the user.

- `{BioProject}_{release_date}_filtered_assemblies.csv` – Contains the filtered assembly data based on the applied metrics and thresholds.
- `{BioProject}_{release_date}_filtered_assemblies_summary_statistics.csv` – Summary statistics for the selected metrics.
- `{BioProject}_{release_date}_filtered_assemblies_info_result.csv` – Additional assembly metadata.
- `{BioProject}_{release_date}_filtered_assemblies_gca_list.csv` – A list of GCAs extracted from the filtered results.

##### Fetch_annotations
Each output file is named using a structured format, where {date_str} corresponds to the execution date, ensuring version control and easy tracking of results over time.
- `annotation_method_summary_{date_str}.csv` – Summary of annotation methods used.
- `yearly_summary_{date_str}.csv` – Yearly summary of annotations, assemblies and Ensembl releases.
- `annotations_{date_str}.csv` – Filtered annotation dataset based on specified criteria.
- `assemblies_{date_str}.csv` – Filtered assembly dataset based on specified criteria.
- `assembly_annotation_status_{date_str}.csv` – Status of assembly annotations, including metadata.
- `assembly_level_summary_{date_str}.csv` – Summary of assembly levels across the dataset.
- `annotation_GCA_update_status_{date_str}.csv` – Status of GCAs of already released annotations.


### GUI in the browser
This currently only works with the GB_metadata_reporting module.
Start the GUI from the command-line:
```sh
streamlit run src/python/ensembl/genes/metadata/app.py
```

![image](https://github.com/user-attachments/assets/c3aa162d-a616-432a-a0ed-7236ef072c8c)
















