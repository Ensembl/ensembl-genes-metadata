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
- Generates a summaries and check for GCA updates.
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
- `--release_date`: Assembly release date to filter by (e.g., 2019-01-01).
- `--output_dir`: Directory to save the CSV output files.


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

Upon execution, the script generates the following CSV files in the specified output directory (if no directory is specified it saves the to the working directory):
##### GB_metadata_bioproject_id
- `{BioProject}_filtered_assemblies.csv` - Filtered assembly data.
- `{BioProject}_filtered_assemblies_summary_statistics.csv` - Summary statistics of selected metrics.
- `{BioProject}_filtered_assemblies_info_result.csv` - Additional assembly information.
- `{BioProject}_filtered_assemblies_gca_list.csv` - List of GCA IDs.


##### Fetch_annotations
- `assembly_level_summary.csv` - Filtered assembly data binned by assembly level.
- `assemblies.csv` - Filtered assembly data.
- `assembly_annotation_status.csv` -Contains the status of assembly annotations, detailing the current state of each assembly.
- `annotation_method_summary.csv` - Contains a summary of annotation methods used in the analysis.
- `annotations.csv` - Contains filtered annotation data based on the specified criteria.
- `yearly_summary.csv` - Number of assemblies and annotations by year.
- `annotation_GCA_update_status.csv` - Contains the latest update status for GCAs associated with annotations, indicating if the latest version of a GCA is annotated.


### GUI in the browser
This currently only works with the GB_metadata_reporting module.
Start the GUI from the command-line:
```sh
streamlit run src/python/ensembl/genes/metadata/app.py
```

![image](https://github.com/user-attachments/assets/c3aa162d-a616-432a-a0ed-7236ef072c8c)
















