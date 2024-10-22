# Genebuild assembly metadata pipeline

The pipelines connects to NCBI API to retrieve the latest Eucaryotic assemblies and helps to manage the correspoding metadata.
The data is stored in a MySQL assembly metadata database located in GB1 server. The pipeline use as main source the datasets NCBI API, download the information, collect and process the relevant information using multiple python scripts. 

Extra information is obtained from the taxonomy NCBI API and dtol API.  

The repository provides the configuration profile, which specify the input variables (mandatory and default) and two execution profiles: slurm and lsf

This pipeline has been tested in the new modular environment and with the latest Nextflow version (version 24.04.3) available in Codon 


# Requirements 

## Ensembl dependencies 


It is recommended that all the repositories are cloned into the same folder. 

| Repository name | branch | url |
|-------------------|-------|----|
| ensembl-genes | main | https://github.com/Ensembl/ensembl-genes.git |
|-------------------|-------|----|

# Running options

## Mandatory

- `--output_dir`: path to the directory where to store the results of the pipeline
- `--enscode`: environmental enscode variable   
- `-profile`: specify which execution profile use (slurm or lsf) 

## Default variables
- `--taxon`: taxonomy ID. Default value is 2759
- `--ncbi_url`: current URL to access 
- `--ncbi_params`: path to NCBI API params (file in data folder)
- `--db_table_conf`: path to db table configuration (file in data folder)
- `--metadata_params`: path to assembly metadata db credentials (file in data folder)
- `--registry_params` : path to assembly registry db credentials (file in data folder)

# Using the pipeline

**Getting help**
```bash
nextflow run ${ENSCODE}/ensembl-genes-metadata/pipeline/assembly_pipeline.nf --help
```

**Executing pipeline**
```bash
nextflow -C ${ENSCODE}/ensembl-genes-metadata/conf/assembly_pipeline.conf \
run ${ENSCODE}/ensembl-genes-metadata/pipeline/assembly_pipeline.nf \
--output_dir /path/to/output/dir/ \
--enscode ${ENSCODE} -profile slurm
```

# Standalone modules










