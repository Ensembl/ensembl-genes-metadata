# ensembl-genes-metadata
This repo contains modules and scripts for maintaning the Ensembl genebuild metadatabase
Currently repo contains a standalone module for adding assembly records to the database
To add a new assembly to the database, you have to ensure that you have the required perl5lIb for running the genebuild annotation pipeline in your path
Also, you need to have biopython installed and set in your path
Lastly, make a copy of the registry.ini file and fill in the required fields

To run the module, type the following command
$ENSCODE/ensembl-hive/scripts/standaloneJob.pl $ENSCODE/ensembl-genes-metadata/AssemblyRegistrySubmission.pm -config_file registry.ini
