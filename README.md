# ensembl-genes-metadata
This repo contains modules and scripts for maintaning the Ensembl genebuild metadatabase
Currently repo contains a Perl and Python standalone module for adding assembly records to the database
To add a new assembly to the database, you have to ensure that you have the required perl5lIb for running the genebuild annotation pipeline in your path
Also, you need to have biopython and the listed python packages (see requirements.txt) installed either locally or centrally and set in your path
Lastly, make a copy of the registry.ini file and fill in the required fields. If running this module as a member of the Metazoa team, set the value of import_type to any of the listed values as reqired ('import_community','import_flybase','import_genbank','import_refseq','import_veupathdb','import_wormbase')

Note: The registry.ini file is the configuration file that contains the genbank accessions for custom registration and the database connection settings. 
      This module can be run be each genebuilder not necessarily the virtual user

To run the Python module, type the following command
$ENSCODE/ensembl-hive/scripts/standaloneJob.pl $ENSCODE/ensembl-genes-metadata/Register_Asm.py -config_file registry.ini -language python3
