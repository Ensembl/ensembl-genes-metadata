# Ensembl-genes-metadata

The EnsEMBL Genebuild Meta Database System

This repo contains modules and scripts for maintaning the Ensembl genebuild metadatabase
Currently, this repo contains both Perl and Python modules and/or scripts for adding assembly records to the database
To add a new assembly to the database, you have to ensure that you have the required PERL5LIB for running the genebuild annotation pipeline in your path

## Requirements
In addition to the listed python packages (see requirements.txt), you will need to have biopython installed either locally or centrally and set in your path, 

### Perl EnsEMBL repositories you need to have

We recommend that you clone all the repositories into one directory
| Repository name | branch | URL|
|-----------------|--------|----|
| ensembl | default | https://github.com/Ensembl/ensembl.git |
| ensembl-hive | default | https://github.com/Ensembl/ensembl-hive.git |
| ensembl-compara | release/98 | https://github.com/Ensembl/ensembl-compara.git |
| ensembl-production | default | https://github.com/Ensembl/ensembl-production.git |
| ensembl-taxonomy | default | https://github.com/Ensembl/ensembl-taxonomy.git |
| ensembl-orm | default | https://github.com/Ensembl/ensembl-orm.git |
| ensembl-killlist | default | https://github.com/Ensembl/ensembl-killlist.git |
| ensembl-datacheck | default | https://github.com/Ensembl/ensembl-datacheck.git |
| ensembl-metadata | default | https://github.com/Ensembl/ensembl-metadata.git |
| ensembl-io | default | https://github.com/Ensembl/ensembl-io.git |

For each of these repository, you will need to install their dependencies using the cpanfile provided in their Git repositories

You can use the [Ensembl git commands](https://github.com/Ensembl/ensembl-git-tools) and run the following command to clone the repositories
```
git ensembl --clone genebuild
```

### Python EnsEMBL repositories you need to have

| Repository name | branch | URL|
|-----------------|--------|----|
| ensembl-genes | default | https://github.com/Ensembl/ensembl-genes.git |

### Python virtual environment

You will need to create the virtual environment:
- `genebuild-metadb` using the requirements.txt file; it needs to be activated for the pipeline to run

### Shell environment

If you are not part of the Ensembl Genebuild team, you will need to set some shell environment variables to avoid having to provide the information to the configuration files. We will assume you are using your home directory
| Variable | Value | Hive configuration parameter | Description |
|----------|-------|------------------------------|-------------|
| ENSCODE | $HOME | -enscode\_root\_dir | Directory path where you cloned all the Perl repositories |
| ENSEMBL\_SOFTWARE\_HOME | $HOME | -software\_base\_path | Directory where pyenv, plenv and linuxbrew are installed |
| LINUXBREW\_HOME | $HOME/.linuxbrew | -linuxbrew\_home\_path | Base directory for your Linuxbrew installation |
| PYTHONPATH | $HOME/ensembl-genes/ensembl\_genes:$HOME/ensembl-hive/wrappers/python3/ | | It needs to be set until the package can be installed properly |

### MySQL

We currently use MySQL databases to store our data. To avoid having to do many changes to the configuration files we recommend having one read-only user and one read-write user. It is also better to use different servers for keeping fail-safe copies of the database.

## Running the EnsEMBL Genebuild Meta Database System (Assembly registry pipeline)

There is a main configuration file, `Bio::EnsEMBL::Pipeline::PipeConfig::AssemblyRegistrationConf`, which will generate a set of analyses to:
- sync genebuild entries between the meta database and the production portal
- check for new eukaryote genomes
- check and update assembly names and refseq accessions
- backup meta database when new new eukaryote genomes become available
- register new genomes
- copy and restore updated database across fail-safe servers
The whole system is explained in more details below

### Initialising the pipeline

You will need to activate the genebuild virtual environment
```
pyenv activate genebuild-metadb
```

#### Filling the main configuration automatically

If you are operating within an environment prepared for Ensembl with the assembly registry you can use the `$ENSCODE/ensembl-genes-metadata/src/perl/Bio/EnsEMBL/Pipeline/PipeConfig/AssemblyRegistryConf.pm`.

You would need to edit `$ENSCODE/ensembl-genes-metadata/config/registry.ini`

Note: The registry.ini file is the configuration file that contains the genbank accessions for custom registration and the database connection settings. 

Then you can run
```
perl $ENSCODE/ensembl-hive/scripts/init_pipeline.pl $ENSCODE/ensembl-genes-metadata/src/perl/Bio/EnsEMBL/Pipeline/PipeConfig/AssemblyRegistryConf.pm -output_path <writable location>
```

### Running the Assembly registry pipeline

To start the pipeline you need the URL to your pipeline database which will be provided when running the init\_Pipeline.pl script. If you initialised the pipeline automatically, you need to look at the command file created in your `working_dir` directory to retrieve the information.
```
export EHIVE_URL=mysql://readwrite_user:password@host:port/dbname
```

You can now start the pipeline with
```
perl $ENSCODE/ensembl-hive/scripts/beekeeper.pl -url $EHIVE_URL -loop
```

If you only want to run some analyses, you can run
```
perl $ENSCODE/ensembl-hive/scripts/beekeeper.pl -url $EHIVE_URL -loop -analyses_pattern 1..5
```

### Monitoring the pipeline

#### GuiHive

To follow the pipeline steps, it is better to use GuiHive, a graphical interface to ensembl-hive, which allows you to change parameters, debug your problems and much more https://github.com/Ensembl/guiHive

#### What to do when the main pipeline fails

You should first look at the job tab to know the reason of the failure
* Insufficient memory: you can either use a different resource or add a new one more suited to your needs
* Error in the code: I'm afraid you will need to do proper debugging

Once you are happy with your fix, you would need to reset the jobs with
```
perl $ENSCODE/ensembl-hive/scripts/beekeeper.pl -url $EHIVE_URL -reset_failed_jobs
```
and restart the pipeline

#### How can I debug a job

By default ensembl-hive redirect all output to `/dev/null` unless you used some logging parameters.

You will need to run the problematic job with runWorker. First you will need to retrieve the job id using GuiHive or the pipeline database. Then you can run
```
perl $ENSCODE/ensembl-hive/scripts/runWorker.pl -url $EHIVE_URL -debug 1 -job_id XX
```

Using a higher value for `-debug` is usually not useful as it is mostly seen as a boolean flag.

## The different parts of the EnsEMBL Assembly Registration System

### Sync meta database

#### What it does

It will reference the production portal to identify all live databases linked to the current Ensembl release. These databases are the queried to obtain meta information such as the genebuild method, genebuild completion date and genome accession. The retrieved meta data is checked against all pending genebuilds and when a match is found, the pending genebuild status is updated. 

#### Notifications

None

#### Caveats

When a database appears to be present both in Rapid and Main sites, the released server is set to Main

### Check for new assembly

#### What it does

It checks the public archives for any new eukaryotic genome since the last meta database update. If it finds one, an output file containing all new accessions and meta database configuration settings is created.

#### Notifications

None

#### Caveats

None

### Refseq and Assembly name check

#### What it does

This checks for any changes in the Refseq and assembly name values of existing assemblies.

#### Notifications

All updated assemblies are reported via Slack 

#### Caveats

None

### Assembly registration

#### What it does

All new assemblies obtained earlier are recorded in the database. 

#### Notifications

All new assemblies are reported via Slack 

#### Caveats

Depending on whether the new assembly comes from an already existing species or not, there may be need to create and assign a unique prefix and stable id space range to the assembly record.

### Meta database sync

#### What it does

When the meta database gets updated, copies are taken and restored across fail-safe servers 

#### Notifications

None

#### Caveats

None

### Running the Transcriptomic assessment pipeline

To start the pipeline you need the URL to your pipeline database which will be provided when running the init\_Pipeline.pl script. If you initialised the pipeline automatically, you need to look at the command file created in your `working_dir` directory to retrieve the information.
```
export EHIVE_URL=mysql://readwrite_user:password@host:port/dbname
```

You can now start the pipeline with
```
perl $ENSCODE/ensembl-hive/scripts/beekeeper.pl -url $EHIVE_URL -loop
```

If you only want to run some analyses, you can run
```
perl $ENSCODE/ensembl-hive/scripts/beekeeper.pl -url $EHIVE_URL -loop -analyses_pattern 1..5
```
## The different parts of the EnsEMBL Transcriptomic Assessment System

### Fetch candidate assemblies

#### What it does

It queries the meta database to retrieve all assemblies with contig_N50 > 100000 with total gap length less than 30% of the genome. 

#### Notifications

None

#### Caveats

When multiple assemblies exist for one species, it returns a list of one assembly per species. 

### Transcriptomic data check

#### What it does

It queries the ENA for any available transcriptomic data per species. 

#### Notifications

It returns a list of species with available data and those without

#### Caveats

When no data is found at species level, it retries at the genus level

### Fetch reads

#### What it does

It downloads both short read and long read data where available per species

#### Notifications

None

#### Caveats

None

### Prepare assembly for alignment

#### What it does

The genomes of all assemblies earlier identified are downloaded and indexed, ready for alignment

#### Notifications

None

#### Caveats

None

###  Subsampling of reads

#### What it does

All fastq files are subsampled randomly to a size of 50000 reads and the initial large files are erased.

#### Notifications

None

#### Caveats

None

### Read validation

#### What it does

The subsampled files are validated to ensure the files confrom to standard fastq requirements. 
Also, the per base sequence quality is tested using Fastqc

#### Notifications

None

#### Caveats

Reads not meeting the validation criteria or failing the per base sequence quality tests are discarded

### Read alignment

#### What it does

The reads that pass the prior validation steps are aligned against their corresponding genomes. Long reads are aligned using Minimap and short reads are aligned using Star.

#### Notifications

None

#### Caveats

None

### Read classification

#### What it does

Using arbitrary criteria such as percentage mapping quality, per base sequence quality, total read count per sample, the samples are classed as good, weak and unusable.
An assembly with five or more good samples gets a green status.
An assembly with more than one good or weak sample, with total read count greater than 100000000 gets an amber status.
An assembly not meeting either of the above gets a red status.

#### Notifications

None

#### Caveats

None

