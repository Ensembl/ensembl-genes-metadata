=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package Bio::EnsEMBL::Pipeline::PipeConfig::TranscriptomicRegistryConf;

use strict;
use warnings;
use File::Spec::Functions;
use Bio::EnsEMBL::Hive::Version 2.5;
use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_analysis_settings);
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');

sub default_options {
  my ($self) = @_;
  return {
    # inherit other stuff from the base class
    %{ $self->SUPER::default_options() },
    user => $ENV{'GBUSER_R'},
    password => $ENV{'GBPASS'},
    user_w => $ENV{'GBUSER'},
    'pipe_db_name' => 'transcriptomic_assessment_pipeline',
    #hash to hold pipeline db settings
    'pipeline_db' => {
                       -dbname => $self->o('pipe_db_name'),
                       -host   => $ENV{GBS1},
                       -port   => $ENV{GBP1},
                       -user   => $self->o('user_w'),
                       -pass   => $self->o('password'),
                       -driver => 'mysql',
    },
    #hash to hold registry settings
    'output_db' => {
                       -dbname => $ENV{REG_DB},
                       -host   => $ENV{GBS1},
                       -port   => $ENV{GBP1},
                       -user   => $self->o('user_w'),
                       -pass   => $self->o('password'),
                       -driver => 'mysql',
    },
######################################################
#
# Variable settings- You change these!!!
#
######################################################
########################
# Misc setup info
########################

    'dbowner'                   => '' || $ENV{EHIVE_USER} || $ENV{USER},
    'user_r'                    => 'ensro', # read only db user
    'password_r'                => $ENV{GBPASS},
    'registry_server'           => $ENV{GBS1}, # host for general output dbs
    'pipe_db_port'              => $ENV{GBP1}, # port for pipeline host
    'registry_port'             => $self->o('pipe_db_port'), # port for general output db host
    'output_path'               => $ENV{meta_database_dir}, # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses
    'genome_path'               => catfile($self->o('output_path'), 'genomes'),
    'assembly_registry_host'    => $self->o('registry_server'),
    'assembly_registry_port'    => $self->o('registry_port'),
     samtools_path              => catfile($self->o('binary_base'), 'samtools'), #You may need to specify the full path to the samtools binary
     star_path                  => catfile($self->o('binary_base'), 'STAR'),
     fastqc                     => '/hps/software/users/ensembl/genebuild/FastQC/fastqc',
     validator_prog             => '/hps/software/users/ensembl/genebuild/fastQValidator/bin/fastQValidator',
    'rnaseq_dir'                => catdir($self->o('output_path'), 'genomes'),
    'rnaseq_ftp_base'           => 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/',
    'get_assembly_script_path'  => $ENV{ENSEMBL_GENES_META} . '/src/perl/scripts/get_assembly.pl',
    'tissue_dict_path'          => '/nfs/production/flicek/ensembl/genebuild/meta_database/',
    'minimap2_path'             => catfile($self->o('binary_base'), 'minimap2'),
    'species_list'                => catdir($self->o('output_path'), 'species.csv'),
    'species_list_nr'                => catdir($self->o('output_path'), 'nr_species.csv'),
    
    # What Read group tag would you like to group your samples
    # by? Default = ID
    read_group_tag => 'SM',
    read_id_tag    => 'ID',

    use_threads => 3,
    star_threads    => 12,    

    'summary_file_delimiter' => '\t',            # Use this option to change the delimiter for your summary data file
    'short_csv_table'      => 'csv_short_read',
    'long_csv_table'      => 'csv_long_read',
    'read_length_table'      => 'read_length',
    'read_count_table'      => 'read_count',
    'alignment_stats'      => 'alignment_stats',
    'filename_tag' => 'filename',
    #'rnaseq_search_file'        => '',
    ####################################################################
    # This is just an example based on the file snippet shown below.  It
    # will vary depending on how your data looks.
    ####################################################################
    short_read_columns => ['SM','ID','is_paired','taxon_id','filename','DS','url','md5','species','source_id'],
    file_columns => ['SM', 'species', 'accession' ],
    minimap_columns => ['iid', 'species', 'accession' ],
    assembly_columns => ['chain', 'version', 'ass_name', 'sp_name'],
    db_columns => ['sp_id', 'ass_id', 'chain', 'version', 'clade', 'contigN50', 'ass_lev', 'ass_name', 'sp_name', 'genus'],
    long_read_columns => ['SM','filename','taxon_id','DS','url','md5', 'species','source_id'],
    
    # Regular expression to allow FastQ files to be correctly paired,
    # for example: file_1.fastq and file_2.fastq could be paired using
    # the expression "\S+_(\d)\.\S+".  Need to identify the read number
    # in brackets; the name the read number (1, 2) and the
    # extension.
    pairing_regex => '\S+_(\d)\.\S+',
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # No option below this mark should be modified
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #
    ########################################################
    # URLs for retrieving the INSDC contigs and RefSeq files
    ########################################################
    'ncbi_base_ftp'           => 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all',
    'insdc_base_ftp'          => $self->o('ncbi_base_ftp').'/#expr(substr(#assembly_accession#, 0, 3))expr#/#expr(substr(#assembly_accession#, 4, 3))expr#/#expr(substr(#assembly_accession#, 7, 3))expr#/#expr(substr(#assembly_accession#, 10, 3))expr#/#assembly_accession#_#assembly_name#',
    'assembly_ftp_path'       => $self->o('insdc_base_ftp'),
    
    ########################
    # Extra db settings
    ########################

    'num_tokens' => 10,
    mysql_dump_options => '--max_allowed_packet=400MB',
    
    'registry_db' => {
      -host   => $self->o('registry_server'),
      -port   => $self->o('registry_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -dbname => $ENV{REG_DB},
      -driver => $self->o('hive_driver'),
    },
  };
}
sub pipeline_create_commands {
    my ($self) = @_;

    my $sr_tables;
    my $lr_tables;
    my %small_columns = (
        paired => 1,
        read_length => 1,
        is_13plus => 1,
        is_mate_1 => 1,
        );
    # We need to store the values of the csv file to easily process it. It will be used at different stages
    foreach my $key (@{$self->default_options->{'short_read_columns'}}) {
        if (exists $small_columns{$key}) {
            $sr_tables .= $key.' SMALLINT UNSIGNED NOT NULL,'
        }
        elsif ($key eq 'DS' or $key eq 'url') {
            $sr_tables .= $key.' VARCHAR(255) NOT NULL,'
        }
        elsif ($key eq 'SM' or $key eq 'CN') {
            $sr_tables .= $key.' VARCHAR(300) NOT NULL,'
        }
        elsif ($key eq 'taxon_id' or $key eq 'source_id') {
            $sr_tables .= $key.' INT(11) NOT NULL,'
        }
        else {
            $sr_tables .= $key.' VARCHAR(50) NOT NULL,'
        }
    }
    $sr_tables .= ' KEY(SM), KEY(ID)';

    foreach my $key (@{$self->default_options->{'long_read_columns'}}) {
        if (exists $small_columns{$key}) {
            $lr_tables .= $key.' SMALLINT UNSIGNED NOT NULL,'
        }
        elsif ($key eq 'DS' or $key eq 'url') {
            $lr_tables .= $key.' VARCHAR(255) NOT NULL,'
        }
        elsif ($key eq 'SM' or $key eq 'CN') {
            $lr_tables .= $key.' VARCHAR(300) NOT NULL,'
        }
        elsif ($key eq 'taxon_id' or $key eq 'source_id') {
            $lr_tables .= $key.' INT(11) NOT NULL,'
        }
        else {
            $lr_tables .= $key.' VARCHAR(50) NOT NULL,'
        }
    }
    $lr_tables .= ' KEY(SM), KEY(filename)';
    return [
      'mkdir -p '.$self->o('rnaseq_dir'),
       # inheriting database and hive tables' creation
        @{$self->SUPER::pipeline_create_commands},
        $self->db_cmd( 'CREATE TABLE ' . $self->o('short_csv_table') . " ($sr_tables)" ),
        $self->db_cmd( 'CREATE TABLE ' . $self->o('long_csv_table') . " ($lr_tables)" ),
        $self->db_cmd( 'CREATE TABLE ' . $self->o('read_count_table') . ' (' .
        'fastq varchar(50) NOT NULL,' .
        'read_count int(50) NOT NULL,' .
        'species_id int(11) NOT NULL,' .
        'read_length int(11) NOT NULL,' .
        'PRIMARY KEY (fastq,species_id))' ),
        $self->db_cmd( 'CREATE TABLE ' . $self->o('alignment_stats') . ' (' .
        'accession varchar(50) NOT NULL,' .
        'ID varchar(50) NOT NULL,' .
        'perc_mapped int(5) NOT NULL,' .
        'perc_paired int(5) NOT NULL,' .
        'PRIMARY KEY (accession,ID))' ),
    ];
}
sub pipeline_wide_parameters {
  my ($self) = @_;
  return {
  }
}
sub pipeline_analyses {
   my ($self) = @_;
    return [
    	{ 
			-logic_name => 'fetch_assembly',
			-module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
		   	-parameters => {
			  cmd => "perl ".$self->o('get_assembly_script_path') . " -file " . $self->o('species_list'),
			},
			-input_ids => [{}],
			 -rc_name => 'default',
			-flow_into => {
			  1 => {'create_csv_download_jobs' => {'inputfile' => $self->o('species_list_nr')}},
			},
			
		 },
		 {
			-logic_name => 'create_csv_download_jobs',
			-module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
			-parameters => {
			  column_names => $self->o('db_columns'),
			  delimiter => $self->o('summary_file_delimiter'),
			},
			-rc_name => 'default',
			-flow_into => { 
				         1      => ['backup_original_csv'],
				        '2->A' => ['download_short_read_csv','download_long_read_csv'],
					'A->1' => {'create_species_alignment_job' => {'inputfile' => $self->o('species_list')}},
			}, 
		},
		{
			-logic_name => 'backup_original_csv',
			-module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
			-parameters => {
				cmd => 'mv '.'#inputfile#'.' '.'#inputfile#'.'_orig_bak',
			},
                         -rc_name => 'default',
        
                },


		 { 
			-logic_name => 'download_short_read_csv',
			-module     => 'Bio::EnsEMBL::Pipeline::Runnable::FetchTranscriptomicData',
			-rc_name => '1GB',
			-parameters => {
                         species_list => $self->o('species_list'),
			 output_dir  => $self->o('output_path'),
			 tissue_dict  => $self->o('tissue_dict_path') . 'rnaseq_tissues_list.txt',
                         pipe_db => $self->o('pipe_db_name'),
			},
			-rc_name => 'default',
			-flow_into  => {
                         1 => WHEN ('#flow# == 1' => {'create_assembly_download_jobs' => { 'species' => '#sp_name#'}},
                                    '#flow# == 1' => {'create_sr_download_jobs' => {'csv_file' => '#csvfile#', 'species' => '#sp_name#'}},
                              ELSE {'download_genus_short_read_csv' => {'sp_name' => '#sp_name#', 'taxonomy' => '#sp_id#', 'sp_id' => '#genus#', 'genus' => 1, 'ass_name' => '#ass_name#', 'chain' => '#chain#', 'version' => '#version#'}}),
                         
                       },
                        
			 },
			 {
			 	 -logic_name => 'download_long_read_csv',
				 -module     => 'Bio::EnsEMBL::Pipeline::Runnable::FetchTranscriptomicData',
			 	 -rc_name => '1GB',
			 	 -parameters => {
			 	 output_dir  => $self->o('output_path'),
                                 species_list => $self->o('species_list'),
			 	 tissue_dict  => $self->o('tissue_dict_path') . 'rnaseq_tissues_list.txt',
			 	 read_type => 'isoseq',
                                 pipe_db => $self->o('pipe_db_name'),
			 	 },

			 	 -flow_into  => {
                         1 => WHEN ('#flow# == 1' => {'create_lr_download_jobs' => {'csv_file' => '#csvfile#', 'species' => '#sp_name#'}},
                              ELSE {'download_genus_long_read_csv' => {'sp_name' => '#sp_name#', 'taxonomy' => '#sp_id#', 'sp_id' => '#genus#', 'genus' => 1, 'ass_name' => '#ass_name#', 'chain' => '#chain#', 'version' => '#version#'}}),
                         
                       },
			 },
			 { 
				-logic_name => 'download_genus_short_read_csv',
				-module     => 'Bio::EnsEMBL::Pipeline::Runnable::FetchTranscriptomicData',
				-rc_name => '1GB',
				-parameters => {
                                 species_list => $self->o('species_list'),
				 output_dir  => $self->o('output_path'),
				 tissue_dict  => $self->o('tissue_dict_path') . 'rnaseq_tissues_list.txt',
				},
				-rc_name => 'default',
				-flow_into  => {
                         1 => WHEN ('#flow# == 1' => {'create_assembly_download_jobs' => { 'species' => '#sp_name#'}},
                                    '#flow# == 1' => {'create_sr_download_jobs' => {'csv_file' => '#csvfile#', 'species' => '#sp_name#'}},
                              ELSE ['no_rnaseq_data']),
                         
                       },
                        
			 },
			 {
			 	 -logic_name => 'download_genus_long_read_csv',
				 -module     => 'Bio::EnsEMBL::Pipeline::Runnable::FetchTranscriptomicData',
			 	 -rc_name => '1GB',
			 	 -parameters => {
			 	 	 output_dir  => $self->o('output_path'),
			 	 	 tissue_dict  => $self->o('tissue_dict_path') . 'rnaseq_tissues_list.txt',
			 	 	 read_type => 'isoseq',
                                         species_list => $self->o('species_list'),
			 	 },
	
			 	 -flow_into  => {
					 1 => WHEN ('#flow# == 1' => {'create_lr_download_jobs' => {'csv_file' => '#csvfile#', 'species' => '#sp_name#'}},
						  ELSE ['no_isoseq_data']),
							 
					 },
			 },
                         {
                                 -logic_name => 'create_assembly_download_jobs',
                                 -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
                                 -parameters => {
                                   inputcmd => 'grep -w #species# ' . $self->o('species_list') . ' | cut -f 3,4,8,9',
                                   column_names => $self->o('assembly_columns'),
                                   delimiter => '\t',
                                   
                                 },
                                 -rc_name => '5GB',
                                 -flow_into => {
                                         2 => ['download_assembly_data'],
                                 },
                         },

			 {
			 	 -logic_name => 'download_assembly_data',
			 	 -module     => 'Bio::EnsEMBL::Pipeline::Runnable::AssemblyLoading',	
			 	 -parameters => {
			 	 	 assembly_name => '#ass_name#',
			 	 	 species_name => '#sp_name#',
			 	 	 'assembly_accession' => '#chain#' . '.' . '#version#',
			 	 	 output_dir => $self->o('output_path'),
			 	 	 'full_ftp_path' => $self->o('insdc_base_ftp'),
			 	 	 'genome_path' => $self->o('genome_path'),
			 	 	 minimap2_path => $self->o('minimap2_path'),
                                         star_path => $self->o('star_path'),
			 	 },
			 	 -rc_name => '10GB',
			 	 -flow_into => {
			 	 	 -1 => ['download_assembly_data_20GB'],
			 	 },
			 }, 
			 {
			 -logic_name => 'download_assembly_data_20GB',
			 -module     => 'Bio::EnsEMBL::Pipeline::Runnable::AssemblyLoading',
			 -parameters => {
			 	assembly_name => '#ass_name#',
			 	species_name => '#sp_name#',
			 	'assembly_accession' => '#chain#' . '.' . '#version#',
			    output_dir => $self->o('output_path'),
			    'full_ftp_path' => $self->o('insdc_base_ftp'),
			    'genome_path' => $self->o('genome_path'),
			    minimap2_path => $self->o('minimap2_path'),
                            star_path => $self->o('star_path'),
			},
			 -rc_name => '20GB',
			 -flow_into => {
			 	 -1 => ['download_assembly_data_40GB'],
			 }
		 }, 
		 {
			 -logic_name => 'download_assembly_data_40GB',
			  -module     => 'Bio::EnsEMBL::Pipeline::Runnable::AssemblyLoading',
			 -parameters => {
			 	assembly_name => '#ass_name#',
			 	species_name => '#sp_name#',
			 	'assembly_accession' => '#chain#' . '.' . '#version#',
			    output_dir => $self->o('output_path'),
			    'full_ftp_path' => $self->o('insdc_base_ftp'),
			    'genome_path' => $self->o('genome_path'),
			    minimap2_path => $self->o('minimap2_path'),
                            star_path => $self->o('star_path'),
			},
			 -rc_name => '15GB',
		 }, 
		 	{
			 	 -logic_name => 'no_rnaseq_data',
				  -module     => 'Bio::EnsEMBL::Pipeline::Runnable::RNASeqReport',
			 	 -rc_name => 'default',
			 }, 
			 {
			 	 -logic_name => 'no_isoseq_data',
				 -module     => 'Bio::EnsEMBL::Pipeline::Runnable::RNASeqReport',
			 	 -rc_name => 'default',
			 },
			 {
			 	 -logic_name => 'create_sr_download_jobs',
			 	 -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
			 	 -parameters => {
			 	 	 inputfile => '#csv_file#',
			 	 	 species => '#species#',
			 	 	 column_names => $self->o('short_read_columns'),
			 	 	 delimiter => '\t',
			 	 },
			 	 -rc_name => 'default',
			 	 -flow_into => {
                                   2 => ['download_sr_fastq'],
			 	 },
			 },
			 {  
			-logic_name => 'download_sr_fastq',
			-module     => 'Bio::EnsEMBL::Pipeline::Runnable::FetchReads',
			-parameters =>{
			  ftp_base_url => $self->o('rnaseq_ftp_base'),
			  input_dir => $self->o('output_path'),
                          pipe_db => $self->o('pipe_db_name'),
			  
			},
			 -rc_name => '1GB',
			 -flow_into => {
			 	 2 => {'subsample_fastq' => {'iid' => '#filename#', 'taxon_id' => '#taxon_id#', 'species' => '#species#', 'is_paired' => '#is_paired#','source_id' => '#source_id#'}},
			 	 
			 }
		 },
		
		{
			-logic_name => 'create_lr_download_jobs',
			-module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
			-parameters => {
			  inputfile => '#csv_file#',
			  species => '#species#',
			  column_names => $self->o('long_read_columns'),
			  delimiter => '\t',
			},
			 -rc_name => 'default',
			-flow_into => {
                          2 => ['download_lr_isoseq'],
			},
		},
			
	
		 {  
			-logic_name => 'download_lr_isoseq',
			-module     => 'Bio::EnsEMBL::Pipeline::Runnable::FetchReads',
			-parameters =>{
			  ftp_base_url => $self->o('rnaseq_ftp_base'),
			  input_dir => $self->o('output_path'),
			  samtools_path => $self->o('samtools_path'),
			  decompress => 1,
			  create_faidx => 1,
                          pipe_db => $self->o('pipe_db_name'),
			},
			 -rc_name => '1GB',
			 -flow_into => {
			 	 2 => {'subsample_fastq' => {'iid' => '#filename#', 'species' => '#species#', 'taxon_id' => '#taxon_id#', 'is_paired' => '#is_paired#', 'source_id' => '#source_id#'}},
			 	 
			 }
		},
		{  
			-logic_name => 'subsample_fastq',
			-module     => 'Bio::EnsEMBL::Pipeline::Runnable::SubSampleReads',
			-parameters =>{
			  input_dir => $self->o('output_path').'#species#'. "/fastq/",
			  subsample_size => 50000,
                          pipe_db => $self->o('pipe_db_name'),
                          species => '#species#',
			  
			  
			},
			 -rc_name => '3GB',
			 -flow_into => {
			 	 -1 => ['subsample_fastq_6GB'],
			 	 2 => {'validate_fastq' => {'iid' => '#renamed_file#', 'species' => '#species#', 'is_paired' => '#is_paired#', 'taxon_id' => '#taxon_id#', 'source_id' => '#source_id#',}},
			 }
		},
		
		{  
			-logic_name => 'subsample_fastq_6GB',
			-module     => 'Bio::EnsEMBL::Pipeline::Runnable::SubSampleReads',
			-parameters =>{
			  input_dir => $self->o('output_path').'#species#'. "/fastq/",
			  subsample_size => 50000,
                          pipe_db => $self->o('pipe_db_name'),
                          species => '#species#',
			},
			 -rc_name => '6GB',
			 -flow_into => {
			 	 2 => {'validate_fastq' => {'iid' => '#renamed_file#', 'species' => '#species#', 'is_paired' => '#is_paired#', 'taxon_id' => '#taxon_id#',  'source_id' => '#source_id#',}},
			 }
		},
		
		{  
			-logic_name => 'validate_fastq',
			-module     => 'Bio::EnsEMBL::Pipeline::Runnable::SampleValidateReads',
			-parameters =>{
			  input_dir => $self->o('output_path').'#species#'. "/fastq/",
			  validator => $self->o('validator_prog'),	  
			},
			 -rc_name => '3GB',
			 -flow_into  => {
                         1 => WHEN ('#flow# == 1' => ['quality_check_6GB'],
                              ELSE ['failed_validation']),
                         
                       },
		},
		
		{
			 -logic_name => 'failed_validation',
			 -module => 'Bio::EnsEMBL::Pipeline::Runnable::RNASeqReport',	
			 -rc_name => 'default',
		 },
		 
		{  
			-logic_name => 'quality_check_6GB',
			-module     => 'Bio::EnsEMBL::Pipeline::Runnable::FastQCheck',
			-parameters =>{
			  input_dir => $self->o('output_path').'#species#'. "/fastq/",
			  quality_checker => $self->o('fastqc'). " -t 6",
                          pipe_db => $self->o('pipe_db_name'),
			},
			 -rc_name => '6GB',
			 -flow_into  => {
                        
                         -3 => ['quality_check_9GB'],
                         1 => WHEN ('#flow# == 1' => ['failed_quality_check'],
                             ),
                       },
		},
		
		{  
			-logic_name => 'quality_check_9GB',
			-module     => 'Bio::EnsEMBL::Pipeline::Runnable::FastQCheck',
			-parameters =>{
			  input_dir => $self->o('output_path').'#species#'. "/fastq/",
			  quality_checker => $self->o('fastqc'). " -t 9",
                          pipe_db => $self->o('pipe_db_name'),
			},
			 -rc_name => '9GB',
			 -flow_into  => {
                       
                         1 => WHEN ('#flow# == 1' => ['failed_quality_check'],
                             ),
                       },
		},
		
		{
			 -logic_name => 'failed_quality_check',
			 -module => 'Bio::EnsEMBL::Pipeline::Runnable::RNASeqReport',	
			 -rc_name => 'default',
		 },

                {
                        -logic_name => 'create_species_alignment_job',
                        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
                        -parameters => {
                          inputfile => $self->o('species_list'),
                          column_names => $self->o('db_columns'),
                          delimiter => '\t',
                        },
                         -rc_name => 'default',
                        -flow_into => {
                          2 => {'fan_sr_data' =>  {'accession' => '#chain#' . '.' . '#version#', 'species' => '#sp_name#', 'sp_id' => '#sp_id#', 'ass_id' => '#ass_id#'},'fan_lr_data' =>  {'accession' => '#chain#' . '.' . '#version#', 'species' => '#sp_name#', 'sp_id' => '#sp_id#', 'ass_id' => '#ass_id#'}},
                        },
                },

		{
			-logic_name => 'fan_sr_data',
			-module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
			-parameters => {
                                pipe_db_conn => $self->o('pipeline_db'),
				cmd => 'if [[ $(#db_query#) -gt 0 ]]; then exit 0; else exit 42;fi',
				return_codes_2_branches => {'42' => 2},
                                db_query => 'mysql -h #expr(#pipe_db_conn#->{-host})expr# -P #expr(#pipe_db_conn#->{-port})expr# -u #expr(#pipe_db_conn#->{-user})expr# -p#expr(#pipe_db_conn#->{-pass})expr# #expr(#pipe_db_conn#->{-dbname})expr# -NB -e "SELECT COUNT(*) FROM csv_short_read WHERE taxon_id = #sp_id#"',
			},
        	-rc_name    => 'default',
                -flow_into => {
                            '1->A' => ['create_star_jobs'],
                             'A->1' => ['classify_rnaseq_data'],

                          },
        },
        
        {
			-logic_name => 'fan_lr_data',
			-module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
                -parameters => {
                                pipe_db_conn => $self->o('pipeline_db'),
                                cmd => 'if [[ $(#db_query#) -gt 0 ]]; then exit 0; else exit 42;fi',
                                return_codes_2_branches => {'42' => 2},
                                db_query => 'mysql -h #expr(#pipe_db_conn#->{-host})expr# -P #expr(#pipe_db_conn#->{-port})expr# -u #expr(#pipe_db_conn#->{-user})expr# -p#expr(#pipe_db_conn#->{-pass})expr# #expr(#pipe_db_conn#->{-dbname})expr# -NB -e "SELECT COUNT(*) FROM csv_long_read WHERE taxon_id = #sp_id#"',
                        },
        	-rc_name    => 'default',
                -flow_into => {
                            '1->A' => ['generate_minimap2_jobs'],
                             'A->1' =>['classify_isoseq_data'],

                          },
        },
        
	          
                  {  
                    -logic_name => 'create_star_jobs',
                    -module     => 'Bio::EnsEMBL::Pipeline::Runnable::CreateStarJobs',
                    -parameters => {   
                      input_dir         => $self->o('output_path'),
                      sample_column    => $self->o('read_group_tag'),
                      sample_id_column => $self->o('read_id_tag'),
                      filename_column  => $self->o('filename_tag'),
                      pipe_db => $self->o('pipe_db_name'),
                      target_batch_size => 1,
                      csvfile_table    => $self->o('short_csv_table'),
                      column_names     => $self->o('file_columns'),
                    },
                    -rc_name   => 'default',
                    -flow_into => {
                      2 => ['star'],
                    },
                  },

                  {
                    -logic_name => 'star',
                    -module     => 'Bio::EnsEMBL::Pipeline::Runnable::StarAlign',
                    -parameters => {
                      disconnect_jobs    => 1,
                      input_dir          => $self->o('output_path').'#species#'. "/fastq/",
                      output_dir         => $self->o('output_path').'#species#/#accession#'.'/alignment',
                      short_read_aligner => $self->o('star_path'),
                      genome_dir         => catfile($self->o('output_path'),'#species#/#accession#/'),
                      num_threads        => $self->o('star_threads'),
                      species            => '#species#',
                      pipe_db => $self->o('pipe_db_name'),
                    },
                    -flow_into => {
                       -1 => ['star_85GB'],
                     },
                    -rc_name => '65GB_star',
                  },

                  {
                    -logic_name => 'star_85GB',
                    -module     => 'Bio::EnsEMBL::Pipeline::Runnable::StarAlign',
                    -parameters => {
                      disconnect_jobs    => 1,
                      input_dir          => $self->o('output_path').'#species#'. "/fastq/",
                      output_dir => $self->o('output_path').'#species#/#accession#'.'/alignment',
                      short_read_aligner => $self->o('star_path'),
                      genome_dir         => catfile($self->o('output_path'),'#species#/#accession#/'),
                      num_threads        => $self->o('star_threads'),
                      species            => '#species#',
                      pipe_db => $self->o('pipe_db_name'),
                    },
                    -rc_name => '85GB_star',
                  },

		  
		  {
		-logic_name => 'generate_minimap2_jobs',
                  -module     => 'Bio::EnsEMBL::Pipeline::Runnable::CreateMinimapJobs',
                  -parameters => {
                     pipe_db => $self->o('pipe_db_name'),
                     column_names     => $self->o('minimap_columns'),
		  },
		    -rc_name => 'default',
		    -flow_into => {
                      2 => ['minimap2'],
		    },
          },


          {
          	  -logic_name => 'minimap2',
          	  -module     => 'Bio::EnsEMBL::Pipeline::Runnable::Minimap',
          	  -parameters => {
                         genome_file => catfile($self->o('output_path').'#species#/#accession#/'),
                         minimap2_path => $self->o('minimap2_path'),
                         input_dir => $self->o('output_path').'#species#'. "/fastq/",
                         output_dir => $self->o('output_path').'#species#/#accession#'.'/alignment',
                         pipe_db => $self->o('pipe_db_name'),
                         samtools => $self->o('samtools_path'),
                       	},
              -rc_name => '15GB',
              -flow_into => {
                           -1 => ['minimap2_himem'],
              },
          },
          
          {
          	  -logic_name => 'minimap2_himem',
          	  -module     => 'Bio::EnsEMBL::Pipeline::Runnable::Minimap',
          	  -parameters => {
                         genome_file => catfile($self->o('output_path').'#species#/#accession#/'),
                         minimap2_path => $self->o('minimap2_path'),
                         input_dir => $self->o('output_path').'#species#'. "/fastq/",
                         output_dir => $self->o('output_path').'#species#/#accession#'.'/alignment',
                         samtools => $self->o('samtools_path'),
                         pipe_db => $self->o('pipe_db_name'),
                       	},
              -rc_name => '30GB',
              
          },

		  { 
		    -logic_name => 'classify_rnaseq_data',
		    -module     => 'Bio::EnsEMBL::Pipeline::Runnable::ClassifyTranscriptomicData',
		   	-parameters => {
			  pipe_db => $self->o('pipe_db_name'), 
			  },
			 -flow_into => {
		       1 => ['delete_output_files'],
		     },
			  -rc_name => '5GB',
			 
			
		  },
		  { 
		    -logic_name => 'classify_isoseq_data',
		    -module     => 'Bio::EnsEMBL::Pipeline::Runnable::ClassifyIsoseqData',
		   	-parameters => {
			  pipe_db => $self->o('pipe_db_name'),
			  },
			  -flow_into => {
		       1 => ['delete_output_files'],
		     },
		    -rc_name => '5GB',
		  },
		  {
		  	  -logic_name => 'delete_output_files',
		  	  -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
		  	  -parameters => {
                      cmd => 'rm -r '.$self->o('output_path').'#species# ',
                     },
              -rc_name    => 'default',
               -max_retry_count => 1,
               -wait_for => ['classify_isoseq_data','classify_rnaseq_data'],
		  
          },

	];
}
sub resource_classes {
  my $self = shift;

  return {
  	'default' => { LSF => $self->lsf_resource_builder('production', 900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '1GB' => { LSF => $self->lsf_resource_builder('production', 1000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'registry_server'}], [$self->default_options->{'num_tokens'}])},
    '3GB' => { LSF => $self->lsf_resource_builder('production', 3000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '5GB' => { LSF => $self->lsf_resource_builder('production', 5000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '6GB' => { LSF => $self->lsf_resource_builder('production', 6000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '9GB' => { LSF => $self->lsf_resource_builder('production', 6000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '10GB' => { LSF => $self->lsf_resource_builder('production', 10000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '15GB' => { LSF => $self->lsf_resource_builder('production', 15000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '20GB' => { LSF => $self->lsf_resource_builder('production', 15000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '30GB' => { LSF => $self->lsf_resource_builder('production', 30000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'registry_server'}], [$self->default_options->{'num_tokens'}])},
    '50GB' => { LSF => $self->lsf_resource_builder('production', 50000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'registry_server'}], [$self->default_options->{'num_tokens'}])},
    '40GB' => { LSF => $self->lsf_resource_builder('production', 40000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'registry_server'}], [$self->default_options->{'num_tokens'}])},
    '85GB_star'    => { LSF => $self->lsf_resource_builder( 'production', 85000, undef, undef, ( $self->default_options->{'star_threads'} + 1 ) ) },
    '65GB_star'    => { LSF => $self->lsf_resource_builder( 'production', 65000, undef, undef, ( $self->default_options->{'star_threads'} + 1 ) ) },  
  }
}
sub hive_capacity_classes {
  my $self = shift;

  return {
           'hc_low'    => 200,
           'hc_medium' => 500,
           'hc_high'   => 1000,
         };
}

1;
