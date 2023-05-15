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

package rnaseq_registry_codon_conf;

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
    user => 'ensro',
    password => 'ensembl',
    user_w => 'ensadmin',
    #hash to hold pipeline db settings
    'pipeline_db' => {
                       -dbname => 'transcriptomic_registry_model_org',
                       -host   => $ENV{GBS4},
                       -port   => $ENV{GBP4},
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
    'output_path'               => '/nfs/production/flicek/ensembl/genebuild/do1/transcriptomic_registry/', # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses
    'genome_path'               => '/hps/nobackup/flicek/ensembl/genebuild/transcriptomic_registry/genomes/',
    'assembly_registry_host'    => $self->o('registry_server'),
    'assembly_registry_port'    => $self->o('registry_port'),
     samtools_path              => catfile($self->o('binary_base'), 'samtools'), #You may need to specify the full path to the samtools binary
     star_path                  => catfile($self->o('binary_base'), 'STAR'),
     bwa_local                  => catfile($self->o('binary_base'), 'bwa'),
     fastqc                     => '/nfs/production/flicek/ensembl/genebuild/do1/FastQC/fastqc',
     validator_prog             => '/nfs/production/flicek/ensembl/genebuild/do1/fastQValidator/bin/fastQValidator',
    'rnaseq_dir'                => catdir($self->o('output_path'), 'rnaseq'),
    'input_dir'                 => catdir($self->o('rnaseq_dir'),'input'),
    'output_dir'                => catdir($self->o('rnaseq_dir'),'output'),
    'rnaseq_ftp_base'           => 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/',
    'get_assembly_script_path'  => '/homes/do1/assembly_registry_codon/get_candidate_assembly.pl',
    'merge_script_path'         => '/homes/do1/perl_scripts/merge_bam_files.pl',
    'tissue_dict_path'          => '/nfs/production/flicek/ensembl/genebuild/do1/',
    'minimap2_path'             => catfile($self->o('binary_base'), 'minimap2'),
    
    # What Read group tag would you like to group your samples
    # by? Default = ID
    read_group_tag => 'SM',
    read_id_tag    => 'ID',

    use_threads => 3,
    star_threads    => 12,    

    'summary_file_delimiter' => '\t',            # Use this option to change the delimiter for your summary data file
    'summary_csv_table'      => 'csv_data',
    'read_length_table'      => 'read_length',
    'filename_tag' => 'filename',
    #'rnaseq_search_file'        => '',
    ####################################################################
    # This is just an example based on the file snippet shown below.  It
    # will vary depending on how your data looks.
    ####################################################################
    short_read_columns => ['SM', 'ID', 'is_paired', 'taxon_id', 'filename', 'is_mate_1', 'read_length', 'is_13plus', 'CN', 'PL', 'DS', 'fastq_ftp', 'fastq_md5'],
    file_columns => [ 'SM', 'ID', 'is_paired', 'taxon_id', 'filename', 'is_mate_1', 'read_length', 'is_13plus', 'CN', 'PL', 'DS', 'url', 'md5sum' ],
    db_columns => ['sp_id', 'ass_id', 'chain', 'version', 'clade', 'contigN50', 'ass_lev', 'ass_name', 'sp_name', 'genus'],
    long_read_columns => ['SM','filename','is_paired','taxon_id','DS','url','md5sum'],
    
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

    my $tables;
    my %small_columns = (
        paired => 1,
        read_length => 1,
        is_13plus => 1,
        is_mate_1 => 1,
        );
    # We need to store the values of the csv file to easily process it. It will be used at different stages
    foreach my $key (@{$self->default_options->{'short_read_columns'}}) {
        if (exists $small_columns{$key}) {
            $tables .= $key.' SMALLINT UNSIGNED NOT NULL,'
        }
        elsif ($key eq 'DS' or $key eq 'url') {
            $tables .= $key.' VARCHAR(255) NOT NULL,'
        }
        elsif ($key eq 'SM' or $key eq 'CN') {
            $tables .= $key.' VARCHAR(138) NOT NULL,'
        }
        else {
            $tables .= $key.' VARCHAR(50) NOT NULL,'
        }
    }
    $tables .= ' KEY(SM), KEY(ID)';

    return [
      'mkdir -p '.$self->o('rnaseq_dir'),
       # inheriting database and hive tables' creation
        @{$self->SUPER::pipeline_create_commands},
        $self->db_cmd( 'CREATE TABLE ' . $self->o('csvfile_table') . " ($tables)" ),

        $self->db_cmd( 'CREATE TABLE ' . $self->o('read_length_table') . ' (' .
        'fastq varchar(50) NOT NULL,' .
        'read_length int(50) NOT NULL,' .
        'PRIMARY KEY (fastq))' ),
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
			-logic_name => 'get_candidate_assembly',
			-module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
		   	-parameters => {
			  cmd => "perl ".$self->o('get_assembly_script_path') . " -file " . $self->o('output_path') . "final_candidate_assembly.csv",
			},
			-input_ids => [{}],
			 -rc_name => 'default',
			-flow_into => {
			  1 => {'create_csv_download_jobs' => {'inputfile' => $self->o('output_path') .'final_candidate_assembly.csv'}},
			},
			
		 },
		 {
			-logic_name => 'create_csv_download_jobs',
			-module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
			-parameters => {
			  column_names => $self->o('db_columns'),
			  delimiter => '\t',
			},
			-rc_name => 'default',
			-flow_into => { 
				            1      => ['backup_original_csv'],
				            '2->A' => ['download_short_read_csv','download_long_read_csv'],
							'A->1' => {'create_sample_merge_jobs' => {'inputfile' => '#inputfile#'}},
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
			-module     => 'DownloadCsvENA',
			-rc_name => '1GB',
			-parameters => {
			 output_dir  => $self->o('output_path'),
			 tissue_dict  => $self->o('tissue_dict_path') . 'rnaseq_tissues_list.txt',
			},
			-rc_name => 'default',
			-flow_into  => {
                         1 => WHEN ('#flow# == 1' => ['download_assembly_data'],#{'assembly_name' => '#ass_name#', 'assembly_accession' => '#accession#', 'species_name' => '#sp_name#'},
                                    '#flow# == 1' => {'create_sr_download_jobs' => {'csv_file' => '#csvfile#', 'species' => '#sp_name#'}},
                              ELSE {'download_genus_short_read_csv' => {'sp_name' => '#sp_name#', 'sp_id' => '#genus#', 'genus' => 1, 'ass_name' => '#ass_name#', 'chain' => '#chain#', 'version' => '#version#'}}),
                         
                       },
                        
			 },
			 {
			 	 -logic_name => 'download_long_read_csv',
			 	 -module     => 'DownloadCsvENA',
			 	 -rc_name => '1GB',
			 	 -parameters => {
			 	 output_dir  => $self->o('output_path'),
			 	 tissue_dict  => $self->o('tissue_dict_path') . 'rnaseq_tissues_list.txt',
			 	 read_type => 'isoseq',
			 	 },

			 	 -flow_into  => {
                         1 => WHEN ('#flow# == 1' => ['download_assembly_data'],
                         	         '#flow# == 1' => {'create_lr_download_jobs' => {'csv_file' => '#csvfile#', 'species' => '#sp_name#'}},
                              ELSE {'download_genus_long_read_csv' => {'sp_name' => '#sp_name#', 'sp_id' => '#genus#', 'genus' => 1, 'ass_name' => '#ass_name#', 'chain' => '#chain#', 'version' => '#version#'}}),
                         
                       },
			 },
			 { 
				-logic_name => 'download_genus_short_read_csv',
				-module     => 'DownloadCsvENA',
				-rc_name => '1GB',
				-parameters => {
				 output_dir  => $self->o('output_path'),
				 tissue_dict  => $self->o('tissue_dict_path') . 'rnaseq_tissues_list.txt',
				},
				-rc_name => 'default',
				-flow_into  => {
                         1 => WHEN ('#flow# == 1' => ['download_assembly_data'],#{'assembly_name' => '#ass_name#', 'assembly_accession' => '#accession#', 'species_name' => '#sp_name#'},
                                    '#flow# == 1' => {'create_sr_download_jobs' => {'csv_file' => '#csvfile#', 'species' => '#sp_name#'}},
                              ELSE ['no_rnaseq_data']),
                         
                       },
                        
			 },
			 {
			 	 -logic_name => 'download_genus_long_read_csv',
			 	 -module     => 'DownloadCsvENA',
			 	 -rc_name => '1GB',
			 	 -parameters => {
			 	 	 output_dir  => $self->o('output_path'),
			 	 	 tissue_dict  => $self->o('tissue_dict_path') . 'rnaseq_tissues_list.txt',
			 	 	 read_type => 'isoseq',
			 	 },
	
			 	 -flow_into  => {
					 1 => WHEN ('#flow# == 1' => ['download_assembly_data'],
					 	        '#flow# == 1' => {'create_lr_download_jobs' => {'csv_file' => '#csvfile#', 'species' => '#sp_name#'}},
						  ELSE ['no_isoseq_data']),
							 
					 },
			 },
			 {
			 	 -logic_name => 'download_assembly_data',
			 	 -module => 'AssemblyLoading',	
			 	 -parameters => {
			 	 	 assembly_name => '#ass_name#',
			 	 	 species_name => '#sp_name#',
			 	 	 'assembly_accession' => '#chain#' . '.' . '#version#',
			 	 	 output_dir => $self->o('output_path'),
			 	 	 'full_ftp_path' => $self->o('insdc_base_ftp'),
			 	 	 'genome_path' => $self->o('genome_path'),
			 	 	 short_read_aligner => $self->o('bwa_local'),
			 	 	 minimap2_path => $self->o('minimap2_path'),
                                         star_path => $self->o('star_path'),
			 	 },
			 	 -rc_name => '5GB',
			 	 -flow_into => {
			 	 	 -1 => ['download_assembly_data_10GB'],
			 	 },
			 }, 
			 {
			 -logic_name => 'download_assembly_data_10GB',
			 -module => 'AssemblyLoading',	
			 -parameters => {
			 	assembly_name => '#ass_name#',
			 	species_name => '#sp_name#',
			 	'assembly_accession' => '#chain#' . '.' . '#version#',
			    output_dir => $self->o('output_path'),
			    'full_ftp_path' => $self->o('insdc_base_ftp'),
			    'genome_path' => $self->o('genome_path'),
			    short_read_aligner => $self->o('bwa_local'),
			    minimap2_path => $self->o('minimap2_path'),
                            star_path => $self->o('star_path'),
			},
			 -rc_name => '10GB',
			 -flow_into => {
			 	 -1 => ['download_assembly_data_15GB'],
			 }
		 }, 
		 {
			 -logic_name => 'download_assembly_data_15GB',
			 -module => 'AssemblyLoading',	
			 -parameters => {
			 	assembly_name => '#ass_name#',
			 	species_name => '#sp_name#',
			 	'assembly_accession' => '#chain#' . '.' . '#version#',
			    output_dir => $self->o('output_path'),
			    'full_ftp_path' => $self->o('insdc_base_ftp'),
			    'genome_path' => $self->o('genome_path'),
			    short_read_aligner => $self->o('bwa_local'),
			    minimap2_path => $self->o('minimap2_path'),
                            star_path => $self->o('star_path'),
			},
			 -rc_name => '15GB',
		 }, 
		 	{
			 	 -logic_name => 'no_rnaseq_data',
			 	 -module => 'RNASeqReport',	
			 	 -rc_name => 'default',
			 }, 
			 {
			 	 -logic_name => 'no_isoseq_data',
			 	 -module => 'RNASeqReport',	
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
			 	 2 => {'download_sr_fastq' => {'iid' => '#filename#', 'species' => '#species#', 'is_paired' => '#is_paired#', 'taxon_id' => '#taxon#', 'fastq_ftp' => '#fastq_ftp#', 'fastq_md5' => '#fastq_md5#'}},
			 	 },
			 },
			 {  
			-logic_name => 'download_sr_fastq',
			-module     => 'HiveDownloadRNASeqFastqs',
			-parameters =>{
			  ftp_base_url => $self->o('rnaseq_ftp_base'),
			  input_dir => $self->o('output_path').'#species#'. "/",
			  
			},
			 -rc_name => '1GB',
			 -flow_into => {
			 	 2 => {'subsample_fastq' => {'iid' => '#iid#', 'species' => '#species#', 'is_paired' => '#is_paired#'}},
			 	 
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
			  2 => {'download_lr_isoseq' => {'iid' => '#filename#', 'species' => '#species#', 'is_paired' => '#is_paired#', 'taxon_id' => '#taxon#', 'fastq_ftp' => '#fastq_ftp#', 'fastq_md5' => '#fastq_md5#'}},
			},
		},
			
	
		 {  
			-logic_name => 'download_lr_isoseq',
			-module     => 'HiveDownloadRNASeqFastqs',
			-parameters =>{
			  ftp_base_url => $self->o('rnaseq_ftp_base'),
			  input_dir => $self->o('output_path').'#species#'. "/",
			  samtools_path => $self->o('samtools_path'),
			  decompress => 1,
			  create_faidx => 1,
                          pipe_db_conn => $self->o('pipeline_db'),
			},
			 -rc_name => '1GB',
			 -flow_into => {
			 	 2 => {'subsample_fastq' => {'iid' => '#iid#', 'species' => '#species#', 'is_paired' => '#is_paired#'}},
			 	 
			 }
		},
		{  
			-logic_name => 'subsample_fastq',
			-module     => 'SubSampleReads',
			-parameters =>{
			  input_dir => $self->o('output_path').'#species#'. "/fastq/",
			  subsample_size => 50000,
			  
			  
			},
			 -rc_name => '3GB',
			 -flow_into => {
			 	 -1 => ['subsample_fastq_6GB'],
			 	 2 => {'validate_fastq' => {'iid' => '#renamed_file#', 'species' => '#species#', 'is_paired' => '#is_paired#'}},
			 }
		},
		
		{  
			-logic_name => 'subsample_fastq_6GB',
			-module     => 'SubSampleReads',
			-parameters =>{
			  input_dir => $self->o('output_path').'#species#'. "/fastq/",
			  subsample_size => 50000,
			},
			 -rc_name => '6GB',
			 -flow_into => {
			 	 2 => {'validate_fastq' => {'iid' => '#renamed_file#', 'species' => '#species#', 'is_paired' => '#is_paired#'}},
			 }
		},
		
		{  
			-logic_name => 'validate_fastq',
			-module     => 'SampleValidateReads',
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
			 -module => 'RNASeqReport',	
			 -rc_name => 'default',
		 },
		 
		{  
			-logic_name => 'quality_check_6GB',
			-module     => 'FastQCheck',
			-parameters =>{
			  input_dir => $self->o('output_path').'#species#'. "/fastq/",
			  quality_checker => $self->o('fastqc'). " -t 6",
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
			-module     => 'FastQCheck',
			-parameters =>{
			  input_dir => $self->o('output_path').'#species#'. "/fastq/",
			  quality_checker => $self->o('fastqc'). " -t 9",
			},
			 -rc_name => '9GB',
			 -flow_into  => {
                       
                         -3 => ['quality_check_9GB'],
                         1 => WHEN ('#flow# == 1' => ['failed_quality_check'],
                             ),
                       },
		},
		
		{
			 -logic_name => 'failed_quality_check',
			 -module => 'RNASeqReport',	
			 -rc_name => 'default',
		 },

			 {
			-logic_name => 'create_sample_merge_jobs',
			-module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
			-parameters => {
			  inputfile => '#inputfile#',
			  column_names => $self->o('db_columns'),
			  delimiter => '\t',
			},
			 -rc_name => 'default',
			-flow_into => {
			  2 => {'fan_sr_merge_csv' =>  {'accession' => '#chain#' . '.' . '#version#', 'species' => '#sp_name#', 'sp_id' => '#sp_id#', 'ass_id' => '#ass_id#'},'fan_lr_merge_csv' =>  {'accession' => '#chain#' . '.' . '#version#', 'species' => '#sp_name#', 'sp_id' => '#sp_id#', 'ass_id' => '#ass_id#'}},
			},
		},
		{
			-logic_name => 'fan_sr_merge_csv',
			-module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
			-parameters => {
				cmd => 'if compgen -G "'.catfile($self->o('output_path').'#species#'. "/", 'fastq/*_sr.csv"').' > /dev/null; then exit 0; else exit 42;fi',
				return_codes_2_branches => {'42' => 2},
			},
        	-rc_name    => 'default',
        	-flow_into  => { '1' => ['merge_sr_sample_csv'] },
        },
        
        {
			-logic_name => 'fan_lr_merge_csv',
			-module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
			-parameters => {
				cmd => 'if compgen -G "'.catfile($self->o('output_path').'#species#'. "/", 'fastq/*_lr.csv"').' > /dev/null; then exit 0; else exit 42;fi',
				return_codes_2_branches => {'42' => 2},
			},
        	-rc_name    => 'default',
        	-flow_into  => { '1' => ['merge_lr_sample_csv'] },
        },
        
		{ 
		    -logic_name => 'merge_sr_sample_csv',
		    -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
		   	-parameters => {
			  
			  cmd => 'cat '.catfile($self->o('output_path').'#species#'. "/", 'fastq/*_sr.csv').' > '.$self->o('output_path').'#species#'. "/".'bwa_file_list.txt',
			  },
			  -rc_name => 'default',
			  -flow_into => {
			    '1->A' => ['create_bwa_jobs','create_star_jobs'],
					'A->1' => {'classify_rnaseq_data' => {'species' => '#species#', 'sp_id' => '#sp_id#', 'ass_id' => '#ass_id#'}},
			 
			  },
			
		  },
		  
		  { 
		    -logic_name => 'merge_lr_sample_csv',
		    -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
		   	-parameters => {
			  
			  cmd => 'cat '.catfile($self->o('output_path').'#species#'. "/", 'fastq/*_lr.csv').' > '.$self->o('output_path').'#species#'. "/".'minimap_file_list.txt',
			  },
			  -rc_name => 'default',
			  -flow_into => {
			    '1->A' => ['generate_minimap2_jobs'],
					'A->1' => {'classify_isoseq_data' => {'species' => '#species#', 'sp_id' => '#sp_id#', 'ass_id' => '#ass_id#'}},
			 
			  },
			
		  },
	          
                  {  
                    -logic_name => 'create_star_jobs',
                    -module     => 'HiveCreateStarJobs',
                    -parameters => {   
                      input_dir         => $self->o('output_path'),
                      sample_column    => $self->o('read_group_tag'),
                      sample_id_column => $self->o('read_id_tag'),
                      filename_column  => $self->o('filename_tag'),
                      csvfile_table    => $self->o('summary_csv_table'),
                      column_names     => $self->o('file_columns'),
                    },
                    -rc_name   => 'default',
                    -flow_into => {
                      2 => ['star'],
                    },
                  },

                  {
                    -logic_name => 'star',
                    -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveStar',
                    -parameters => {
                      disconnect_jobs    => 1,
                      input_dir          => $self->o('output_path').'#species#'. "/fastq/",
                      output_dir => $self->o('output_path').'#species#'.'/alignment',
                      short_read_aligner => $self->o('star_path'),
                      genome_dir         => catfile($self->o('genome_path').'#species#'. "/",  '#accession#'),,
                      num_threads        => $self->o('star_threads'),
                    },
                    -flow_into => {
                       -1 => ['star_15GB'],
                     },
                    -rc_name => '5GB_star',
                  },

                  {
                    -logic_name => 'star_15GB',
                    -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveStar',
                    -parameters => {
                      disconnect_jobs    => 1,
                      input_dir => $self->o('output_path').'#species#'. "/fastq/",
                      output_dir => $self->o('output_path').'#species#'.'/alignment',
                      short_read_aligner => $self->o('star_path'),
                      genome_dir         => catfile($self->o('genome_path').'#species#'. "/",  '#accession#'),
                      num_threads        => $self->o('star_threads'),
                    },
                    -rc_name => '15GB_star',
                  },

		  {  
		    -logic_name => 'create_bwa_jobs',
		    -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
		    -parameters => {
		      inputfile => $self->o('output_path').'#species#'. "/".'bwa_file_list.txt',
		      column_names => ['mate1', 'mate2'],
		      delimiter => '\t',
		    },
		    -rc_name => 'default',
		    -flow_into => {
		      2 => {'bwa_paired_alignment' => {'species' => '#species#', 'read1' => '#mate1#', 'read2' => '#mate2#', 'accession' => '#accession#'}},
		    },
			
		   },
	
		   { 
		     
		     -logic_name => 'bwa_paired_alignment',
		     -module     => 'BWA2BAM',
		     -parameters => {
		       short_read_aligner => $self->o('bwa_local'),
		       input_dir => $self->o('output_path').'#species#'. "/fastq/",
		       genome_file => catfile($self->o('genome_path').'#species#'. "/",  '#accession#'),
		       output_dir => $self->o('output_path').'#species#'.'/alignment',
		       samtools => $self->o('samtools_path'),
		     },
		     -rc_name    => '5GB',
		     -flow_into => {
		       -1 => ['bwa_paired_alignment_15GB'],
		     },
		  },
		  
		  { 
		     
		     -logic_name => 'bwa_paired_alignment_15GB',
		     -module     => 'BWA2BAM',
		     -parameters => {
		       short_read_aligner => $self->o('bwa_local'),
		       input_dir => $self->o('output_path').'#species#'. "/fastq/",
		       genome_file => catfile($self->o('genome_path').'#species#'. "/",  '#accession#'),
		       output_dir => $self->o('output_path').'#species#'.'/alignment',
		       samtools => $self->o('samtools_path'),
		     },
		      -flow_into => {
		       -1 => ['bwa_paired_alignment_30GB'],
		     },
		     -rc_name    => '15GB',
			
		  },
		  
		  { 
		     
		     -logic_name => 'bwa_paired_alignment_30GB',
		     -module     => 'BWA2BAM',
		     -parameters => {
		       short_read_aligner => $self->o('bwa_local'),
		       input_dir => $self->o('output_path').'#species#'. "/fastq/",
		       genome_file => catfile($self->o('genome_path').'#species#'. "/",  '#accession#'),
		       output_dir => $self->o('output_path').'#species#'.'/alignment',
		       samtools => $self->o('samtools_path'),
		     },
		     -flow_into => {
		       -1 => ['bwa_paired_alignment_50GB'],
		     },
		     -rc_name    => '30GB',
			
		  },
		  
		  { 
		     
		     -logic_name => 'bwa_paired_alignment_50GB',
		     -module     => 'BWA2BAM',
		     -parameters => {
		       short_read_aligner => $self->o('bwa_local'),
		       input_dir => $self->o('output_path').'#species#'. "/fastq/",
		       genome_file => catfile($self->o('genome_path').'#species#'. "/",  '#accession#'),
		       output_dir => $self->o('output_path').'#species#'.'/alignment',
		       samtools => $self->o('samtools_path'),
		     },
		     -rc_name    => '50GB',
			
		  },
		  
		  {
		  	  -logic_name => 'generate_minimap2_jobs',
               -rc_name      => '2GB',
               -flow_into => {
                        2 => {'minimap2' => {'input_file' => '#fastq_file#','iid' => '#iid#'}},
                      },
             -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
             -parameters => {
		      inputfile => $self->o('output_path').'#species#'. "/".'minimap_file_list.txt',
		      column_names => ['mate1', 'mate2'],
		      delimiter => '\t',
		    },
		    -rc_name => 'default',
		    -flow_into => {
		      2 => {'minimap2' => {'species' => '#species#', 'read1' => '#mate1#', 'accession' => '#accession#'}},
		    },
          },


          {
          	  -logic_name => 'minimap2',
          	  -module     => 'Minimap2',
          	  -parameters => {
                         genome_file => catfile($self->o('genome_path').'#species#'. "/",  '#accession#'),
                         minimap2_path => $self->o('minimap2_path'),
                         input_dir => $self->o('output_path').'#species#'. "/fastq/",
                         output_dir => $self->o('output_path').'#species#'.'/alignment',
                         samtools => $self->o('samtools_path'),
                       	},
              -rc_name => '15GB',
              -flow_into => {
                       	   -1 => {'minimap2_himem' => {'species' => '#species#', 'read1' => '#mate1#', 'accession' => '#accession#'}},
              },
          },
          
          {
          	  -logic_name => 'minimap2_himem',
          	  -module     => 'Minimap2',
          	  -parameters => {
                         genome_file => catfile($self->o('genome_path').'#species#'. "/",  '#accession#'),
                         minimap2_path => $self->o('minimap2_path'),
                         input_dir => $self->o('output_path').'#species#'. "/fastq/",
                         output_dir => $self->o('output_path').'#species#'.'/alignment',
                         samtools => $self->o('samtools_path'),
                       	},
              -rc_name => '30GB',
              
          },

		  { 
		    -logic_name => 'classify_rnaseq_data',
		    -module     => 'ClassifyTranscriptomicData',
		   	-parameters => {
		   	  flagstats => $self->o('output_path').'#species#'.'/alignment/flagstats.txt',
		   	  csv_file => $self->o('output_path').'#species#'.'/'.'#species#_rnaseq'.'.csv',
		   	  read_length_file => $self->o('output_path').'#species#'.'/fastq/fastq_read_length.csv',
		   	  fastq_report_file => $self->o('output_path').'#species#'.'/fastq/fastqc_report.txt',
			  
			  },
			 -flow_into => {
		       1 => ['delete_output_files'],
		     },
			  -rc_name => '5GB',
			 
			
		  },
		  { 
		    -logic_name => 'classify_isoseq_data',
		    -module     => 'ClassifyTranscriptomicData',
		   	-parameters => {
		   	  flagstats => $self->o('output_path').'#species#'.'/alignment/flagstats.txt',
		   	  csv_file => $self->o('output_path').'#species#'.'/'.'#species#_isoseq'.'.csv',
		   	  read_length_file => $self->o('output_path').'#species#'.'/fastq/fastq_read_length.csv',
		   	  fastq_report_file => $self->o('output_path').'#species#'.'/fastq/fastqc_report.txt',
			  
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
                      cmd => 'rm -r '.$self->o('output_path').'#species# ' .$self->o('genome_path').'#species#',
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
    '30GB' => { LSF => $self->lsf_resource_builder('production', 30000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'registry_server'}], [$self->default_options->{'num_tokens'}])},
    '50GB' => { LSF => $self->lsf_resource_builder('production', 50000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'registry_server'}], [$self->default_options->{'num_tokens'}])},
    '40GB' => { LSF => $self->lsf_resource_builder('production', 40000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'registry_server'}], [$self->default_options->{'num_tokens'}])},
    '5GB_star'    => { LSF => $self->lsf_resource_builder( 'production', 5000, undef, undef, ( $self->default_options->{'star_threads'} + 1 ) ) },
    '15GB_star'    => { LSF => $self->lsf_resource_builder( 'production', 15000, undef, undef, ( $self->default_options->{'star_threads'} + 1 ) ) },  
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
