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

package Bio::EnsEMBL::Pipeline::PipeConfig::AssemblyRegistryConf;
use strict;
use warnings;
use feature 'say';
use File::Spec::Functions;
use POSIX qw/strftime/;
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use parent ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');
use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Hive::Version 2.4;

sub default_options {
  my ($self) = @_;
  return {
    # inherit other stuff from the base class
    %{ $self->SUPER::default_options() }, 
    user => $ENV{'USER_R'},
    password => $ENV{'GBPASS'},
    user_w => $ENV{'GBUSER'},
    #hash to hold pipeline db settings
    'pipeline_db' => {
       -dbname => 'registry_assembly_pipe',
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
    
    'output_path' => $self->param('output_path') || $ENV{meta_database_dir},
      
  }
}
 
sub pipeline_analyses {
	
  my ($self) = @_;
  return [
    {
      -logic_name => 'sync_meta_database_with_production',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
         cmd => 'python '.$ENV{ENS_GENES_META} . '/src/python/ensembl/genes/metadata/sync_meta_database_production_db.py',
      },
      -rc_name => 'default',
      -flow_into => { 1 => ['check_for_updates_to_meta_database'],},  
      -input_ids => [{}],
    },
    {
      -logic_name => 'check_for_updates_to_meta_database',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
         cmd => 'python ' . $ENV{ENS_GENES_META} . '/src/python/ensembl/genes/metadata/check_for_new_assemblies.py --reg_path '.$self->o('output_path'),
      },
      -rc_name    => 'default',
      -flow_into => {
        1 => ['refseq_accession_assembly_name_update'],
        },
    },
    {
      -logic_name => 'refseq_accession_assembly_name_update',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
         cmd => 'python ' . $ENV{ENS_GENES_META} . '/src/python/ensembl/genes/metadata/check_and_update_refseq.py',
      },
      -rc_name => 'default',
      -flow_into => {
        1 => ['any_update'],
	 },
    }, 
    {			
      -logic_name => 'any_update',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name => 'default',
      -parameters => {
          cmd => 'if grep GCA_ "' .catfile($self->o('output_path'), 'assemblies_to_register.ini"').' > /dev/null; then exit 0; else exit 42;fi',
          return_codes_2_branches => {'42' => 2},
       },
      -flow_into => {1 => ['backup_db'],
         	     2 =>  ['no_update']
                 },
      },
    {
      -logic_name => 'no_update',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'python ' . $ENV{ENS_GENES_META} . '/src/python/ensembl/genes/metadata/no_update.py --msg #msg# ',
        msg => 'Nothing to update at this time ' .strftime('%Y-%m-%d',localtime),

      },
    },
    {
      -logic_name => 'backup_db',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'python ' . $ENV{ENS_GENES_META} . '/src/python/ensembl/genes/metadata/registry_backup.py --backup_file #output_file# --port '.$ENV{GBP1}. ' --server '. $ENV{GBS1} . ' --dbname '. $ENV{REG_DB} . ' --user '. $ENV{GBUSER_R},
        output_file => $self->o('output_path') . '/registry_db_bak/' . 'registry_bkup_'.strftime('%Y-%m-%d',localtime), . '.sql'

      },
      -flow_into => { 1 => ['register_assembly'],},
    },
    {
      -logic_name => 'register_assembly',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'python ' . $ENV{ENS_GENES_META} . '/src/python/ensembl/genes/metadata/register_assembly.py --config_file #config_file#',
        'config_file' => $self->o('output_path') . 'assemblies_to_register.ini',      
      },
      -rc_name => 'default',
      -flow_into => { 1 => ['sync_db'],},
     },
     {		
       -logic_name => 'sync_db',
       -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
       -parameters => {
        cmd => 'python ' . $ENV{ENS_GENES_META} . '/src/python/ensembl/genes/metadata/sync_db.py --backup_file #output_file# --port '.$ENV{GBP1}. ' --server '. $ENV{GBS1} . ' --dbname '. $ENV{REG_DB} . ' --user '. $ENV{GBUSER} . ' --p '. $ENV{GBPASS},
        output_file => $self->o('output_path') . 'registry_dump.sql',
      
      },
       -rc_name => 'default',
     },	
  ];
}
 

sub resource_classes {
	my ($self) = @_;
	return {
	  'default' => { LSF => $self->lsf_resource_builder('production', 3000)},
	};
}
  
sub pipeline_wide_parameters {
	my ($self) = @_;
	return {
  	};
}

sub pipeline_create_commands {
	my ($self) = @_;
	return [
	       'mkdir -p '.$self->param('output_path').'/registry_db_bak',
    # inheriting database and hive tables' creation
		@{$self->SUPER::pipeline_create_commands},
	];
}
1;
