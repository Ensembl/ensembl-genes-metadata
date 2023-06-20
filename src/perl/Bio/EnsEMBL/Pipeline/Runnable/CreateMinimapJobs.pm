=head1 LICENSE
 
Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=cut

package Bio::EnsEMBL::Pipeline::Runnable::CreateMinimapJobs;

use strict;
use warnings;
use File::Spec::Functions;
use feature 'say';
use Bio::EnsEMBL::Pipeline::Runnable::TranscriptomicRegistryAdaptor;
use parent ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');

=head2 param_defaults

 Arg [1]    : None
 Description: Default parameters
                'compression_ratio' => 3,
                'target_batch_size' => 10000000000,
 Returntype : Hashref
 Exceptions : None

=cut

sub param_defaults {
  my $self = shift;

  return {
    %{$self->SUPER::param_defaults},
    'compression_ratio' => 3,
    'target_batch_size' => 10000000000,
  }
}

=head2 fetch_input

 Arg [1]    : None
 Description: Creates input id based on a custom table 'csv_long_read' in the hive database
              It will generate the parameters for Minimap based on the data for each file
              It stores the input ids in 'inputlist'
 Returntype : None
 Exceptions : None

=cut

sub fetch_input {
  my $self = shift;
  my $registry_adaptor = new Bio::EnsEMBL::Pipeline::Runnable::TranscriptomicRegistryAdaptor(
        -user   => $ENV{GBUSER},
        -dbname => $self->param('pipe_db'),
        -host   => $ENV{GBS1},
        -port   => $ENV{GBP1},
        -pass   => $ENV{GBPASS},
        -driver => 'mysql',
    );
   my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -user   => $ENV{GBUSER},
        -dbname => $ENV{REG_DB},
        -host   => $ENV{GBS1},
        -port   => $ENV{GBP1},
        -pass   => $ENV{GBPASS},
        -driver => $ENV{GBDRIVER},
   );
  #Retrieve only reads that are considered good for alignment to genome
  my $fastqc_entry = $dba->dbc->prepare("select * from fastqc where species_id = ? and per_base_seq_quality in (?,?)");
  my %fqrep;
  $fastqc_entry->bind_param(1,$self->param('sp_id'));
  $fastqc_entry->bind_param(2,'PASS');
  $fastqc_entry->bind_param(3,'WARN');
  if ($fastqc_entry->execute){
    while (my @entry = $fastqc_entry->fetchrow_array()){
      $fqrep{$entry[1]} = lc($entry[2]);
    }
  }  
  my $output_ids = [];
  my $csv_long_read = $registry_adaptor->dbc->prepare("select filename from csv_long_read where taxon_id = ?");
  $csv_long_read->bind_param(1,$self->param('sp_id'));
  if ($csv_long_read->execute){
    while (my @reads = $csv_long_read->fetchrow_array()){
      if (!exists $fqrep{$reads[0]}){
        next;
      }
      push (@{$output_ids}, [$reads[0],$self->param('species'),$self->param('accession')]);
    }
  }
  if (!(@$output_ids)){
      $self->complete_early("No long read data has passed the test for alignment");
    }
  $self->param('inputlist',$output_ids);
}

1;
