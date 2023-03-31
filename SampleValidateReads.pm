=head1 LICENSE

Copyright [2018-2019] EMBL-European Bioinformatics Institute

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

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadSequences

=head1 SYNOPSIS


=head1 DESCRIPTION

Base module to load data into Hive custom tables. Modules should override
create_row_data which will be called before loading data into the table.
the accession field will be return in an arrayref on branch 2

=cut

package SampleValidateReads;

use strict;
use warnings;
use feature 'say';
use File::Spec::Functions;
use File::Basename;
use POSIX;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

=head2 fetch_input

 Arg [1]    : None
 Description: Check that sequence_file and table_name exists. Load the module
              needed to parse the file
 Returntype :
 Exceptions :

=cut

sub fetch_input {
  my $self = shift;

  $self->param_required('validator');
  $self->param_required('input_dir');
  $self->param_required('iid');
  $self->param_required('is_paired');
    
  
}


sub run {
  my ($self) = @_;
  
  #command to validate read format
  my $cmd = $self->param('validator') . " --file " . $self->param('input_dir') . $self->param('iid');
  my $validation_result = system($cmd);
  
  #check if file meets fastq formatting standards
  if ($validation_result == 0){
    $self->param('flow', 1);
  }
  else{
    $self->param('flow', 0);
  }
  
}


sub write_output {
  my ($self) = @_;
  
}

1;
