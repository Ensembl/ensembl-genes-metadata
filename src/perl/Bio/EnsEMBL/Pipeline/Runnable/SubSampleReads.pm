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

Bio::EnsEMBL::Pipeline::Runnable::SubSampleReads

=head1 SYNOPSIS


=head1 DESCRIPTION

Module to generate random samples off the downloaded fastq files. These subsamples are then alligned against the genome
=cut

package Bio::EnsEMBL::Pipeline::Runnable::SubSampleReads;

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

  $self->param_required('subsample_size');
  $self->param_required('input_dir');
  $self->param_required('iid');
  $self->param_required('is_paired');
}


sub run {
  my ($self) = @_;
  my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -user   => $ENV{GBUSER},
        -dbname => $self->param('pipe_db'),
        -host   => $ENV{GBS1},
        -port   => $ENV{GBP1},
        -pass   => $ENV{GBPASS},
        -driver => $ENV{GBDRIVER},
    );
  my $sth = $dba->dbc->prepare("insert into read_count values(?,?,?,?)");
  my $fastq_file_name = $self->param('iid');
  my $readcount_file = catfile($self->param('input_dir'),'readcount.txt');
  my $sample_size = $self->param('subsample_size');
  (my $subsample_name = $fastq_file_name) =~ s/fastq.gz/fq/;
  $self->param('renamed_file', $subsample_name);
  my $subsample = catfile($self->param('input_dir'),$subsample_name);
  my $fastq_file = catfile($self->param('input_dir'),$fastq_file_name);
  say "Counting the reads in ".$fastq_file_name;

  my $read_length = 0;
  my $count = 0;

  open(PH, 'zcat '.$fastq_file.' |')
    or $self->throw('Could not open '.$fastq_file);
  while (my $line = <PH>) {
    if ($count%4 == 1){
      $read_length = length($line) if (length($line) > $read_length);
    }
    ++$count;
  }
  open(DH, ">>$readcount_file");
  unless($count) {
    $self->throw("Failed to count reads in the following file:\n".$fastq_file);
  }
  #store read counts in file
  say DH $fastq_file_name . "\t$count";
  say "Found ".$count." reads";
  $sth->bind_param(1,$fastq_file_name);
  $sth->bind_param(2,$count);
  $sth->bind_param(3,$self->param('taxon_id'));
  $sth->bind_param(4,$read_length);
  if ($sth->execute){
   say "Read count for $fastq_file_name stored successfully";
  }
  else{
    $self->throw("Failed to store read length for file $fastq_file_name");
  }
  #check if read count is more than subsample size
   if($count > $sample_size) {
    # #take a random sample of the reads. use same random seed to ensure same pairs are extracted
     `seqtk sample -s100 $fastq_file $sample_size > $subsample`;
   }
   else{
   #say "no need to sample, just set filename accordingly";
    `zcat $fastq_file > $subsample`;
    $subsample =~ s/fastq.gz/fq/;
   }
   #`rm $fastq_file`;
}

sub create_faidx {
  my ($self,$path,$fastq) = @_;

  if(-e $path.'/'.$fastq.'.fai') {
    $self->warning("You selected faidx, but the faidx exists, so will not try and create");
    return;
  }

  my $cmd = $self->param_required('samtools_path').' faidx '.$path.'/'.$fastq;
  my $faidx_res = system($cmd);
  if($faidx_res) {
    $self->throw("Failed to index file. Command:\n".$cmd);
  }
}

sub write_output {
  my ($self) = @_;
 
  my $fastq = $self->param('renamed_file');
  my $path = $self->param('input_dir');
  if($self->param('create_faidx')) {
  	  $self->create_faidx($path,$fastq);
  }
  $self->dataflow_output_id({iid=>$self->param('renamed_file'),species=>$self->param('species'),is_paired=>$self->param('is_paired')}, "2");
}

1;
