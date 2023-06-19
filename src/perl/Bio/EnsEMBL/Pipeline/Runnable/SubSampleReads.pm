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
  #$self->param_required('rnaseq_summary_file');
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
  my $fastq_read_length = $self->param('input_dir') . "/fastq_read_count.csv";
  open(DH, ">>$fastq_read_length");
  my $fastq_file_name = $self->param('iid');
  my $sample_size = $self->param('subsample_size');
  (my $subsample_name = $fastq_file_name) =~ s/fastq.gz/fq/;
  say "Subsample name is $subsample_name";
  $self->param('renamed_file', $subsample_name);
  #my $subsample = catfile('/nfs/production/flicek/ensembl/genebuild/do1/non_transcriptomic_registry/',$self->param('species'),'fastq',$subsample_name);
  #`mkdir $temp_path -p`;
  my $subsample = catfile($self->param('input_dir'),$subsample_name);
  my $fastq_file = catfile($self->param('input_dir'),$fastq_file_name);
  #unless(-e $fastq_file)  {
   #$self->complete_early("Transcriptomic data has been tested before so I will not download again");
   #$fastq_file = catfile($temp_path,$fastq_file_name);
   #$subsample = $subs;
   #$self->throw("Could not find fastq file on path. Path:\n".$fastq_file);
  #}
  say "Counting the reads in ".$fastq_file_name;
#  my $count_cmd = 'zcat '.$fastq_file.' | echo $((`wc -l`/4))';
 # my $count_reads = `$count_cmd`;
  #chomp($count_reads);

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
  # check if read is paired ended or not
  if ($self->param('is_paired')){
    #for paired ended reads, process both reads together when you find the first pair by substituting the suffix to 2
     #if ($subsample_name =~ m/_1/){
      # say "Subsample name is $subsample_name";
        # (my $mate = $subsample_name) =~ s/_1/_2/;
        # $subsample_name = $subsample_name ."\t$mate";
         #$self->param('paired',$subsample_name);
         #$self->param('flag',1); #flag to indicate the pair has now been processed, so best to skip it next time
       # 
      # }
      # my @csv = split(/_/,$subsample_name);
       #$subsample_name = $csv[0] . '_new.csv';
       #$subsample_name = catfile($self->param('input_dir'),$subsample_name);
       #$self->param('csv',$subsample_name);
  }
  else{
    #deal with single ended reads
    #$subsample_name = $subsample_name ."\t0";#set the pair to 0
    #$self->param('single',$subsample_name);
    #$self->param('flag',2);#this is to indicate single ended read processed
    #my @csv = split(/\./,$subsample_name);#extract read id from file name
   # $subsample_name = $csv[0] . '_new.csv';
   # $subsample_name = catfile($self->param('input_dir'),$subsample_name);
   # $self->param('csv',$subsample_name);
  }
  # say "Subsample name is $subsample_name";
  # my @csv = split(/_/,$subsample_name);
   #$subsample_name = $csv[0] . '_new.csv';
   #$subsample_name = catfile($self->param('input_dir'),$subsample_name);
   #remove downloaded fastq file
   `rm $fastq_file`;
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
  #my $output_file = $self->param('csv');
  #my $single_file = $self->param('single_csv');
  
  #if ($self->param('flag') == 1){
     #unless(open(PAIRED,">".$output_file)) {
    #   $self->throw("Failed to open the csv file for writing paired entries. Path used:\n".$output_file);
   #  }
   #  say PAIRED $self->param('paired');
   #  close PAIRED;
 #  }
 #  elsif ($self->param('flag') == 2){
   #  unless(open(SINGLE,">".$output_file)) {
    #   $self->throw("Failed to open the csv file for writing single entries. Path used:\n".$output_file);
     #}
     #say SINGLE $self->param('single');
    # close SINGLE;
   #}
  # else{}
}

1;
