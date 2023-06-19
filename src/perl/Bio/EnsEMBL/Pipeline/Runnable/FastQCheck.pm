=head1 LICENSE

Copyright [2018-2019] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
:128

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

package Bio::EnsEMBL::Pipeline::Runnable::FastQCheck;

use strict;
use warnings;
use feature 'say';
use File::Spec::Functions;
use File::Basename;
use POSIX;
use TranscriptomicRegistryAdaptor;

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

  $self->param_required('quality_checker');
  $self->param_required('input_dir');
  $self->param_required('iid');
  $self->param_required('is_paired');
  $self->param_required('species');
  $self->param('source_id',7215);
    
  
}


sub run {
  my ($self) = @_;
  my $fastqc = $self->param('quality_checker');
  my $file = $self->param('input_dir') . $self->param('iid');
  my $output_path = $self->param('input_dir');
  $file =~ s/fastq.gz/fq/;
  #command to run fastqc
  `$fastqc $file -q --outdir $output_path --extract`;
  my @run = split(/\./,$self->param('iid'));
  my $summary = $self->param('input_dir') . $run[0] .'_fastqc/summary.txt';
  
  #get run id without suffix for naming file to contain paired read ids
  my @csv_file = split(/_/,$run[0]);
  
  $self->param('csvfile', $csv_file[0]);
  my $cmd = "";
  my %report; 
  #read the summary report file and get the required status value
  open T, ("<$summary");
  while (<T>){
    if ($_ !~ m/Basic Statistics/){
      #store all report lines in hash
      my @status = split(/\t/,$_);
      $report{$status[1]} = $self->param('iid') . "\t" . $status[0];
    }
  }
  close (T);
  if ($self->param('is_paired') == 1 || $self->param('is_paired') == 0){
    $cmd = $self->param('input_dir') . $run[0] . '_fastqc ' . $self->param('input_dir') . $run[0] . '_fastqc.html ' . $self->param('input_dir') . $run[0] . '_fastqc.zip ' . $self->param('csvfile') .  '_sr.csv';;
  }
  else{
    $cmd = $self->param('input_dir') . $run[0] . '_fastqc ' . $self->param('input_dir') . $run[0] . '_fastqc.html ' . $self->param('input_dir') . $run[0] . '_fastqc.zip ' . $self->param('input_dir') . $self->param('csvfile') .  '_lr.csv';
  } 
  $self->param('db_report',\%report);
}


sub write_output {
  my ($self) = @_;
  my $file = $self->param('input_dir') . $self->param('iid');
  my $output_file = "";
  my $fqcrep = $self->param('input_dir') .'fastqc_report.txt';
  if ($self->param('is_paired') == 1 || $self->param('is_paired') == 0){
  	  $output_file = $self->param('input_dir') . $self->param('csvfile') . "_sr.csv";
  }
  else{
  	  $output_file = $self->param('input_dir') . $self->param('csvfile') . "_lr.csv";
  }
  open QC, (">>$fqcrep");
  my @report;
  my $registry_adaptor = new TranscriptomicRegistryAdaptor(
        -user   => $ENV{GBUSER},
        -dbname => $ENV{REG_DB},
        -host   => $ENV{GBS1},
        -port   => $ENV{GBP1},
        -pass   => $ENV{GBPASS},
        -driver => 'mysql',
    );
   my ($qc_rep,$rep_line);
   my $sth = $registry_adaptor->dbc->prepare("insert ignore into fastqc values (?,?,?,?,?,?,?,?,?,?,?,?)");
   $sth->bind_param(1,$self->param('taxon_id'));
   while (($rep_line, $qc_rep) = each ($self->param('db_report'))){
    @report = split(/\t/,$qc_rep);
    if ($rep_line =~ /Per base sequence quality/){
      $sth->bind_param(2,$report[0]);
      $sth->bind_param(3,$report[1]);
      #For now we still want to write to report file
      say QC $report[0] . "\t" . $report[1];
      if ($report[1] !~ m/pass|warn/i) {
        $self->param('flow', 1);
      }
    }
    elsif ($rep_line =~ /Per sequence quality scores/){
      $sth->bind_param(4,$report[1]);
    }
    elsif ($rep_line =~ /Per base sequence content/){
      $sth->bind_param(5,$report[1]);
    }
    elsif ($rep_line =~ /Per sequence GC content/){
      $sth->bind_param(6,$report[1]);
    }
    elsif ($rep_line =~ /Per base N content/){
      $sth->bind_param(7,$report[1]);
    }
    elsif ($rep_line =~ /Sequence Length Distribution/){
      $sth->bind_param(8,$report[1]);
    }
    elsif ($rep_line =~ /Sequence Duplication Levels/){
      $sth->bind_param(9,$report[1]);
    }
    elsif ($rep_line =~ /Overrepresented sequences/){
      $sth->bind_param(10,$report[1]);
    }
    elsif ($rep_line =~ /Adapter Content/){
      $sth->bind_param(11,$report[1]);
    }
    else{}
   }
   $sth->bind_param(12,$self->param('source_id'));
   if ($sth->execute){
     say "Data quality information stored for fastq file $report[0]";
   }
   else{
     $self->throw("Failed to store data quality information for file $report[0]");
   }
  if ($self->param('is_paired') < 2){
  	  
     my $read = $self->param('iid');
     $read =~ s/fastq.gz/fq/;
     if ($read =~ m/_1/){
       unless(open(PAIRED,">".$output_file)) {
         $self->throw("Failed to open the csv file for writing paired entries. Path used:\n".$output_file);
       }
      #setting paired information;
        (my $mate = $read) =~ s/_1/_2/;
        $read = $read ."\t$mate";
        unless($read =~ m/\.fq$/) {
          $self->throw("Can't write to csv file. Check value:\n".$read);
        }
        say PAIRED $read;
        close PAIRED;
     }
   }
   else{
   	   say "file is $file";
   	 if(-e $file.'.fai') {
   	 }
     unless(open(SINGLE,">".$output_file)) {
       $self->throw("Failed to open the csv file for writing single entries. Path used:\n".$output_file);
     }
     say SINGLE $self->param('iid');
     close SINGLE;
   }
   
  
}

1;
