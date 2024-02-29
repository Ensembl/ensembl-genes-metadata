=head1 LICENSE

 Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 Copyright [2016-2024] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::StarAlign

=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Pipeline::Runnable::StarAlign->new( );

$runnableDB->fetch_input();
$runnableDB->run();

=head1 DESCRIPTION

This module uses Star to align fastq to a genomic sequence

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Pipeline::Runnable::StarAlign;

use warnings;
use strict;
use feature 'say';
use File::Basename;
use File::Spec::Functions qw(catfile);
use Bio::EnsEMBL::Analysis::Runnable::Samtools;
use Bio::EnsEMBL::Analysis::Runnable::Star;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 param_defaults

 Arg [1]    : None
 Description: Default parameters
 Returntype : Hashref
                threads => 1,
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    threads => 1,
    samtools => 'samtools',
    samtools_use_threading => 1,
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Fetch parameters for STAR
 Returntype : None
 Exceptions : Throws if 'filename' does not exist
              Throws if 'fastqpair' does not exist
              Throws if 'wide_short_read_aligner' is not defined

=cut

sub fetch_input {
  my ($self) = @_;
  unless (-e $self->param('genome_dir') . 'Genome'){
    $self->throw("Genome file ".$self->param('genome_dir') . 'Genome'." not found\n");
  }
  unless (-d $self->param('output_dir')){
     my $dir = $self->param('output_dir');
     `mkdir $dir -p`;
  }
  my $input_ids = $self->param('SM');
  $self->say_with_header('Found '.scalar(@$input_ids).' input ids');
  foreach my $input_id (@$input_ids) {
    my $sample_id = $input_id->{'ID'};
    $sample_id  =~ s/^\s+//;
    $self->say_with_header("Processing sample: $sample_id");
    my $files = $input_id->{'files'};
    my $file1 = ${$files}[0];
    my $file2 = ${$files}[1];
    $file1 =~ s/^\s+//; $file2 =~ s/^\s+//;
    $self->say_with_header("Found file: $file1");
    my $filepath1 = $self->param('input_dir').'/'.$file1;
    $self->throw("Fastq file ".$filepath1." not found\n") unless ( -e $filepath1 );

    my $filepath2 = "";
    if($file2) {
      $self->say_with_header("Found paired file: $file2");
      $filepath2 = $self->param('input_dir').'/'.$file2;
      $self->throw("Fastq file ".$filepath2." not found\n") unless ( -e $filepath2 );
    }

    my $program = $self->param('short_read_aligner');
    $self->throw("Star program not defined in analysis\n") unless (defined $program);

    my $runnable = Bio::EnsEMBL::Analysis::Runnable::Star->new
    (
     -analysis       => $self->create_analysis,
     -program        => $program,
     -options        => $self->param('short_read_aligner_options'),
     -outdir         => $self->param('output_dir'),
     -genome_dir     => $self->param('genome_dir'),
     -genome     => $self->param('genome_dir'),
     -sample_name    => $sample_id,
     -fastq          => $filepath1,
     -fastqpair      => $filepath2,
     -threads        => $self->param('num_threads'),
    );
    if ($self->param_is_defined('rg_lines')) {
      $runnable->rg_lines($self->param('rg_lines'));
    }
    else {
      $runnable->rg_lines("ID:$sample_id");
    }
    $self->runnable($runnable);
  }
}


=head2 write_output

 Arg [1]    : None
 Description: Check that the output log file does not contain any error message.
              Check that the output BAM file is not truncated by running samtools flagstat.
              If no error nor truncation is found, dataflow the name of the resulting BAM file on branch 1 via 'filename'
 Returntype : None
 Exceptions : Throw if an error is found in the output log file or the BAM file is truncated.

=cut

sub write_output {
  my ($self) = @_;

  my $output_files = $self->output;
  foreach my $output_file (@$output_files) {
    my $output_file_basename = basename($output_file);
    my $output_file_dirname = dirname($output_file);
    my @output_file_basename_split = split('_',$output_file_basename,2);
    my $srr = shift(@output_file_basename_split);
    my $log_file = catfile($output_file_dirname,$srr.'_Log.out');
    my $log_file_ok = 1;
    my ($alignment_stats,$paired_stats);
    open(LOGFILE,$log_file) or die("Log file $log_file could not be opened.");
    while (my $string = <LOGFILE>) {
      if ($string =~ /Unexpected block structure/ or
          $string =~ /Possible output corruption/) {
        $log_file_ok = 0;
	last;
      }
    }
    close(LOGFILE) or die("Log file $log_file could not be closed.");

    if (!$log_file_ok) {
      $self->throw("'Unexpected block structure' or 'Possible output corruption' found in the log file $log_file");
    }

    if (-e $output_file) {
      my $samtools = Bio::EnsEMBL::Analysis::Runnable::Samtools->new(
                     -program => $self->param('samtools'),
                     -use_threading => $self->param('samtools_use_threading')
                     );
      $samtools->flagstat($output_file);# this will throw if the file is truncated
    } else {
      $self->throw("The output file $output_file does not exist. Cannot run samtools flagstat $output_file");
    }
    my $mapped = "samtools flagstat $output_file | awk -F \"[(|\%]\"  'NR == 5 {print \$2}'";
    $alignment_stats = `$mapped`;
    if ($alignment_stats =~ m/N\/A : N\/A\)/){
      $alignment_stats = 0;
    }
    
    $alignment_stats =~ s/\R\z//;
    my $paired = "samtools flagstat $output_file | awk -F \"[(|\%]\"  'NR == 9 {print \$2}'";
    $paired_stats = `$paired`;
    $paired_stats =~ s/\R\z//;
     if ($paired_stats =~ m /N\/A : N\/A\)/){
      $paired_stats = 0;
    }
    my $stats_file = catfile($self->param('output_dir'),'sr_flagstats.txt');
    open ID, (">>$stats_file");
    print ID $srr . "\t$alignment_stats\t$paired_stats\n";
    my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -user   => $ENV{GBUSER},
        -dbname => $self->param('pipe_db'),
        -host   => $ENV{GBS1},
        -port   => $ENV{GBP1},
        -pass   => $ENV{GBPASS},
        -driver => $ENV{GBDRIVER},
    );
     my $sth = $dba->dbc->prepare("insert into alignment_stats values (?,?,?,?)");
    $sth->bind_param(1,$self->param('accession'));
    $sth->bind_param(2,$srr);
    $sth->bind_param(3,$alignment_stats);
    $sth->bind_param(4,$paired_stats);
   
    if ($sth->execute){
            say "Mapping stats stored for alignment file $srr";
            `rm $output_file`;
          }
          else{
            $self->throw("Failed to store mapping stats for alignment $srr");
          }
    $self->say_with_header("Output file: ".$output_file);
    $self->dataflow_output_id([{'iid' => $output_file}], $self->param('_branch_to_flow_to'));
  }
}

1;
