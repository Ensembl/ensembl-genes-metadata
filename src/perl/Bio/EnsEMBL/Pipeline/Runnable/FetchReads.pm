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

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::FetchReads

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Pipeline::Runnable::FetchReads;

use warnings;
use strict;
use feature 'say';
use POSIX;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub param_defaults {
  my ($self) = @_;

  return {
          %{$self->SUPER::param_defaults},
    decompress => 0,
    create_faidx => 0,
  }
} 

sub exit_code_test {
  my ($self,$wget_cmd) = @_;

  my $res = $self->run_system_command($wget_cmd);
  if ($res) {
    $res >>= 8;
    $self->warning("wget died with error code $res");
    return 0;
  } else {
    return 1;
  }

}

#sub to check if transcriptomic data has been tested before
sub is_rnaseq_classified{
  my ($self) = @_;
    my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -user   => $ENV{GBUSER},
        -dbname => $ENV{REG_DB},
        -host   => $ENV{GBS1},
        -port   => $ENV{GBP1},
        -pass   => $ENV{GBPASS},
        -driver => $ENV{GBDRIVER},
    );
    my $status = 0;
    my $species_id = $_[1];
    my @ID = split (/_/, $_[0]);
    my $sth;
    if ($_[2] == 1){
      $sth = $dba->dbc->prepare("select * from tested_read_data where ID = ? and taxon_id = ?");
    }
    else{
      $sth = $dba->dbc->prepare("select * from tested_read_data where ID = ? and taxon_id = ?");
    }
    $sth->bind_param(1,$ID[0]);
    $sth->bind_param(2,$species_id);
    
    #if short read data has been tested before for this species, set value = 1 and we don't re-download it again
    if ($sth->execute){
      while (my @result = $sth->fetchrow_array()){
	say "File returned is ",@result;
        $status = 1;
      }
    }
    
    return $status;
}

sub decompress {
  my ($self,$path,$fastq) = @_;
  my $cmd = 'gunzip '.$path.'/'.$fastq;

  # Remove this in case indexing in the code block after this one
  if($fastq =~ s/\.gz$//) {
    my $gunzip_res = system($cmd);
    if($gunzip_res) {
      $self->throw("Failed to decompress file. Command:\n".$cmd);
    }
  } else {
    $self->warning("You selected decompress, but the file did not have a .gz extension, so will not try and decompress");
  }

  # Update these in case the extension was removed
  $self->param('iid',$fastq);
  $self->param('fastq_file',$fastq);

  return($fastq);
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

sub fetch_input {
  my ($self) = @_;
  $self->param_required('filename'); 
  $self->param_required('taxon_id');
  my $rnaseq_data_status = is_rnaseq_classified($self->param('filename'),$self->param('taxon_id'),$self->param('is_paired'));
  if ($rnaseq_data_status == 1) {
    $self->complete_early("Transcriptomic data has been tested before so I will not download again");
  } 
  
}

sub write_output {
  my ($self) = @_;
  my $ftp_base_url = $self->param('ftp_base_url');
  my $fastq = $self->param('filename');
  my $path = catdir($self->param('input_dir'),$self->param('species'),"/fastq/");
  my @file;
  unless (-d $path){
  	  `mkdir -p $path`;
  }
  my $srr;
   
  if ($fastq =~ m/_/){
    $srr = (split /_/, $fastq)[0];
  }
  else{
    $srr = (split /\./, $fastq)[0];
  }
  my $first = substr $srr, 0, 6;
  my $second = "";
  if (length($srr) > 10){ 
   $second = '0'.(substr $srr, -2, 2);
  }else{
   $second = '00'.(substr $srr, -1, 1);
  }
  say "About to download file ", "$ftp_base_url/$first/$srr/$fastq",  '-P', $path;
  my $second_a = '00'.(substr $srr, -1, 1);
  my $second_b ='0'.(substr $srr, -2, 2);

  my $exit_code = 0;
  my $fastq_ftp = $self->param('url');
  my $wget_cmd_list = [['wget', '-c', '-qq', "$fastq_ftp", '-P', $path],['wget', '-c', '-qq', "$fastq_ftp", '-P', $path],['wget', '-c', '-qq', "$fastq_ftp", '-P', $path]];
  foreach my $wget_cmd (@$wget_cmd_list){
    $exit_code = $self->exit_code_test($wget_cmd);
    if ($exit_code){
      last;
    }
  }
  if (!$exit_code){
    if (-e $path.'/'.$fastq) {
      $self->run_system_command(['rm',"$path/$fastq"]);
    }
    $self->throw("Failed to download $fastq");
  }
  unless(-e $path.'/'.$fastq) {
     #$self->complete_early("Transcriptomic data has been tested before so I will not download again");
    $self->throw("Did not find the fastq file on the expected path. Path:\n".$path."/".$fastq);
  }
 #Add downloaded fastq files to csv table in pipeline db
  my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -user   => $ENV{GBUSER},
        -dbname => $self->param('pipe_db'),
        -host   => $ENV{GBS1},
        -port   => $ENV{GBP1},
        -pass   => $ENV{GBPASS},,
        -driver => 'mysql',
    );
    my $sth;
    my ($pair) = $self->param('filename') =~ /\S+_(\d)\.\S+/;
    my $source = '';
    if (!defined($self->param('source_id'))){
          $source = $self->param('taxon_id');
        }
        else{
          $source = $self->param('source_id');
        }
    if (defined($pair) and ($self->param('is_paired') ==  1)){
            if (($pair == 1) or ($self->param('is_paired') == 0)) {
                $self->param('is_mate_1',1);
            }
            else {
                $self->param('is_mate_1',0);
            }
    	$sth = $dba->dbc->prepare("insert  into csv_short_read values (?,?,?,?,?,?,?,?,?,?)");
    	$sth->bind_param(1,$self->param('SM'));
    	$sth->bind_param(2,$self->param('ID'));
    	$sth->bind_param(3,$self->param('is_paired'));
    	$sth->bind_param(4,$self->param('taxon_id'));
    	$sth->bind_param(5,$self->param('filename'));
    	$sth->bind_param(6,$self->param('DS'));
    	$sth->bind_param(7,$self->param('url'));
    	$sth->bind_param(8,$self->param('md5'));
    	$sth->bind_param(9,$self->param('species'));
        $sth->bind_param(10,$source);
        
     }
     else{
       #Long read data being processed
        my $subsample = $self->param('filename');
        $subsample =~ s/fastq.gz/fq/;
        $sth = $dba->dbc->prepare("insert into csv_long_read values (?,?,?,?,?,?,?,?)");
        $sth->bind_param(1,$self->param('SM'));
        $sth->bind_param(2,$self->param('filename'));
        $sth->bind_param(3,$self->param('taxon_id'));
        $sth->bind_param(4,$self->param('DS'));
        $sth->bind_param(5,$self->param('url'));
        $sth->bind_param(6,$self->param('md5'));
        $sth->bind_param(7,$self->param('species'));
        $sth->bind_param(8,$self->param('source_id'));
     }
    if ($sth->execute){
        print("Downloaded fastq files registered in csv table");
    }
    else{
     	print("Failed to register downloaded fastq file");
    }
  if ($self->param('is_paired')){
    $self->dataflow_output_id({iid=>$self->param('filename'),species=>$self->param('species'),is_paired=>$self->param('is_paired')}, "2");
  }
  else{
    #To handle long reads
    $self->dataflow_output_id({iid=>$self->param('filename'),species=>$self->param('species'),is_paired=>2}, "2");
  }
}

1;
