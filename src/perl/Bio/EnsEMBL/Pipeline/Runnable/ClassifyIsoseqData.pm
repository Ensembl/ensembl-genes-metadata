#!/usr/bin/env perl
#
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License

#This script is designed to classify transcriptomic data following the traffic light system
#Once classification is complete, the results are stored in the database

#
package Bio::EnsEMBL::Pipeline::Runnable::ClassifyIsoseqData;


use strict;
use warnings;
use Getopt::Long;
use feature 'say';
use POSIX;
use TranscriptomicRegistryAdaptor;
use List::Util qw( min max );
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor;
use Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor;
use File::Basename;
use Data::Dumper;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub fetch_input{
	my ($self) = @_;
	
	#$self->param_required('flagstats');
#	$self->param_required('csv_file');
#	$self->param_required('read_length_file');
#	$self->param_required('fastq_report_file');
	$self->param_required('sp_id');
        $self->param_required('accession');
        $self->param_required('pipe_db');
  	$self->param_required('ass_id');
        say "Values are ",$self->param('ass_id');
        say "Values are ",$self->param('accession');
	$self->param('acc',$self->param('accession'));
}

sub run{
  my ($self) = @_;
  my $registry_dba = new Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor(
  -host    => 'mysql-ens-genebuild-prod-1',
  -port    => 4527,
  -user    => 'ensro',
  -dbname  => 'gb_assembly_registry',
  -pass    => '',
  -driver  => 'mysql',);
   my $ass = $self->param('ass_id');

  #my $sql = `gb2 gb_assembly_registry -e "SELECT CONCAT(chain,'.',version) FROM assembly WHERE assembly_id=$ass"`;
  #my @ass_id = split(/\n/,$sql);
  #my $accession = $self->param('ass_id');
  my $species = $self->param('sp_id');
  my $ass = $self->param('ass_id');
  my $db = $self->param('pipe_db'); 
  my $accession = $self->param('accession');
  #call function to retrieve alignment and stores details
  &store_alignment_details($ass,$species,$db,$accession);
 
  #begin classification of data
  
  &classify_rnaseq($species,$ass);
  sub classify_rnaseq{
    my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -user   => $ENV{GBUSER},
        -dbname => $ENV{REG_DB},
        -host   => $ENV{GBS1},
        -port   => $ENV{GBP1},
        -pass   => $ENV{GBPASS},
        -driver => $ENV{GBDRIVER},
    );
    my $status = ""; my %readcount;
    my $gcnt = 0; my $acnt = 0; my $ass_id = ""; my $sum_readcnt = 0;
    my $sth_summary = $dba->dbc->prepare("insert ignore into isoseq_data_summary (species_id, SM, total_read_count, avg_mapped_identity, status) values (?, ?, ? , ?, ?)");
    #my $sth_assembly = $dba->dbc->prepare("update assembly set rnaseq_data = ? where taxonomy = ? and assembly_id = ?");
    my $sth_csv = $dba->dbc->prepare("insert ignore into long_read_data (species_id, ID, SM, filename, url, md5sum, data_source, read_length) values (?,?,?,?,?,?,?,?)");
    #For long reads, we check for both read quality and mapping percentage to determine the sample classification
    my $sql = `mysql -h $ENV{GBS1} -P $ENV{GBP1} -u $ENV{GBUSER} -p$ENV{GBPASS} $ENV{REG_DB} -e \"select SM, round(sum(perc_mapped)/count(*),0) as avg_perc_mapped, isoseq_data_extra.data_quality, sum(read_count) as number_of_reads, isoseq_data_extra.assembly_id from isoseq_data_extra, isoseq_data_run_stats where project_id in (select project_id from rnaseq_info where species_id = $_[0]) and isoseq_data_extra.ID = isoseq_data_run_stats.ID and isoseq_data_extra.species_id = $_[0] and isoseq_data_extra.data_quality != 'fail' group by SM\"`;
    my @row = split (/\n/, $sql);
    for (my $y = 1; $y < @row; $y+=1){
      my @entry = split(/\t/,$row[$y]);
      $ass_id = $_[1];
      #check how many samples pass the good, weak and usable tests by looking at readcount, avg perc mapped, avg perc paired, average readlength and quality of read
        #code for testing
     if ($entry[1] > 70 && $entry[2] eq 'pass'){
        say "Sample ", $entry[0], " is good";
        $status = "good";
        $gcnt++;
        $readcount{$entry[0]} = $entry[3];#adding read count by sample name
        my $sql = `mysql -h $ENV{GBS1} -P $ENV{GBP1} -u $ENV{GBUSER} -p$ENV{GBPASS} $ENV{REG_DB} -e \"select ID,url,md5sum,data_source,calc_read_length from isoseq_data_extra where species_id = $_[0] and SM = '$entry[0]'\"`;
        my @record = split (/\n/, $sql);
        say "Total read retrieved is ",scalar(@record); 
        for (my $n = 1; $n < @record; $n+=1){
          my @line = split(/\t/,$record[$n]); 
            #here we deal with long reads
            $sth_csv->bind_param(1,$_[0]);
            $sth_csv->bind_param(2,$line[0]);
            $sth_csv->bind_param(3,$entry[0]);
            $sth_csv->bind_param(4,basename($line[1]));
            $sth_csv->bind_param(5,$line[1]);
            $sth_csv->bind_param(6,$line[2]);
            $sth_csv->bind_param(7,$line[3]);
            $sth_csv->bind_param(8,$line[4]);
            
            #Add eligible record to csv table
            if ($sth_csv->execute){
              say "Isoseq data successfully entered into csv table";
            }
            else{
              throw("Could not write to isoseq csv table");
            }
        }
      }
      #  #code for testing
      elsif (($entry[1] > 50 && $entry[2] eq 'pass') || ($entry[1] > 50 && $entry[2] eq 'warn')){
        say "Sample ", $entry[0], " is weak";
        $status = "weak";
        $acnt++;
        $readcount{$entry[0]} = $entry[3];#adding read count by sample name
        my $sql = `mysql -h $ENV{GBS1} -P $ENV{GBP1} -u $ENV{GBUSER} -p$ENV{GBPASS} $ENV{REG_DB} -e \"select ID,url,md5sum,data_source,calc_read_length from isoseq_data_extra where species_id = $_[0] and SM = '$entry[0]'\"`;
        my @record = split (/\n/, $sql);
        say "Total read retrieved is ",scalar(@record);
        for (my $n = 0; $n < @record; $n+=1){
          my @line = split(/\t/,$record[$n]);
            #here we deal with long ended reads
            $sth_csv->bind_param(1,$_[0]);
            $sth_csv->bind_param(2,$line[0]);
            $sth_csv->bind_param(3,$entry[0]);
            $sth_csv->bind_param(4,basename($line[1]));
            $sth_csv->bind_param(5,$line[1]);
            $sth_csv->bind_param(6,$line[2]);
            $sth_csv->bind_param(7,$line[3]);
            $sth_csv->bind_param(8,$line[4]);
            
            #Add eligible record to csv table
            if ($sth_csv->execute){
              say "Isoseq data successfully entered into csv table";
            }
            else{
              throw("Could not write to isoseq csv table");
            }
            
         }
      }
      else{
        
        say "Sample ", $entry[0], " is unusable $entry[1] $entry[2]";
        $status = "unusable";
      }
      $sth_summary->bind_param(1,$_[0]);
      $sth_summary->bind_param(2,$entry[0]);
      $sth_summary->bind_param(3,$entry[3]);
      $sth_summary->bind_param(4,$entry[1]);
      $sth_summary->bind_param(5,$status);
      if ($sth_summary->execute){
            say "Rnaseq summary successfully entered";
      }
      else{
        throw("Could not write to rnaseq summary table");
      }
         
    }
  }


  sub store_alignment_details{
  	my ($self) = @_;
    my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -user   => $ENV{GBUSER},
        -dbname => $ENV{REG_DB},
        -host   => $ENV{GBS1},
        -port   => $ENV{GBP1},
        -pass   => $ENV{GBPASS},
        -driver => $ENV{GBDRIVER},
   );
   
   my $pipe_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -user   => $ENV{GBUSER},
        -dbname => $_[2],
        -host   => $ENV{GBS1},
        -port   => $ENV{GBP1},
        -pass   => $ENV{GBPASS},
        -driver => $ENV{GBDRIVER},
   );

    my $taxonomy_adaptor = Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor->new(
              -user   => $ENV{GBUSER_R},
              -pass   => '',
              -dbname => 'ncbi_taxonomy_109',#$ENV{TAX_DB},
              -host   => 'mysql-ens-mirror-1',#$ENV{GBS1},
              -port   => 4240#$ENV{GBP1}
             );
    my $node_adaptor = $taxonomy_adaptor->get_TaxonomyNodeAdaptor();
    my $taxon_node = $node_adaptor->fetch_by_taxon_id($_[1]);
    my $rank = $taxon_node->rank();
    say "Taxonomy passed into alignment is ",$taxon_node->rank();
    #if ($rank eq 'subspecies'){
    #This is to allow for cases where rank could be anything other than a species
    #The query that normally fetches data from ENA would always fetch data from subordinate taxa where available 
    if ($rank ne 'species'){
      foreach my $ancestor ( @{ $node_adaptor->fetch_ancestors($taxon_node)}){
        if ($ancestor->rank eq 'species'){
          say "My rank is Species";    
          $rank = "species";
          last;
        }
      }
    }
    my %runs;
    my %mapped;
    my %calc;
    my %proj;
    my %readcnt;
    my %readlength;
    my %fqrep;
    my %rnaseq_info;
    my $species_id = $_[1];
   my $registry_adaptor = new TranscriptomicRegistryAdaptor(
        -user   => $ENV{GBUSER},
        -dbname => $db,
        -host   => $ENV{GBS1},#$self->param('pipe_host'),
        -port   => $ENV{GBP1},#$self->parma('pipe_port'),
        -pass   => $ENV{GBPASS},#$ENV{GBPASS},
        -driver => 'mysql',#$ENV{GBDRIVER},
    );
   say "Species id fed to long read is $species_id";
   #Retrieve metrics from pipe db and store in transcriptomic registry tables
   my $csv_entry = $pipe_db->dbc->prepare("select * from csv_long_read where taxon_id=?");
   my $fastqc_entry = $dba->dbc->prepare("select * from fastqc where species_id = ?");
   my $read_count_entry = $pipe_db->dbc->prepare("select * from read_count where species_id = ?");
   my $alignment_stats = $pipe_db->dbc->prepare("select * from alignment_stats where accession = ?");
 
       my $sth_info = $dba->dbc->prepare("select project_id,species_id from rnaseq_info");
       my $sth_rnaseq = $dba->dbc->prepare("insert into rnaseq_info (species_id, project_id, data_source) values (?, ?, ?)");
       my $sth_extra = $dba->dbc->prepare("insert ignore into isoseq_data_extra values (?,?,?,?,?,?,?,?,?,?,?,?)");
       my $sth_stats = $dba->dbc->prepare("insert ignore into isoseq_data_run_stats values (?, ?, ?, ?)");
       if ($sth_info->execute){
          while (my @entry = $sth_info->fetchrow_array()){
            $rnaseq_info{$entry[0] . '_' . $entry[1]} = $entry[1];
          }
        }
    my $run_id;
    $read_count_entry->bind_param(1,$species_id);
    if ($read_count_entry->execute){
          while (my @entry = $read_count_entry->fetchrow_array()){
            if ($entry[0] =~ m/_/){
               my @val = split(/_/,$entry[0]);
               $run_id = $val[0];
            }
            else{
               my @val = split(/\./,$entry[0]);
               $run_id = $val[0];
            }
            $readcnt{$run_id} = $entry[1];
            $readlength{$run_id} = $entry[3];
         }
    }
    $fastqc_entry->bind_param(1,$species_id);
    if ($fastqc_entry->execute){
          while (my @entry = $fastqc_entry->fetchrow_array()){
            if ($entry[1] =~ m/_/){
               my @val = split(/_/,$entry[1]);
               $run_id = $val[0];
            }
            $fqrep{$run_id} = lc($entry[2]);
         }
    }
    say "Accession is ",$_[3];
    my $accession = $_[3];
    $alignment_stats->bind_param(1,$accession);
    if ($alignment_stats->execute){
          while (my @entry = $alignment_stats->fetchrow_array()){
            say "Alignment here is $entry[1]";
            $mapped{$entry[1]} = $entry[2];
         }
    }
    #get read entry from csv table
    my %csv;
    my @results;
    $csv_entry->bind_param(1,$species_id);
    if ($csv_entry->execute){
      say "CSV record retrieval successful";
    }
    #Loop through each record and register accordingly
    while (my @entry = $csv_entry->fetchrow_array()){
      push @results, join("\t ", @entry);
    }
   foreach my $result (@results) {
      $result =~ s/[^[:ascii:]|[:alnum:]]//g;
      my @arr = split(/\t/,$result);
      $arr[1] =~ s/fastq.gz/fq/;
      $arr[1] =~ s/^\s+|\s+$//g;
      $csv{$arr[1]} = $result;
      say "read is $arr[1]";
      my @tem = split(/,/,$arr[3]);
      my $run_id = "";
      say "species id is ", length($arr[2]) . " " . $arr[1];
      if ($arr[1] =~ m/_/){
        my @val = split(/_/,$arr[1]);
        $run_id = $val[0];
      }
      else{
        my @val = split(/\./,$arr[1]);
        $run_id = $val[0];
      }
      #If read was not tested in alignment, skip
      if (!exists $fqrep{$run_id}){
        next;
      }
      if (!exists $runs{$run_id}){
        my $read_source = "";
        $runs{$run_id} = 1;
        say "Readcount is ",$readcnt{$run_id};
        if ($readcnt{$run_id} !~ /\d+/){
          say "run does not exist";
          next;
        }
        my $proj_spec_id = $tem[0] . '_' . $_[1];
        say "Combined primary id is $proj_spec_id";
        if ((exists($rnaseq_info{$proj_spec_id}) && ($rnaseq_info{$proj_spec_id} !~ m/$_[1]/)) or (!exists($rnaseq_info{$proj_spec_id}))){
          $sth_rnaseq->bind_param(1,$_[1]);
          say "Species tested is $_[1]";
          $sth_rnaseq->bind_param(2,$tem[0]);
        }
        if ($species_id != $arr[7]){#where species id differs with entry in csv file, check that species is subspecies. If so, data level remains species
          if ($rank eq 'subspecies'){
            $sth_rnaseq->bind_param(3,'species');
            $read_source = 'species';
          }
          else{
            $sth_rnaseq->bind_param(3,'genus');
            $read_source = 'genus';
          }
        }
        else{
          $sth_rnaseq->bind_param(3,'species');
          $read_source = 'species';
        }
        $sth_rnaseq->bind_param(3,$read_source);
        say "Project tested is $tem[0]";
        if (!exists $rnaseq_info{$proj_spec_id}){
           say "Species again tested is $species_id";
           if ($sth_rnaseq->execute){
            $rnaseq_info{$proj_spec_id} = $species_id;
          }
          else{
            throw("Could not write mapped data to rnaseq info table");
          }
        }
        $arr[3] = substr($arr[3],0,250);
        $arr[3] =~ tr/:\t/ /;
        $arr[0] = substr($arr[0],0,50);
        my $is_mate_1 = 0;
        say "Run id being used is $run_id";
        say "length of read is ",$readlength{$run_id};
        chomp($fqrep{$run_id});
        $sth_extra->bind_param(1,$tem[0]);
        $sth_extra->bind_param(2,$run_id);
        $sth_extra->bind_param(3,$_[0]);
        $sth_extra->bind_param(4,$arr[0]);
        $sth_extra->bind_param(5,$readcnt{$run_id});
        $sth_extra->bind_param(6,$arr[3]);
        $sth_extra->bind_param(7,$species_id);
        $sth_extra->bind_param(8,$fqrep{$run_id});
        $sth_extra->bind_param(9,$arr[4]);
        $sth_extra->bind_param(10,$arr[5]);
        $sth_extra->bind_param(11,$read_source);
        $sth_extra->bind_param(12,$readlength{$run_id});
        if ($sth_extra->execute){
          
        }
        else{
          throw("Could not write mapped to rnaseq extra table");
        }
        $sth_stats->bind_param(1,$_[0]);
        $sth_stats->bind_param(2,'test');
        say "Run ID is $run_id";
        say "Value from mapped hash is ",$mapped{$run_id};
        #checks when percentage identity is not set, assign 0
        if ($mapped{$run_id} =~ m/N/ || $mapped{$run_id} eq ''){
          $sth_stats->bind_param(3,0);
        }else{
          $sth_stats->bind_param(3,$mapped{$run_id});
        }
        $sth_stats->bind_param(4,$run_id);
        
        if ($sth_stats->execute){
          
        }
        else{
          throw("Could not write mapped to rnaseq data run table");
        }
        
      }
      
    }
    
  }

}

sub write_output{
  my ($self) = @_;
}



1;
