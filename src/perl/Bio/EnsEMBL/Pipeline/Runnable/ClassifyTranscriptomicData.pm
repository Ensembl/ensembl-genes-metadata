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
package Bio::EnsEMBL::Pipeline::Runnable::ClassifyTranscriptomicData;

use strict;
use warnings;
use Getopt::Long;
use feature 'say';
use POSIX;
use List::Util qw( min max );
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor;
use File::Basename;
use Data::Dumper;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub fetch_input{
	my ($self) = @_;
	
	
}

sub run{
  my ($self) = @_;
  my $accession = $self->param('accession'); 
  my $species = $self->param('sp_id');
  my $ass = $self->param('ass_id');
  my $db = $self->param('pipe_db');
  
  #call function to retrieve alignment and stores details
  &store_alignment_details($species,$ass,$db,$accession);
 
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
    my $sth_rnaseq = $dba->dbc->prepare("insert ignore into rnaseq_data_summary (species_id, SM, total_read_count, avg_mapped_identity, status, avg_paired_identity) values (?, ?, ?, ? , ?, ?)");
    my $sth_assembly = $dba->dbc->prepare("update assembly set rnaseq_data = ? where taxonomy = ? and assembly_id = ?");
    my $sth_csv = $dba->dbc->prepare("insert ignore into short_read_data (species_id, ID, SM, filename, is_paired, read_length, url, md5sum, data_source) values (?,?,?, ?, ?, ?, ?, ?, ?)");
    my $sql = `mysql -h $ENV{GBS1} -P $ENV{GBP1} -u $ENV{GBUSER} -p$ENV{GBPASS} $ENV{REG_DB} -e \"select sum(read_count) as number_of_reads, SM, round(sum(perc_mapped)/count(*),0) as avg_perc_mapped, rnaseq_data_extra.data_quality, rnaseq_data_extra.assembly_id, round(sum(calc_read_length)/count(*),0) as avg_read_length, round(sum(perc_paired)/count(*),0) as avg_perc_paired from rnaseq_data_extra, rnaseq_data_run_stats where project_id in (select project_id from rnaseq_info where species_id = $_[0]) and rnaseq_data_extra.ID = rnaseq_data_run_stats.ID and rnaseq_data_extra.species_id = $_[0] and rnaseq_data_extra.data_quality != 'fail' and rnaseq_data_extra.assembly_id = $_[1] and rnaseq_data_run_stats.assembly_id = $_[1]  group by SM\"`;
    my @row = split (/\n/, $sql);
    say "Record count = ",scalar(@row);
    for (my $y = 1; $y < @row; $y+=1){
      my @entry = split(/\t/,$row[$y]);
      $ass_id = $_[1];
      #check how many samples pass the good, weak and usable tests by looking at readcount, avg perc mapped, average readlength and quality of read
        #code for testing
     if ($entry[0] > 50000000 && $entry[2] > 70 && $entry[3] eq 'pass' && $entry[5] >= 90){
        say "Sample ", $entry[1], " is good";
        say "Sample ", $entry[1], " is good $entry[0] $entry[2] $entry[3] $entry[5]";
        $status = "good";
        $readcount{$entry[1]} = $entry[0];#adding read count by sample name
        $gcnt++;
        my $sql = `mysql -h $ENV{GBS1} -P $ENV{GBP1} -u $ENV{GBUSER} -p$ENV{GBPASS} $ENV{REG_DB} -e \"select ID,is_paired,calc_read_length,url,md5sum,data_source from rnaseq_data_extra where species_id = $_[0] and SM = '$entry[1]'\"`;
        my @record = split (/\n/, $sql);
        
        for (my $n = 1; $n < @record; $n+=1){
        
          my @line = split(/\t/,$record[$n]); 
          my $rcnt = 0;
          if ($line[1] == 1){
            #deal with is_paired eneded reads
            my $read1 = $line[0] . '_1.fastq.gz' . ';' . $line[0] . '_2.fastq.gz';
            
            my @query = ([$_[0],$line[0],$entry[1],$read1,$line[1],$line[2],$line[3],$line[4],$line[5]]);
           foreach my $rec ( @query ) { 
             if ($sth_csv->execute( @$rec )){
             }
           }
            
          }
            
        }
      }
      #  #code for testing
      elsif (($entry[0] > 30000000 && $entry[2] > 50 && $entry[3] eq 'pass' && $entry[5] >= 70) || ($entry[0] > 30000000 && $entry[2] > 50 && $entry[3] eq 'warn' && $entry[5] >= 70)){
        say "Sample ", $entry[1], " is weak";
        say "Sample ", $entry[1], " is weak $entry[0] $entry[2] $entry[3] $entry[5]";
        $status = "weak";
        $acnt++;
        $readcount{$entry[1]} = $entry[0];#adding read count by sample name
        my $sql = `mysql -h $ENV{GBS1} -P $ENV{GBP1} -u $ENV{GBUSER} -p$ENV{GBPASS} $ENV{REG_DB} -e \"select ID,is_paired,calc_read_length,url,md5sum,data_source from rnaseq_data_extra where species_id = $_[0] and SM = '$entry[1]'\"`;
        my @record = split (/\n/, $sql);
        say "Record is ",@record;
        for (my $n = 1; $n < @record; $n+=1){
          my @line = split(/\t/,$record[$n]); 
          my $rcnt = 0;
          if ($line[1] == 1){
            #deal with is_paired eneded reads
            my $read1 = $line[0] . '_1.fastq.gz' . ';' . $line[0] . '_2.fastq.gz';
            $sth_csv->bind_param(1,$_[0]);
            $sth_csv->bind_param(2,$line[0]);
            $sth_csv->bind_param(3,$entry[1]);
            $sth_csv->bind_param(4,$read1);
            $sth_csv->bind_param(5,$line[1]);
            $sth_csv->bind_param(6,$line[2]);
            $sth_csv->bind_param(7,$line[3]);
            $sth_csv->bind_param(8,$line[4]);
            $sth_csv->bind_param(9,$line[5]);
             if ($sth_csv->execute){
               say "We entered successfully";
             }
            
          }
            
        }
      }
      else{
        
        say "Sample ", $entry[1], " is unusable $entry[0] $entry[2] $entry[3] $entry[5]";
        $status = "unusable";
      }
      #Store data summary
      $sth_rnaseq->bind_param(1,$_[0]);
      $sth_rnaseq->bind_param(2,$entry[1]);
      $sth_rnaseq->bind_param(3,$entry[0]);
      $sth_rnaseq->bind_param(4,$entry[2]);
      $sth_rnaseq->bind_param(5,$status);
      $sth_rnaseq->bind_param(6,$entry[6]);
      if ($sth_rnaseq->execute){
            say "Rnaseq summary successfully entered";
      }
      else{
        throw("Could not write to rnaseq summary table");
      }
         
    }
   
    if ($gcnt >= 5){#if  5 or more samples meet the good criteria, then set species green - automate annotation
      $sth_assembly->bind_param(1,'green');
      $sth_assembly->bind_param(2,$_[0]);
      $sth_assembly->bind_param(3,$ass_id);
      
    }
    elsif ($gcnt < 5){#if not enough good samples, get read counts for all good and weak samples
      while( my( $key, $value ) = each %readcount){
        $sum_readcnt += $value;#summing up read counts for all species specific good and weak samples
        $acnt++;
      }
      #checking to see that all amber reads sum up to at least 100,000,000 reads;
      if (($acnt > 1) && ($sum_readcnt >= 100000000)){#if sum of all good and weeak samples > 100,000,000 and there are > 1 samples, then assign sample group amber
        $sth_assembly->bind_param(1,'amber');
        $sth_assembly->bind_param(2,$_[0]);
        $sth_assembly->bind_param(3,$ass_id);
        
      }
      else{#otherwise assign red
        $sth_assembly->bind_param(1,'red');
        $sth_assembly->bind_param(2,$_[0]);
        $sth_assembly->bind_param(3,$ass_id);
        
      }
    }
    else{}
    if ($sth_assembly->execute){
      say "Rnaseq classification successfully entered";	  
    }
    else{
      throw("Could not write to assembly table");
    }
   my $rnaseq_data_source;
   #get summary source of transcriptomic data
   my $data_source = $dba->dbc->prepare("select distinct(data_source) from rnaseq_info where species_id = ?");
   $data_source->bind_param(1,$_[0]);
   if ($data_source->execute){
     if ($data_source->rows > 1){
       $rnaseq_data_source = 'hybrid';
     }
     else{
       $rnaseq_data_source = $data_source->fetchrow();
     }
   }
   #Update assembly table with source of transcriptomic data
   my $sth_data_source = $dba->dbc->prepare("update assembly set rnaseq_data_source = ? where taxonomy = ?");
   $sth_data_source->bind_param(1,$rnaseq_data_source);
   $sth_data_source->bind_param(2,$_[0]);
   if ($sth_data_source->execute){
     say "data source updated in assembly table";
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
              -dbname => 'ncbi_taxonomy_109',
              -host   => $ENV{MIRROR},
              -port   => $ENV{MPORT}
             );
    my $node_adaptor = $taxonomy_adaptor->get_TaxonomyNodeAdaptor();
    my $taxon_node = $node_adaptor->fetch_by_taxon_id($_[0]);
    my $species_rank = $taxon_node->rank();
    my %runs;
    my %mapped;
    my %is_paired;
    my %calc;
    my %proj;
    my %readcnt;
    my %fqrep;
    my %readlength;
    my %rnaseq_info1;
    my $species_id = $_[0];
   #Query to retrieve temp data stored in pipeline db
   my $csv_entry = $pipe_db->dbc->prepare("select * from csv_short_read where taxon_id=?");
   my $fastqc_entry = $dba->dbc->prepare("select * from fastqc where species_id = ?");
   my $read_count_entry = $pipe_db->dbc->prepare("select * from read_count where species_id = ?");
   my $alignment_stats = $pipe_db->dbc->prepare("select * from alignment_stats where accession = ?");
   
   #Query to write to meta database
   my $sth_info = $dba->dbc->prepare("select project_id,species_id from rnaseq_info");
   my $sth_rnaseq = $dba->dbc->prepare("insert into rnaseq_info (species_id, project_id, data_source) values (?, ?, ?)");
   my $sth_extra = $dba->dbc->prepare("insert ignore into rnaseq_data_extra values (?,?,?,?,?,?,?,?,?,?,?,?)");
   my $sth_stats = $dba->dbc->prepare("insert ignore into rnaseq_data_run_stats values (?, ?, ?, ?, ?)");
   my $sth_tested_reads = $dba->dbc->prepare("insert ignore into tested_read_data (ID, taxon_id) values (?, ?)"); 
   #Execute query to fetch all existing records per species by project id
   if ($sth_info->execute){
      while (my @entry = $sth_info->fetchrow_array()){
         $rnaseq_info1{$entry[0] . '_' . $entry[1]} = $entry[1];
      }
   }

   #Get readcount per fastq for species
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
   
   #Get fastqc report per run id
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
   
   #Get alignment report per run
   $alignment_stats->bind_param(1,$_[3]);
   if ($alignment_stats->execute){
       while (my @entry = $alignment_stats->fetchrow_array()){
            $mapped{$entry[1]} = $entry[2];# . "\t" . $entry[3];
            $is_paired{$entry[1]} = $entry[3];
       }
   }

   #Get list of fastq files used for alignment from csv table
   my %csv;
   my @pair;
   my @results;
   $csv_entry->bind_param(1,$species_id);
   if ($csv_entry->execute){
      say "CSV record retrieval successful";
   }
   while (my @entry = $csv_entry->fetchrow_array()){
      push @results, join("\t ", @entry);
      $csv{$entry[4]} = join("\t ", @entry);
   }
    foreach my $result (@results) {
      $result =~ s/[^[:ascii:]|[:alnum:]]//g;
      my @arr = split(/\t/,$result);
      $arr[1]=~ s/^\s+|\s+$//g;$arr[9]=~ s/^\s+|\s+$//g;$arr[5]=~ s/^\s+|\s+$//g;
      if ($arr[4] =~ /_1.fastq.gz/){
        (my $mate = $arr[4]) =~ s/_1/_2/;
        chomp($mate);
        $mate =~ s/^\s+|\s+$//g;
        @pair = split(/\t/,$csv{$mate});
        
      }
      else{next;}
      $arr[4] =~ s/fastq.gz/fq/;
      $arr[4] =~ s/^\s+|\s+$//g;
      $csv{$arr[4]} = $result;
      my @tem = split(/,/,$arr[5]);
      if ($arr[4] =~ m/_/){
        my @val = split(/_/,$arr[4]);
        $run_id = $val[0];
      }
      else{
        my @val = split(/\./,$arr[4]);
        $run_id = $val[0];
      }
      #If read was not tested in alignment, skip
      if (!exists $fqrep{$run_id}){
        next;
      }
      if (!exists $runs{$run_id}){
        my $read_source = "";
        $runs{$run_id} = 1;
        if ($readcnt{$run_id} !~ /\d+/){
          say "run does not exist";
          next;
        }
        $tem[0] =~ s/^\s+|\s+$//g;
        my $proj_spec_id = $tem[0] . '_' . $_[0];
        if ((exists($rnaseq_info1{$proj_spec_id}) && ($rnaseq_info1{$proj_spec_id} !~ m/$_[0]/)) or (!exists($rnaseq_info1{$proj_spec_id}))){
          $sth_rnaseq->bind_param(1,$_[0]);
          $sth_rnaseq->bind_param(2,$tem[0]);
        }
        
        #Check for the source of data
        $arr[9] =~ s/^\s+|\s+$//g;
        $node_adaptor = $taxonomy_adaptor->get_TaxonomyNodeAdaptor();
        $taxon_node = $node_adaptor->fetch_by_taxon_id($arr[9]);
	#Because this particular species has some taxonomy related issue, we have to manually set its value
        if (!($taxon_node) and $arr[9] == 3016450){
          $taxon_node = $node_adaptor->fetch_by_taxon_id(9539);
        }
        else{
          $taxon_node = $node_adaptor->fetch_by_taxon_id($taxon_node->parent_id);
        }
        my $parent_rank = $taxon_node->rank();
        if ($species_id != $arr[9]){#where species id differs with entry in csv file, check that species is subspecies. If so, data level remains species
          if ($species_rank eq 'subspecies'){
            $sth_rnaseq->bind_param(3,'species');
            $read_source = 'species';
          }
          elsif ($parent_rank eq 'species'){
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
        if (!exists $rnaseq_info1{$proj_spec_id}){
           if ($sth_rnaseq->execute){
            $rnaseq_info1{$proj_spec_id} = $species_id;
          }
          else{
            throw("Could not write mapped data to rnaseq info table");
          }
        }
        $arr[5] = substr($arr[5],0,250);
        $arr[5] =~ tr/:\t/ /;
        $arr[0] = substr($arr[0],0,50);
        my $is_mate_1 = 0;
        chomp($fqrep{$arr[1]});
        $pair[7] =~ s/^\s+|\s+$//g;
        $pair[6] =~ s/^\s+|\s+$//g;
        $tem[0] =~ s/^\s+|\s+$//g;
        $sth_extra->bind_param(1,$tem[0]);
        $sth_extra->bind_param(2,$arr[1]);
        $sth_extra->bind_param(3,$_[1]);
        $sth_extra->bind_param(4,$arr[0]);
        $sth_extra->bind_param(5,$arr[2]);
        $sth_extra->bind_param(6,$readcnt{$arr[1]});
        $sth_extra->bind_param(7,$readlength{$arr[1]});
        $sth_extra->bind_param(8,$species_id);
        $sth_extra->bind_param(9,$fqrep{$arr[1]});
        $sth_extra->bind_param(10,$arr[6].';'.$pair[6]);
        $sth_extra->bind_param(11,$arr[7].';'.$pair[7]);
        $sth_extra->bind_param(12,$read_source);
        if ($sth_extra->execute){
          
        }
        else{
          throw("Could not write fastq details to rnaseq extra table");
        }
        #Storing tested reads
        $sth_tested_reads->bind_param(1,$arr[1]);
        $sth_tested_reads->bind_param(2,$species_id);
        
        if ($sth_tested_reads->execute){

        }
        else{
          throw("Could not store tested read details");
        }
        #Processing each individual run stats
        $sth_stats->bind_param(1,$_[1]);
        $sth_stats->bind_param(2,'full');
        #checks when percentage identity is not set, assign 0
        if ($mapped{$run_id} =~ m/N/ || $mapped{$run_id} eq ''){
          $sth_stats->bind_param(3,0);
          $sth_stats->bind_param(4,0);
        }else{
          $sth_stats->bind_param(3,$mapped{$run_id});
          $sth_stats->bind_param(4,$is_paired{$run_id});
        }
        $sth_stats->bind_param(5,$run_id);
        
        if ($sth_stats->execute){
        }
        else{
          throw("Could not write alignment stats to rnaseq data run table");
        }
        
      }
      
    }
    
  }

}

sub write_output{
  my ($self) = @_;
}



1;
