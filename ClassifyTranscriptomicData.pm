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
package ClassifyTranscriptomicData;

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
	
	$self->param_required('flagstats');
	$self->param_required('csv_file');
	$self->param_required('read_length_file');
	$self->param_required('fastq_report_file');
	$self->param_required('sp_id');
	$self->param_required('ass_id');
	
}

sub run{
  my ($self) = @_;
 
  my $species = $self->param('sp_id');
  my $ass = $self->param('ass_id');
  my $csv = $self->param('csv_file');
  my $flagstats = $self->param('flagstats');
  my $readl = $self->param('read_length_file');
  my $fastq_report = $self->param('fastq_report_file');
  
  #call function to retrieve alignment and stores details
  &store_alignment_details($csv,$flagstats,$readl,$ass,$species,$fastq_report);
 
  #begin classification of data
  
  &classify_rnaseq($species,$ass);
  sub classify_rnaseq{
    my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -user   => $ENV{GBUSER},
        -dbname => $ENV{REG_DB},
        -host   => $ENV{GBS2},
        -port   => $ENV{GBP2},
        -pass   => $ENV{GBPASS},
        -driver => $ENV{GBDRIVER},
    );
    my $status = ""; my %readcount;
    my $gcnt = 0; my $acnt = 0; my $ass_id = ""; my $sum_readcnt = 0;
    my $sth_rnaseq = $dba->dbc->prepare("insert ignore into rnaseq_data_summary1 (species_id, SM, total_read_count, avg_mapped_identity, status, avg_paired_identity) values (?, ?, ?, ? , ?, ?)");
    my $sth_assembly = $dba->dbc->prepare("update assembly set rnaseq_data = ? where species_id = ?");
    my $sth_csv = $dba->dbc->prepare("insert ignore into csv_data1 (species_id, ID, SM, filename, is_paired, read_length, is_mate_1, is_13plus, CN, PL, DS, url, md5sum, data_level) values (?,?,?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
    my $sql = `mysql -h $ENV{GBS2} -P $ENV{GBP2} -u $ENV{GBUSER} -p$ENV{GBPASS} $ENV{REG_DB} -e \"select sum(read_count) as number_of_reads, SM, round(sum(perc_mapped)/count(*),0) as avg_perc_mapped, rnaseq_data_extra1.data_quality, rnaseq_data_extra1.assembly_id, round(sum(calc_read_length)/count(*),0) as avg_read_length, round(sum(perc_paired)/count(*),0) as avg_perc_paired from rnaseq_data_extra1, rnaseq_data_run_stats1 where project_id in (select project_id from rnaseq_info1 where species_id = $_[0]) and rnaseq_data_extra1.ID = rnaseq_data_run_stats1.ID and rnaseq_data_extra1.species_id = $_[0] and rnaseq_data_extra1.data_quality != 'fail' group by SM\"`;
   
    my @row = split (/\n/, $sql);
    for (my $y = 1; $y < @row; $y+=1){
      my @entry = split(/\t/,$row[$y]);
      $ass_id = $_[1];
      #check how many samples pass the good, weak and usable tests by looking at readcount, avg perc mapped, avg perc paired, average readlength and quality of read
        #code for testing
     if ($entry[0] > 50000000 && $entry[2] > 70 && $entry[3] eq 'pass' && $entry[5] >= 90){
        say "Sample ", $entry[1], " is good";
        $status = "good";
        $readcount{$entry[1]} = $entry[0];#adding read count by sample name
        $gcnt++;
        my $sql = `mysql -h $ENV{GBS2} -P $ENV{GBP2} -u $ENV{GBUSER} -p$ENV{GBPASS} $ENV{REG_DB} -e \"select ID,is_paired,calc_read_length,is_mate_1,CN,PL,DS,url,md5sum,data_level from rnaseq_data_extra1 where species_id = $_[0] and SM = '$entry[1]'\"`;
        my @record = split (/\n/, $sql);
        
        for (my $n = 1; $n < @record; $n+=1){
        
          my @line = split(/\t/,$record[$n]); 
          my $rcnt = 0;
          if ($line[1] == 1){
            #deal with is_paired eneded reads
            my $read1 = $line[0] . '_1.fastq.gz';
            my $read2 = $line[0] . '_2.fastq.gz';
            
            $line[6] =~ s/\\n$//;
            #performing multi row insertion 
            my @query = ([$_[0],$line[0],$entry[1],$read1,$line[1],$line[2],$line[3],0,$line[4],$line[5],$line[6],$line[7],$line[8],$line[9]],[$_[0],$line[0],$entry[1],$read2,$line[1],$line[2],$line[3],0,$line[4],$line[5],$line[6],$line[7],$line[8],$line[9]]);
           foreach my $rec ( @query ) { 
             if ($sth_csv->execute( @$rec )){
             };
            }
            
          }
          else{
            #here we deal with long reads
            $sth_csv->bind_param(1,$_[0]);
            $sth_csv->bind_param(2,$record[0]);
            $sth_csv->bind_param(3,$entry[1]);
            $sth_csv->bind_param(4,$record[0] . '.fastq.gz');
            $sth_csv->bind_param(5,$record[1]);
            $sth_csv->bind_param(6,$record[2]);
            $sth_csv->bind_param(7,$record[3]);
            $sth_csv->bind_param(8,0);
            $sth_csv->bind_param(9,$record[4]);
            $sth_csv->bind_param(10,$record[5]);
            $sth_csv->bind_param(11,$record[6]);
            $sth_csv->bind_param(12,$record[7]);
            $sth_csv->bind_param(13,$record[8]);
            $sth_csv->bind_param(14,$record[9]);
          }
            
        }
      }
      #  #code for testing
      elsif (($entry[0] > 30000000 && $entry[2] > 50 && $entry[3] eq 'pass' && $entry[5] >= 70) || ($entry[0] > 30000000 && $entry[2] > 50 && $entry[3] eq 'warn' && $entry[5] >= 70)){
        say "Sample ", $entry[1], " is weak";
        $status = "weak";
        $acnt++;
        $readcount{$entry[1]} = $entry[0];#adding read count by sample name
        my $sql = `mysql -h $ENV{GBS2} -P $ENV{GBP2} -u $ENV{GBUSER} -p$ENV{GBPASS} $ENV{REG_DB} -e \"select ID,is_paired,calc_read_length,is_mate_1,CN,PL,DS,url,md5sum,data_level from rnaseq_data_extra1 where species_id = $_[0] and SM = '$entry[1]'\"`;
        my @record = split (/\n/, $sql);
        say "Record is ",@record;
        for (my $n = 1; $n < @record; $n+=1){
          my @line = split(/\t/,$record[$n]); 
          my $rcnt = 0;
          my @url = split(/;/,$line[7]);
          my @md5 = split(/;/,$line[8]);
          if ($line[1] == 1){
            #deal with is_paired eneded reads
            my $read1 = $line[0] . '_1.fastq.gz';
            my $read2 = $line[0] . '_2.fastq.gz';
            $line[6] =~ s/\\n$//;
            my @query = ([$_[0],$line[0],$entry[1],$read1,$line[1],$line[2],$line[3],0,$line[4],$line[5],$line[6],$url[0],$md5[0],$line[9]],[$_[0],$line[0],$entry[1],$read2,$line[1],$line[2],$line[3],0,$line[4],$line[5],$line[6],$url[1],$md5[1],$line[9]]);
           foreach my $rec ( @query ) { 
             $sth_csv->execute( @$rec );
            }
            
          }
          else{
            #here we deal with long ended reads
            $sth_csv->bind_param(1,$_[0]);
            $sth_csv->bind_param(2,$record[0]);
            $sth_csv->bind_param(3,$entry[1]);
            $sth_csv->bind_param(4,$record[0] . '.fastq.gz');
            $sth_csv->bind_param(5,$record[1]);
            $sth_csv->bind_param(6,$record[2]);
            $sth_csv->bind_param(7,$record[3]);
            $sth_csv->bind_param(8,0);
            $sth_csv->bind_param(9,$record[4]);
            $sth_csv->bind_param(10,$record[5]);
            $sth_csv->bind_param(11,$record[6]);
            $sth_csv->bind_param(12,$record[7]);
            $sth_csv->bind_param(13,$record[8]);
          }
            
        }
      }
      else{
        
        say "Sample ", $entry[1], " is unusable $entry[0] $entry[2] $entry[3] $entry[5]";
        $status = "unusable";
      }
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
        
      }
      else{#otherwise assign red
        $sth_assembly->bind_param(1,'red');
        $sth_assembly->bind_param(2,$_[0]);
        
      }
    }
    else{}
    if ($sth_assembly->execute){
      say "Rnaseq classification successfully entered";	  
    }
    else{
      throw("Could not write to assembly table");
    }
  }


  sub store_alignment_details{
  	my ($self) = @_;
    my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -user   => $ENV{GBUSER},
        -dbname => $ENV{REG_DB},
        -host   => $ENV{GBS2},
        -port   => $ENV{GBP2},
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
    my $taxon_node = $node_adaptor->fetch_by_taxon_id($_[4]);
    my $species_rank = $taxon_node->rank();
    my %runs;
    my %mapped;
    my %is_paired;
    my %calc;
    my %proj;
    my %readcnt;
    my %fqrep;
    my %rnaseq_info1;
    my $species_id = $_[4];
    open CSV, ("<$_[0]");#csv file
    open FLG, ("<$_[1]");#flagstats
    open LG, ("<$_[2]");#readlength
    open FQREP, ("<$_[5]") or $self->complete_early("No fastq report for this species");#readlength
   
    #not taking readcount from rnaseq meta file because entry is not consistent
    #reading fastq read length file and storing all entries in a hash
    while (<LG>){
      my @val = split(/\t/,$_);
      if ($val[0] =~ /_/){
        my @id = split(/_/, $val[0]);
        $readcnt{$id[0]} = $val[1];
      }
      else{
        my @id = split(/\./, $val[0]);
        $readcnt{$id[0]} = $val[1];
      }
    }
    #reading flagstat file and storing all entries in a hash
    while (<FLG>){
        my @values = split(/\t/,$_);
        $mapped{$values[0]} = $values[1];
        $is_paired{$values[0]} = $values[2];
    }
    #reading fastq report file and storing all entries in a hash
    while (<FQREP>){
        my @values = split(/\t/,$_);
        #check if fastq is is_paired ended
        if ($values[0] =~ m/_/){
          my @d = split(/_/, $values[0]);
          $fqrep{$d[0]} = lc($values[1]);
          say "Fastq report is ",$fqrep{$d[0]};
        }
        else{
          my @d = split(/\./, $values[0]);
          $fqrep{$d[0]} = lc($values[1]);
          say "Fastq report is ",$fqrep{$d[0]};
        }
        
        
    }
       my $sth_info = $dba->dbc->prepare("select project_id,species_id from rnaseq_info1");
       my $sth_rnaseq = $dba->dbc->prepare("insert into rnaseq_info1 (species_id, project_id, data_level) values (?, ?, ?)");
       my $sth_extra = $dba->dbc->prepare("insert ignore into rnaseq_data_extra1 values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
       my $sth_stats = $dba->dbc->prepare("insert ignore into rnaseq_data_run_stats1 values (?, ?, ?, ?, ?)");
       if ($sth_info->execute){
          while (my @entry = $sth_info->fetchrow_array()){
            $rnaseq_info1{$entry[0] . '_' . $entry[1]} = $entry[1];
          }
        }
    my %csv = {};
    #reading csv file and processing all entries
    while (<CSV>){
      say "analysing line $_";
      $_ =~ s/[^[:ascii:]|[:alnum:]]//g;
      say "New line is now $_";
      my @line = split(/\t/,$_);
      $csv{$line[4]} = $_;
    }
    foreach my $read (keys %csv){
      say "read is $read";
      my @line = split(/\t/,$csv{$read});
      say "Length is ",scalar(@line);
      my @tem = split(/,/,$line[10]);
      say "species id is ", length($line[3]) . " " . $line[3];
      if ($read =~ /_1.fastq.gz/){
        (my $mate = $read) =~ s/_1/_2/;
        my @pair = split(/\t/,$csv{$mate});
        chomp($line[12]);chomp($pair[12]);
        $line[11] = $line[11] . ';'. $pair[11];
        $line[12] = $line[12] . ';' .$pair[12];
      }
      else{next;}
      if (!exists $runs{$line[1]}){
        my $read_source = "";
        say "fields required are; ", $line[0], "\t", $line[1], "\t", $line[2], "\t", $line[8], "\t", $tem[0], "\t", $line[9];
        $runs{$line[1]} = 1;
        $calc{$line[1]} = $line[6];
        if ($readcnt{$line[1]} !~ /\d+/){
          say "run does not exist";
          next;
        }
        my $proj_spec_id = $tem[0] . '_' . $_[4];
        if ((exists($rnaseq_info1{$proj_spec_id}) && ($rnaseq_info1{$proj_spec_id} !~ m/$_[4]/)) or (!exists($rnaseq_info1{$proj_spec_id}))){
          $sth_rnaseq->bind_param(1,$_[4]);
          say "Species tested is $_[4]";
          $sth_rnaseq->bind_param(2,$tem[0]);
        }
        $taxon_node = $node_adaptor->fetch_by_taxon_id($line[3]);
        my $source_rank = $taxon_node->rank();
        if ($species_id != $line[3]){#where species id differs with entry in csv file, check that source is subspecies.
          if ($source_rank eq 'subspecies'){
            $sth_rnaseq->bind_param(3,'subspecies');
            $read_source = 'subspecies';
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
        say "Project tested is $tem[0]";
        if (!exists $proj{$tem[0]}){
           say "Species tested is $_[4]";
           if ($sth_rnaseq->execute){
            
          }
          else{
            throw("Could not write mapped data to rnaseq info table");
          }
          $proj{$tem[0]} = 1;
        }
        say "Description before insertion is $line[9]";
        $line[10] = substr($line[10],0,250);
        $line[10] =~ tr/:\t/ /;
        $line[0] = substr($line[0],0,50);
        my $is_mate_1 = 0;
        chomp($fqrep{$line[1]});
        my $regex = '\S+_(\d)\.\S+';
        my ($pair) = $line[4] =~ /$regex/;
            if (($pair == 1) or ($line[2] == 0)) {
                $is_mate_1 = 1;
            }
            else {
                $is_mate_1 = 0;
            }
        $sth_extra->bind_param(1,$tem[0]);
        $sth_extra->bind_param(2,$line[1]);
        $sth_extra->bind_param(3,$_[3]);
        $sth_extra->bind_param(4,$line[0]);
        $sth_extra->bind_param(5,$line[2]);
        $sth_extra->bind_param(6,$is_mate_1);
        $sth_extra->bind_param(7,$line[8]);
        $sth_extra->bind_param(8,$readcnt{$line[1]});
        $sth_extra->bind_param(9,$line[6]);
        $sth_extra->bind_param(10,$calc{$line[1]});
        $sth_extra->bind_param(11,$line[9]);
        $sth_extra->bind_param(12,$line[10]);
        $sth_extra->bind_param(13,$species_id);
        $sth_extra->bind_param(14,$fqrep{$line[1]});
        $sth_extra->bind_param(15,$line[11]);
        $sth_extra->bind_param(16,$line[12]);
        $sth_extra->bind_param(17,$read_source);
        if ($sth_extra->execute){
          
        }
        else{
          throw("Could not write mapped to rnaseq extra table");
        }
        $sth_stats->bind_param(1,$_[3]);
        $sth_stats->bind_param(2,'test');
        #checks when percentage identity is not set, assign 0
        if ($mapped{$line[1]} =~ m/N/ || $mapped{$line[1]} eq ''){
          $sth_stats->bind_param(3,0);
          $sth_stats->bind_param(4,0);
        }else{
          $sth_stats->bind_param(3,$mapped{$line[1]});
          $sth_stats->bind_param(4,$is_paired{$line[1]});
        }
        $sth_stats->bind_param(5,$line[1]);
        
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
