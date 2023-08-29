# Copyright [2018-2020] EMBL-European Bioinformatics Institute
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
# limitations under the License.


package Bio::EnsEMBL::Pipeline::Runnable::TranscriptomicRegistryAdaptor;

use strict;
use warnings;
use feature 'say';
use parent ('Bio::EnsEMBL::DBSQL::DBAdaptor');

sub fetch_all_gca {
  my ($self,$max_version_only) = @_;

  my $sql = "SELECT chain,version FROM assembly order by chain,version";
  my $sth = $self->dbc->prepare($sql);
  $sth->execute();

  my $output_hash;
  while (my ($chain,$version) = $sth->fetchrow_array()) {
    unless($chain =~ /^GCA\_(\d){9}/) {
      next;
    }

    unless($output_hash->{$chain}) {
      $output_hash->{$chain} = [];
    }

    if($max_version_only) {
      if((${$output_hash->{$chain}}[0])) {
        unless(${$output_hash->{$chain}}[0] > $version) {
          ${$output_hash->{$chain}}[0] = $version;
        }
      } else {
        ${$output_hash->{$chain}}[0] = $version;
      }
    } else {
      push(@{$output_hash->{$chain}},$version);
    }
  }

  my $output_array = [];
  foreach my $chain (keys(%{$output_hash})) {
    my $version_array  = $output_hash->{$chain};
    foreach my $version (@{$version_array}) {
      push(@{$output_array},$chain.".".$version);
    }
  }

  return($output_array);
}

sub fetch_samples_readcount_by_taxon_id {
  my ($self,$taxon_id,$type) = @_;
  my @samples;
  my $sql = "SELECT sample_name,total_read_count FROM rnaseq_data_summary WHERE species_id=? and status != 'unusable'";
  my $sth = $self->dbc->prepare($sql);
  $sth->bind_param(1,$taxon_id);
  if ($sth->execute()){
    while (my ($sample_name,$read_count) = $sth->fetchrow_array()){
      unless(($sample_name !~ /^\d+/) && ($read_count =~ /^\d+/)) {
        next;
      }
      push @samples, $sample_name."\t".$read_count;
    }
  }

  else {
    $self->throw("Could not execute query to retrieve samples for taxonomy ".$taxon_id);
  }

  return @samples;
}
sub fetch_short_reads_by_species_id {
  my ($self,$species_id,$type) = @_;
  my @reads;
  my $sql = "SELECT * FROM csv_short_read WHERE taxon_id=?";
  my $sth = $self->dbc->prepare($sql);
  $sth->bind_param(1,$species_id);
  if ($sth->execute()){
    while (my @result = $sth->fetchrow_array()){
      push @reads, join("\t ", @result);
    }
  }

  else {
    $self->throw("Could not execute query to retrieve samples for taxonomy ".$species_id);
  }

  return @reads;
}

sub fetch_long_reads_by_species_id {
  my ($self,$species_id,$type) = @_;
  my @reads;
  my $sql = "SELECT * FROM csv_long_read WHERE taxon_id=?";
  my $sth = $self->dbc->prepare($sql);
  $sth->bind_param(1,$species_id);
  if ($sth->execute()){
    while (my @result = $sth->fetchrow_array()){
      push @reads, join("\t ", @result);
    }
  }

  else {
    $self->throw("Could not execute query to retrieve samples for taxonomy ".$species_id);
  }

  return @reads;
}

sub fetch_sample_name_by_run_id {
  my ($self,$run_id,$type) = @_;
  my @samples;
  my $sql = "SELECT SM FROM short_read_data WHERE ID=? limit 1";
  my $sth = $self->dbc->prepare($sql);
  my $result = '';
  $sth->bind_param(1,$run_id);
  if ($sth->execute()){
    while ($result = $sth->fetchrow_array()){
      unless($result !~ /^\d+/) {
        next;
      }
    }
  }

  else {
    $self->throw("Could not execute query to retrieve sample name for accession ".$run_id);
  }

  return $result;
}

sub fetch_samples_by_taxon_id {
  my ($self,$taxon_id,$type) = @_;
  my @samples;
  my $sql = "SELECT distinct(sample_name) FROM csv_table WHERE species_id=?";
  my $sth = $self->dbc->prepare($sql);
  $sth->bind_param(1,$taxon_id);
  if ($sth->execute()){
    while (my $result = $sth->fetchrow_array()){
      unless($result !~ /^\d+/) {
        next;
      } 
      push @samples, $result;
    }
  }

  else {
    $self->throw("Could not execute query to retrieve samples for taxonomy ".$taxon_id);
  }
  
  return @samples;
}

sub fetch_n50_by_gca {
  my ($self,$chain_version,$type) = @_;

  my ($chain,$version) = $self->split_gca($chain_version);

  my $sql = "SELECT assembly_id FROM assembly WHERE chain=? and version=?";
  my $sth = $self->dbc->prepare($sql);
  $sth->bind_param(1,$chain);
  $sth->bind_param(2,$version);
  $sth->execute();

  my $assembly_id = $sth->fetchrow();
  unless($assembly_id) {
    $self->throw("Could not find assembly id for assembly with chain ".$chain." and version ".$version);
  }

  my $column_type;
  if($type eq 'contig') {
    $column_type = "contig_N50";
  } elsif($type eq 'scaffold') {
    $column_type = "scaffold_N50";
  } else {
    $self->throw("N50 type parameter missing. Expected type to be 'contig' or 'scaffold'");
  }

  $sql = "SELECT ".$column_type." FROM meta WHERE assembly_id=?";
  $sth = $self->dbc->prepare($sql);
  $sth->bind_param(1,$assembly_id);
  $sth->execute();

  my ($n50) = $sth->fetchrow();

  unless($n50) {
    return(0);
  }

  return($n50);
}

sub fetch_gca_by_constraints {
  my ($self,$contig_n50,$scaffold_n50,$total_length,$levels,$max_version_only,$genome_rep) = @_;

  unless($contig_n50) { $contig_n50 = 0;}
  unless($scaffold_n50) { $scaffold_n50 = 0;}
  unless($total_length) { $total_length = 0;}
  unless($levels) { $levels = ['contig','scaffold','chromosome'];}
  unless($genome_rep) { $genome_rep = 'full';}
  my $output_hash;
  foreach my $level (@{$levels}) {
    my $sql;

    $sql = "SELECT chain,version FROM assembly JOIN meta using(assembly_id) WHERE ".
           " contig_N50 >= ? AND total_length >= ? AND assembly_level = ? AND genome_rep = ?";

    unless($level eq 'contig') {
      $sql .= " AND (scaffold_N50 >= ? || scaffold_N50 IS NULL)";
    }

    my $sth = $self->dbc->prepare($sql);
    $sth->bind_param(1,$contig_n50);
    $sth->bind_param(2,$total_length);
    $sth->bind_param(3,$level);
    $sth->bind_param(4,$genome_rep);
    unless($level eq 'contig') {
      $sth->bind_param(5,$scaffold_n50);
    }

    $sth->execute();
    while (my ($chain,$version) = $sth->fetchrow_array()) {
      unless($chain =~ /^GCA\_(\d){9}/) {
        next;
      }

      unless($output_hash->{$chain}) {
        $output_hash->{$chain} = [];
      }

      if($max_version_only) {
        if((${$output_hash->{$chain}}[0])) {
          unless(${$output_hash->{$chain}}[0] > $version) {
            ${$output_hash->{$chain}}[0] = $version;
          }
        } else {
          ${$output_hash->{$chain}}[0] = $version;
        }
      } else {
        push(@{$output_hash->{$chain}},$version);
      }
    } # end while $output_hash
}
  my $output_array = [];
  foreach my $chain (keys(%{$output_hash})) {
    my $version_array  = $output_hash->{$chain};
    foreach my $version (@{$version_array}) {
      push(@{$output_array},$chain.".".$version);
    }
  }

  return($output_array);
}


sub fetch_species_name_by_gca {
  my ($self,$chain_version,$type) = @_;

  my ($chain,$version) = $self->split_gca($chain_version);

  my $sql = "SELECT species_id FROM assembly WHERE chain=? and version=?";
  my $sth = $self->dbc->prepare($sql);
  $sth->bind_param(1,$chain);
  $sth->bind_param(2,$version);
  $sth->execute();

  my $species_id = $sth->fetchrow();
  unless($species_id) {
    $self->throw("Could not find species id for assembly with chain ".$chain." and version ".$version);
  }

  $sql = "SELECT species_name FROM species_space_log WHERE species_id=?";
  $sth = $self->dbc->prepare($sql);
  $sth->bind_param(1,$species_id);
  $sth->execute();

  my ($species_name) = $sth->fetchrow();

  return($species_name);
}


sub fetch_assembly_name_by_gca {
  my ($self,$chain_version,$type) = @_;

  my ($chain,$version) = $self->split_gca($chain_version);

  my $sql = "SELECT assembly_id FROM assembly WHERE chain=? and version=?";
  my $sth = $self->dbc->prepare($sql);
  $sth->bind_param(1,$chain);
  $sth->bind_param(2,$version);
  $sth->execute();

  my $assembly_id = $sth->fetchrow();
  unless($assembly_id) {
    $self->throw("Could not find assembly id for assembly with chain ".$chain." and version ".$version);
  }

  $sql = "SELECT assembly_name FROM meta WHERE assembly_id=?";
  $sth = $self->dbc->prepare($sql);
  $sth->bind_param(1,$assembly_id);
  $sth->execute();

  my ($assembly_name) = $sth->fetchrow();

  return($assembly_name);
}

sub fetch_assembly_path_by_gca {
  my ($self,$chain_version,$type) = @_;

  my ($chain,$version) = $self->split_gca($chain_version);
  my $sql = "SELECT assembly_path FROM assembly WHERE chain=? and version=?";
  my $sth = $self->dbc->prepare($sql);
  $sth->bind_param(1,$chain);
  $sth->bind_param(2,$version);
  $sth->execute();
  my $assembly_path = $sth->fetchrow();
  unless($assembly_path) {
    $self->throw("Could not find assembly path for assembly with chain ".$chain." and version ".$version);
  }

  return($assembly_path);
}

sub fetch_stable_id_prefix_by_gca {
  my ($self,$chain_version,$type) = @_;

  my ($chain,$version) = $self->split_gca($chain_version);

  my $sql = "SELECT species_prefix FROM assembly WHERE chain=? and version=?";
  my $sth = $self->dbc->prepare($sql);
  $sth->bind_param(1,$chain);
  $sth->bind_param(2,$version);
  $sth->execute();
  my $stable_id_prefix = $sth->fetchrow();
  unless($stable_id_prefix) {
    $self->throw("Could not find stable id prefix for assembly with chain ".$chain." and version ".$version);
  }

  return($stable_id_prefix);
}

sub fetch_stable_id_start_by_gca {
  my ($self,$chain_version,$type) = @_;
  my ($chain,$version) = $self->split_gca($chain_version);

  my $sql = "SELECT stable_id_space_id FROM assembly WHERE chain=? and version=?";
  my $sth = $self->dbc->prepare($sql);
  $sth->bind_param(1,$chain);
  $sth->bind_param(2,$version);
  $sth->execute();

  my $stable_id_space = $sth->fetchrow();
  unless($stable_id_space) {
    $self->throw("Could not find stable id space for assembly with chain ".$chain." and version ".$version);
  }

  $sql = "SELECT stable_id_space_start FROM stable_id_space WHERE stable_id_space_id=?";
  $sth = $self->dbc->prepare($sql);
  $sth->bind_param(1,$stable_id_space);
  $sth->execute();

  my ($stable_id_space_start) = $sth->fetchrow();

  return($stable_id_space_start);
}


sub split_gca {
  my ($self,$chain_version) = @_;

  unless($chain_version =~ /^(GCA\_\d{9})\.(\d+)$/) {
    $self->throw("Could not parse versioned GCA. GCA used: ".$chain_version);
  }

  my $chain = $1;
  my $version = $2;

  return($chain,$version);
}


1;
