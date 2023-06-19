#!/usr/bin/env perl

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use feature 'say';
use Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor;
use Getopt::Long;
use Data::Dumper;

#Initialising standard variables to hold assembly meta info
my %clade; 
my %contig;
my %version;
my %accession;
my %final;
my $final_cnt = 0;
my $initial_cnt = 0;
my $input = "";
GetOptions('file:s' => \$input,);

my $candidate_assembly = $input;
my $genus_id = "";
my ($self) = @_;
#setting db adaptor for registry db

my $taxonomy_adaptor = Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor->new(
						-user   => $ENV{GBUSER_R},
						-pass   => '',
						-dbname => 'ncbi_taxonomy_109',
						-host   => 'mysql-ens-mirror-1',
						-port   => 4240
					 );
my $registry_adaptor = new Bio::EnsEMBL::DBSQL::DBAdaptor(
	-user   => $ENV{GBUSER},
	-dbname => $ENV{REG_DB},
	-host   => $ENV{GBS1},
	-port   => $ENV{GBP1},
	-pass   => $ENV{GBPASS},
	-driver => $ENV{GBDRIVER},
);

#Process all non vert assemblies in registry per species
#my $sth = $registry_adaptor->dbc->prepare("select taxonomy, assembly.assembly_id, chain, version, clade, contig_N50, assembly_level, assembly_name, replace(trim(lcase(species_name)), ' ', '_') as species from assembly join meta using (assembly_id)  where genome_rep = ?  and contig_N50 > ? and (rnaseq_data in (?,?,?,?)) and round((30*total_length)/100) > total_gap_length and is_current = 1 and clade not in ('reptiles','mammalia','teleostei','rodentia','teleostei','teleostei','humans','marsupials','teleostei','vertebrates','teleostei','non_vertebrates')");
my $sth = $registry_adaptor->dbc->prepare("select taxonomy, assembly.assembly_id, chain, version, clade, contig_N50, assembly_level, assembly_name, replace(trim(lcase(species_name)), ' ', '_') as species from assembly join meta using (assembly_id)  where genome_rep = ?  and contig_N50 > ? and (rnaseq_data in (?,?,?,?)) and round((30*total_length)/100) > total_gap_length and is_current = 1 and clade = 'teleostei'");
#setting filters for candidate assemblies
$sth->bind_param(1,'full');
$sth->bind_param(2,100000);
$sth->bind_param(3,'available');
$sth->bind_param(4,'red');
$sth->bind_param(5,'amber');
$sth->bind_param(6,'green');
if ($sth->execute){
	#looping through candidate assemblies
  while (my @result = $sth->fetchrow_array()){
	  #remove none standard characters from species name
    $result[8] =~ s/\s+-\s+\w+:\w+$//;
    $result[8] =~ s/[[:space:][:punct:]]+/_/g;
    $result[8] =~ s/_{2,}/_/g;
    $result[8] =~ s/^_|_$//g;
    $result[7] =~ s/\s+-\s+\w+:\w+$//;
    #remove onn standard characters from assembly name
    $result[7] =~ s/[[:space:]]+/_/g;
    $result[7] =~ s/_{2,}/_/g;
    $result[7] =~ s/\+/_/g;
    #get genus id to be used for meta data download
    my $node_adaptor = $taxonomy_adaptor->get_TaxonomyNodeAdaptor();
    my $taxon_node = $node_adaptor->fetch_by_taxon_id($result[0]);
    foreach my $ancestor ( @{ $node_adaptor->fetch_ancestors($taxon_node)}){
      if ($ancestor->rank eq 'species'){
	say "species";
      }
      elsif ($ancestor->rank eq 'genus'){
	$genus_id = $ancestor->taxon_id();
	last;
      }
      else{}
    }
    $final{$result[1]} = join("\t", @result). "\t$genus_id";
  }	
}
else{
	throw("Could not execute query against the registry db");
}
#write filtered candidate assemblies to file
my %unique_list;
my $u = 0;
(my $ulist = $candidate_assembly) =~ s/species.csv/nr_species.csv/;
open FINAL, (">>$candidate_assembly");
#open RUN_LIST, (">>$ulist");
foreach my $assembly (values %final){
	print FINAL $assembly, "\n";
        my @unique_species = split(/\t/,$assembly);
        $unique_list{$unique_species[0]} = $assembly;
	$final_cnt++;
}
say "Unique list is $ulist";
open LIST, (">>$ulist");
foreach my $run_list (values %unique_list){
   print LIST $run_list, "\n";
   $u++;
}
say "Initial list of qualifying assemblies = $final_cnt";say "Unique list for fastq download = $u";

