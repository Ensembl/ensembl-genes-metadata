package Bio::EnsEMBL::Pipeline::Runnable::RNASeqReport;

use strict;
use warnings;
use POSIX;
use List::Util qw( min max );
use feature 'say';
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB'); 

sub fetch_input{
	 my ($self) = @_;
	
}

sub run{
	 my ($self) = @_;
	 
	 
}

sub tta{
	my $tax_user = 'ensro';
	my $tax_pass = '';
	my $tax_db = $ENV{TAX_DB};
	my $tax_host = $ENV{GBS1};
	my $tax_port = $ENV{GBP1};
	
	
	my $tax_dba =  Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor->new(
        -user   => $tax_user,
        -pass   => $tax_pass,
        -dbname => $tax_db,
        -host   => $tax_host,
        -port   => $tax_port);
    
    my $subfamily = 0; my $family = 0; 
   
my $node_adaptor = $tax_dba->get_TaxonomyNodeAdaptor();
# get a node for an ID
my $node = $node_adaptor->fetch_by_taxon_id($_[0]);

# Finding ancestors
foreach my $ancestor ( @{ $node_adaptor->fetch_ancestors($node)}){
	if ($ancestor->rank() eq "subfamily"){
	}elsif ($ancestor->rank() eq "family"){
		$family = $ancestor->taxon_id();
		last;
	}
	else{}
	
}
return $family;
}
sub write_output{
	 my ($self) = @_;
	#print details of species without species or genus level RNASeq data
	my $report = $self->param('output_dir') . "rnaseq_not_found.csv";
	 open D, (">>$report");
	 print D $self->param('chain'),".",$self->param('version'), "\t", $self->param('sp_name'), "\t", $self->param('sp_id'), "\n";
	 
	
}
1;

