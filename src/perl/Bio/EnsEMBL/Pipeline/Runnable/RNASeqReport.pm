package Bio::EnsEMBL::Pipeline::Runnable::RNASeqReport;

use strict;
use warnings;
use POSIX;
use List::Util qw( min max );
use feature 'say';
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB'); 


sub write_output{
	 my ($self) = @_;
	#print details of species without species or genus level RNASeq data
	my $report = $self->param('output_dir') . "rnaseq_not_found.csv";
	 open D, (">>$report");
	 print D $self->param('chain'),".",$self->param('version'), "\t", $self->param('sp_name'), "\t", $self->param('sp_id'), "\n";
	 
	
}
1;

