=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

AssemblyLoading

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Pipeline::Runnable::AssemblyLoading;

use strict;
use warnings;

use Net::FTP;
use Time::Piece;
use File::Fetch;
use File::Temp;
use File::Spec::Functions qw(catfile splitpath catdir);
use File::Path qw(make_path);
use Digest::MD5;
use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError) ;

use Bio::EnsEMBL::Hive::Utils qw(destringify);

use Bio::EnsEMBL::IO::Parser::Fasta;

use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Attribute;
use feature 'say';
use Data::Dumper;
use Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 param_defaults

 Arg [1]    : None
 Description: Default parameters
               toplevel_as_sequence_levels => 1,
               _report_name => '#assembly_accession#_#assembly_name#_assembly_report.txt',
               _md5checksum_name => 'md5checksums.txt',
               _genome_file_name => '#assembly_accession#_#assembly_name#_genomic.fna',
               _genome_zip_ext => '.gz',
               _molecule_relations => {
                 'assembled-molecule' => 'chromosome',
                 'unplaced-scaffold' => 'scaffold',
                 'unlocalized-scaffold' => 'scaffold',
                 'alt-scaffold' => 'scaffold',
               },
               ftp_user => 'anonymous',
               ftp_pass => undef,
               load_non_nuclear => 0,
               _agp_branch => 3,
               _coord_systems => {},
 Returntype : Hashref
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::param_defaults},
    
    # if any DNA slice sequence length is greater than _MAX_SLICE_LENGTH base pairs,
    # all DNA slice sequences will be cut into _SUBSLICE_LENGTH-base-pair long sequences
    _MAX_SLICE_LENGTH => 950000000,  # 950 million base pairs
    _SUBSLICE_LENGTH => 10000000,    # 10 million base pairs
    _exceeded_max_slice_length => 0, # it will be set to 1 when any DNA slice sequence length is greater than _MAX_SLICE_LENGTH
    
    toplevel_as_sequence_levels => 1,
    _report_name => '#assembly_accession#_#assembly_name#_assembly_report.txt',
    _md5checksum_name => 'md5checksums.txt',
    _genome_file_name => '#assembly_accession#_#assembly_name#_genomic.fna',
    _genome_zip_ext => '.gz',
    _molecule_relations => {
      'assembled-molecule' => 'chromosome',
      'unplaced-scaffold' => 'scaffold',
      'unlocalized-scaffold' => 'scaffold',
      'alt-scaffold' => 'scaffold',
    },
    ftp_user => 'anonymous',
    ftp_pass => undef,
    load_non_nuclear => 0,
    _agp_branch => 3,
    _coord_systems => {},
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Retrieve external db ids for INSDC, RefSeq and UCSC. Retrieve the
              attribute information for 'toplevel' and 'karyotype_rank'.
              Download the report file using 'full_ftp_path' and '_report_name'
 Returntype : None
 Exceptions : Throws if 'target_db' is not set
              Throws if 'assembly_accession' is not set
              Throws if 'assembly_name' is not set
              Throws if 'full_ftp_path' is not set

=cut

sub fetch_input {
  my ($self) = @_;

  my $assembly_accession = $self->param_required('assembly_accession');
  my $assembly_name = $self->param_required('assembly_name');
  my $species_name = $self->param_required('species_name');
  
  my $report_dir;
  if ($self->param_is_defined('genome_path')) {
    $report_dir = $self->param('genome_path') . "$species_name/$assembly_accession/";
    $self->param('genome_path', $report_dir);
  }
  else {
    $report_dir = File::Temp->newdir;
    $self->param('genome_path', $report_dir);
  }
  if (!-d $report_dir) {
    make_path($report_dir);
  }
  if (-e $report_dir){
    #`rm -r $report_dir`;
  }
  my $query = "cd $report_dir; /nfs/production/flicek/ensembl/genebuild/do1/datasets download genome accession " . $assembly_accession;
 
  if(system($query)){
  #  `rm -r $report_dir`;
    say "Download failed";
  }
  else{
    $self->param('genome_path', $report_dir);
    say "Download succeeded for genome";
  }
  my $genome_extract = "cd $report_dir; unzip " . $report_dir . 'ncbi_dataset.zip';
  say "Command to extract is $genome_extract";
  if(system($genome_extract)){
    say "Extraction failed";
  }
  else{
    say "Extraction completed successfully";
  }
  # Rename downloaded genome
  my $rename_cmd = "mv " . $report_dir . '/ncbi_dataset/data/' . $assembly_accession . '/*.fna ' . $report_dir . '/' . $assembly_accession . '.fasta'; 
  if(system($rename_cmd)){
    say "File renaming failed";
  }
  else{
    say "File renaming completed successfully";
  } 		
}


=head2 run

 Arg [1]    : None
 Description: 
              Decompress genome and perform indexing
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;
  my $file = "";
  my $genome = $self->param('genome_path'); my $speciesn = $self->param('species_name');

  #index genome
  my $output = $self->param('genome_path').$self->param('assembly_accession').".fasta";
  say "now indexing with STAR";
  my $query = $self->param('star_path') . ' --limitGenomeGenerateRAM 504273583712  --runThreadN 12 --runMode genomeGenerate --outFileNamePrefix ' . $self->param('genome_path') . ' --genomeDir '. $self->param('genome_path') . ' --genomeFastaFiles ' . $output;
  if(system($query)){
    $self->throw("Error indexing genome via Star\nError code: $?\n");
  }
  say "finished creating Star index";
  my $genome_index = $output . '.mmi';
  my $minimap_index = $self->param('minimap2_path'). " -d $genome_index $output";
  say "command to run is $minimap_index";
 if(system($minimap_index)) {
    $self->throw("Error indexing genome via minimap2\nError code: $?\n");
  }  
}


sub write_output {
  my ($self) = @_;

}

1;
