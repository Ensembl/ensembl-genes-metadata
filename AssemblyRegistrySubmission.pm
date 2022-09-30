=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveStar2Introns

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package AssemblyRegistrySubmission;

use warnings;
use strict;
use feature 'say';

use Net::FTP;
use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor;
use POSIX qw(strftime);
use List::MoreUtils qw(any);
use LWP::UserAgent;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

=head2 param_defaults

 Arg [1]    : None
 Description: Default parameters for the analysis
               ena_base_url => 'http://www.ebi.ac.uk/ena/portal/api/search?display=report',
               files_domain => 'domain=read&result=read_run',
               files_fields => 'run_accession,study_accession,experiment_accession,sample_accession,secondary_sample_accession,instrument_platform,instrument_model,library_layout,library_strategy,read_count,base_count,fastq_ftp,fastq_aspera,fastq_md5,library_source,library_selection,center_name,study_alias,experiment_alias,experiment_title,study_title',
               sample_domain => 'domain=sample&result=sample',
               sample_fields => 'accession,secondary_sample_accession,bio_material,cell_line,cell_type,collected_by,collection_date,country,cultivar,culture_collection,description,dev_stage,ecotype,environmental_sample,first_public,germline,identified_by,isolate,isolation_source,location,mating_type,serotype,serovar,sex,submitted_sex,specimen_voucher,strain,sub_species,sub_strain,tissue_lib,tissue_type,variety,tax_id,scientific_name,sample_alias,center_name,protocol_label,project_name,investigation_type,experimental_factor,sample_collection,sequencing_method',
               download_method => 'ftp',
               separator => '\t',
               _read_length => 1,
               _centre_name => 'ENA',
               print_all_info => 0,
               paired_end_only => 1, #by default, module will only add paired-end data to the csv, add "paired_end_only => 0" to pipeline config to include single end data
               read_type => 'short_read',
               instrument_platform => 'ILLUMINA',
               taxon_id_restriction => 0, # Set it to 1 if you are using 'study_accession' which has multiple species
               max_long_read_read_count => 1000000, # if a long read set has more than 1,000,000 reads, discard it. Very crude filter on rax reads
 Returntype : Hashref
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    ena_base_url => 'https://www.ebi.ac.uk/ena/portal/api/search?display=report',#'http://www.ebi.ac.uk/ena/data/warehouse/search?display=report',
    files_domain => 'domain=read&result=read_run',
    files_fields => 'run_accession,study_accession,experiment_accession,sample_accession,secondary_sample_accession,instrument_platform,instrument_model,library_layout,library_strategy,nominal_length,read_count,base_count,fastq_ftp,fastq_aspera,submitted_ftp,fastq_md5,library_source,library_selection,center_name,study_alias,experiment_alias,experiment_title,study_title',
    sample_domain => 'domain=sample&result=sample',
    sample_fields => 'accession,secondary_sample_accession,bio_material,cell_line,cell_type,collected_by,collection_date,country,cultivar,culture_collection,description,dev_stage,ecotype,environmental_sample,first_public,germline,identified_by,isolate,isolation_source,location,mating_type,serotype,serovar,sex,submitted_sex,specimen_voucher,strain,sub_species,sub_strain,tissue_lib,tissue_type,variety,tax_id,scientific_name,sample_alias,center_name,protocol_label,project_name,investigation_type,experimental_factor,sample_collection,sequencing_method',
    download_method => 'ftp',
    separator => '\t',
    _read_length => 1,
    _centre_name => 'ENA',
    print_all_info => 0,
    paired_end_only => 1, #by default, module will only add paired-end data to the csv, add "paired_end_only => 0" to pipeline config to include single end data 
    ftp_base_dir => '/genomes/all/',
    ftphost => "ftp.ncbi.nlm.nih.gov",
    ftpuser => "anonymous",
    ftppassword => "",
  }
}

=head2 _populate_query

 Arg [1]    : String or Arrayref of Strings, it should be a study accession or a taxon id or an array of them
 Arg [2]    : String $format, the format string for the param
 Description: It populates the 'query' parameters with an arrayref of string base on the format Arg[2] and the
              parameter(s) in Arg[1]. You cannot mix taxon_ids and study_accessions.
 Returntype : None
 Exceptions : None

=cut

sub _populate_query {
  my ($self, $params, $format) = @_;

  if (ref($params) eq 'ARRAY') {
    my @queries;
    foreach my $param (@$params) {
      push(@queries, sprintf($format, $param));
    }
    $self->param('query', \@queries);
  }
  else {
    $self->param('query', [sprintf($format, $params)]);

  }

}

=head2 fetch_input

 Arg [1]    : None
 Description: Check that file containing genome accessions in list format and other parameters are parsed as arguments to module.
 Returntype : None
 Exceptions : Throws if no file is parsed or exists

=cut

sub fetch_input {
  my ($self) = @_;

  unless($self->param('config_file')){
    $self->complete_early("Required: Config file containing assembly accessions in an array [GCA_XXX.X,GCA_XXX.X] e.g -config_file xxx");
  }
  unless(-e ($self->param('config_file'))){
    $self->complete_early("Required: Please check that file specified exists");
  }
}

=head2 run

 Arg [1]    : None
 Description: Reads list of accessions and compares with current registry entries. 
              Retrieves assembly meta data for all non-existing assemblies from the archives (i.e., NCBI or ENA)
 Returntype : None
 Exceptions : Throws if cannot open or close the files
              Throws if the accession is not found in the CSV table 'table_name'

=cut

sub run {
  my ($self) = @_;

  my %existing_species_prefix;#hash to store existing species prefixes
  my %existing_species_id;#hash to store existing species id
  my %existing_assemblies; #hash to store existing GCA's
  my $ftp_base_dir = '/genomes/all/';
  my $ftphost = "ftp.ncbi.nlm.nih.gov";
  my $ftpuser = "anonymous";
  my $ftppassword = "";
  my $general_hash;
  my $search_term = "";

  open(IN,$self->param('config_file')) || die("Could not open $self->param('config_file')");
  while (<IN>){
    my $line = $_;

    $line =~ s/\s//g;
    if($line =~ /^(.+)\=(.+)/) {
      my $key = $1;
      my $value = $2;
      say "Found key/value pair: ".$key." => ".$value;

      # Note that in the ini file the key for write user is user_w for clarity, but in the actual
      # hive config it's user (since hive expects a user key). This is just a substitution to account
      # for this
      if($key eq 'user_w') {
        $key = 'user';
        $general_hash->{$key} = $value;
      }
      if($key eq 'user_r') {
        $general_hash->{$key} = $value;
      }
      if($key eq 'password') {
        $general_hash->{$key} = $value;
      }
      if($key eq 'registry_db_host') {
        $general_hash->{$key} = $value;
      }
      if($key eq 'registry_db_port') {
        $general_hash->{$key} = $value;
      }
      if($key eq 'taxonomy_db_host') {
        $general_hash->{$key} = $value;
      }
      if($key eq 'taxonomy_db_port') {
        $general_hash->{$key} = $value;
      }
      if($key eq 'import_type') {
        $general_hash->{$key} = $value;
      }
      else{
        $general_hash->{$key} = $value;
      }
    }elsif($line eq "\n") {
        # Skip
    }else {
      say "Line format not recognised. Skipping line:\n".$line;
    }
  }
  close IN || die("Could not close ". $self->param('config_file'));
  unless($general_hash->{'output_path'}) {
    die("Could not find an output path setting in the config. Expected setting".
        "output_path=/path/to/output/dir/");
  }
  my $output_path = $general_hash->{'output_path'};
  unless(-e $output_path) {
    system("mkdir -p ".$general_hash->{'output_path'});
  }
  #Attempting to register each assembly
  #Create registry adaptor
  my $assembly_registry = new Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor(
    -host    => $general_hash->{'registry_db_host'},
    -port    => $general_hash->{'registry_db_port'},
    -user    => $general_hash->{'user'},
    -pass    => $general_hash->{'password'},
    -dbname  => $general_hash->{'registry_dbname'}
  );
  
  #Get record of all existing assemblies
  my $sql = "SELECT species_prefix, concat(chain,'.', version), taxonomy, annotated_status FROM assembly order by chain,version";
  my $sth = $assembly_registry->dbc->prepare($sql);
  $sth->execute();
  while (my ($species_prefix,$accession,$taxonomy,$annotated_status) = $sth->fetchrow_array()) {
    unless($accession =~ /^GCA\_(\d){9}/) {
      next;
    }
    $search_term = $accession . '_' . $annotated_status;
    $existing_assemblies{$search_term} = 1;
    $existing_species_prefix{$taxonomy} = $species_prefix;
    $existing_species_id{$taxonomy} = 1;
  }
 
  #Verify new assemblies to register have not already been registered
  my $assembly_accessions = $general_hash->{'assembly_accessions'};
  $assembly_accessions =~ /\[(.+)\]/;
  my $accession_string = $1;
  $accession_string =~ s/ //g;
  my @accession_array = split(',',$accession_string);
  unless(scalar(@accession_array)) {
    die("Issue parsing assembly_accessions line. Format expected:\n".
        "assembly_accessions=[GCA_000952055.2,GCA_000164805.2,GCA_001604975.1]");
  }
  foreach my $accession (@accession_array) {
    my $assembly_hash = {};
    chomp($accession);
     my $species_prefix = "";

    say "Processing accession: ".$accession;
    unless($accession =~ /GCA_([\d]{3})([\d]{3})([\d]{3})\.\d+/) {
      die("Found an assembly accession that did not match the regex. Offending accession: ".$accession);
    }
    my $assembly_ftp_path = $ftp_base_dir.'GCA/'.$1.'/'.$2.'/'.$3.'/';
    my $full_assembly_path;
    my $assembly_name;

    my $ftp = Net::FTP->new($ftphost) or die("Can't open $ftphost");
    $ftp->login($ftpuser, $ftppassword) or die("Can't log $ftpuser in");
    $ftp->cwd($assembly_ftp_path);
    my @ftp_dir_contents = $ftp->ls;
    foreach my $entry (@ftp_dir_contents) {
      if($entry =~ /^$accession\_(.+)$/) {
        $full_assembly_path = $assembly_ftp_path.$&."/";
        $assembly_name = $1;
      }
    }
    unless($full_assembly_path && $assembly_name) {
      die("Issue finding ftp path for the following GCA: $accession $assembly_ftp_path");
    }
    #Retrieve metadata related info from assembly report file
    my @assembly_report = parse_assembly_report($ftp,$accession,$assembly_name,$full_assembly_path);
    my ($taxon_id,$assembly_level,$refseq_accession,$wgs_id,$bioproject,$genome_rep,$total_length,$scaffold_cnt,$contig_cnt,$contig_N50,$scaffold_N50,$total_gap_length,$top_level_cnt,$assembly_type,$date,$assembly_ftp,$submitter,$common_name) = @assembly_report;
    #Retrieve taxon id related info from taxonomy db
    my ($parent_id,$clade,$subspecies_name,$species_common_name,$species_name) = $self->get_clade($taxon_id,$output_path);
    #Prioritise common name from taxonomy db over assembly report file entry
    if ((!defined($species_common_name)) && (defined($common_name))){
      $species_common_name = $common_name;
    }
    #Fetch transcriptomic data status from registry
    my $rnaseq_status = "";
    my $rna_cnt = 0;
    my $sql = "SELECT DISTINCT(rnaseq_data) from assembly where species_id = $parent_id";
    my $sth = $assembly_registry->dbc->prepare($sql);
    $sth->execute();
    $rnaseq_status = $sth->fetchrow_array();
    
    if (!defined($rnaseq_status)){
      $rnaseq_status = 'not available';
    }
    if ($rnaseq_status eq 'not available'){#if no transcriptomic data exists, check for new data
      $self->_populate_query($taxon_id, 'tax_tree(%s) AND instrument_platform=ILLUMINA AND library_source=TRANSCRIPTOMIC');
      foreach my $query (@{$self->param('query')}) {
        my $ua = LWP::UserAgent->new;
        $ua->env_proxy;
        my $url = join('&', $self->param('ena_base_url'), 'query="'.$query.'"', $self->param('files_domain'), 'fields='.$self->param('files_fields'));
        $self->say_with_header($url);
        my $response = $ua->get($url);
        my $fastq_file = 'fastq_'.$self->param('download_method');
        if ($response->is_success) {
          my $content = $response->decoded_content(ref => 1);
          if ($content) {
            while ($$content =~ /^(\w+.*)$/mgc) {
              my $line = $1;
              if ($line =~ /^[a-z]/) {
              }
              else {#filtering unwanted read types
                next if ($line =~ / infected | [iIu]mmune| challenge |tomi[zs]ed/);
                next if ($line =~ /[Mm]is\w{0,2}RNA|lncRNA|circRNA|small RNA/);
                $rna_cnt++;
              }
            }
            if ($rna_cnt > 0){
              $rnaseq_status = "available";
            }
          }
        }
      }
    }
    $search_term = $accession . '_' . $general_hash->{'import_type'};
    #Check if assembly already registered 
    if (!exists $existing_assemblies{$search_term}){
      say "Registration in progress for assembly $accession";
      #Assign assembly to the specific project (i.e., DToL, EBP, VGP, etc) it belongs to 
      my $assembly_group = &get_bio($bioproject);
      chomp($assembly_group);
      my ($chain,$version) = $self->split_gca($accession);
      #If species already exists, re-use species prefix
      if (exists $existing_species_prefix{$taxon_id}){
        $species_prefix = $existing_species_prefix{$taxon_id};
        say "Using existing species prefix $species_prefix";
        if (defined($rnaseq_status) && ($rnaseq_status =~ /amber|available|green|not available|red/)){
          #Registering assembly for an existing species using current trascriptomic data status
          register_assembly($chain, $version, $parent_id, $species_prefix, $clade, $taxon_id, $assembly_name, $assembly_level, $wgs_id, $scaffold_N50, $contig_N50, $contig_cnt, $scaffold_cnt, $assembly_type, $refseq_accession, $date, $top_level_cnt, $total_length, $total_gap_length, $assembly_ftp, $subspecies_name, $species_common_name, $genome_rep, $assembly_group, $species_name, $rnaseq_status, $submitter,$assembly_registry);
        }
      }
      else{
        #generate new stable id prefix for species being registered for the first time
        #split species name to find out if it is binomial or trinomial
        my @s_name = split(' ', $species_name);
        my ($first,$second,$third,$prefix,$reg_status,$species_prefix);
        #if array element is > 2, check that species name does not contain commmon name 
        if (scalar(@s_name) > 2){
          $s_name[0] =~ s/[\s\W\d_]//g; $s_name[1] =~ s/[\s\W\d_]//g; $s_name[2] =~ s/[\s\W\d_]//g;
          $first = uc(substr $s_name[0], 0, 1);
          $second = uc(substr $s_name[1], 0, 1);
          $third = uc(substr $s_name[2], 0, 1);
          $prefix = "ENS". $first.$second.$third;
          $reg_status = 3;
          $species_prefix = &iterate_prefix($reg_status, $prefix, \@s_name, \%existing_species_prefix, $parent_id, $accession, $species_name, 0);
        }
        else{#if name is binomial, extract first character from name1 and the first two characters from name2
          $s_name[0] =~ s/[\s\W\d_]//g; $s_name[1] =~ s/[\s\W\d_]//g;
          $first = uc(substr $s_name[0], 0, 1);
          $second = uc(substr $s_name[1], 0, 2);
          $prefix = "ENS". $first.$second;
          $reg_status = 2;
          $species_prefix = &iterate_prefix($reg_status, $prefix, \@s_name, \%existing_species_prefix, $parent_id, $accession, $species_name, 0);
        }
        #updating accession, prefix and species hash to include newly registered species
        my @new_record = split (/\t/, $species_prefix);
        if (@new_record){
          $species_prefix = $new_record[0];
          $existing_species_prefix{$species_prefix} = $new_record[1];
          $existing_assemblies{$new_record[2]} = 1;
          $existing_species_id{$new_record[1]} = $new_record[0];
        }
        #register species for the first time
        register_assembly($chain, $version, $parent_id, $species_prefix, $clade, $taxon_id, $assembly_name, $assembly_level, $wgs_id, $scaffold_N50, $contig_N50, $contig_cnt, $scaffold_cnt, $assembly_type, $refseq_accession, $date, $top_level_cnt, $total_length, $total_gap_length, $assembly_ftp, $subspecies_name, $species_common_name, $genome_rep, $assembly_group, $species_name, $rnaseq_status, $submitter,$assembly_registry);
        say "Registered species is ", $accession . "\t" . $clade . "\t" . $species_name . "\t" . $species_common_name . "\t" . $assembly_level . "\t" . $assembly_name . "\t" . $assembly_group;
        #Update dynamic records with new entries
        $existing_assemblies{$accession} = 1;
        $existing_species_prefix{$taxon_id} = $species_prefix;
        $existing_species_id{$taxon_id} = 1;
      }
    }
    else{
      say "Accession exists, no further registration needed";
    } 
} 

=head2 iterate_prefix

 Arg [1]    : None
 Description: Generate prefix for new species
 Returntype : None
 Exceptions : 

=cut

sub iterate_prefix{
  my $flag = shift;
  my ($pid, @sname, $a, $b, $c, %pre, $taxid, $accession);
  $pid = shift;
  @sname = @{+shift};
  %pre = %{+shift};
  $taxid = shift;
  $accession = shift;
  my $name = shift;
  my $rev_cnt = shift;
  my @temp_name;
  
  #Check if initial prefix generated already exists in database
  if (exists $pre{$pid}){
    #generate another species prefix
    if ($flag == 2){#generate prefix for binomial species name
      $a = uc(substr $sname[0], 0, 1);
      $b = uc(substr $sname[1], 1, 1);
      if (!$b){
        #handle cases where prefix can't be generated because specie has just one character
      }
      else{
      #say "removing conflicting character";
      $sname[1] =~ s/$b//gi;
      $b = uc(substr $sname[1], 0, 2);
      $pid = "ENS".$a.$b;
      }
      undef @temp_name;#initialise array for temporarily holding prefixes
      #store updated species name
      push @temp_name, $sname[0];
      push @temp_name, $sname[1];
    }
    elsif ($flag == 3){#generate prefix for trinomial species name
      #generating another prefix for trinomial
      $a = uc(substr $sname[0], 0, 1);
      $b = uc(substr $sname[1], 0, 1);
      $c = uc(substr $sname[2], 0, 1);
      if ((!$c) or (length($sname[2]) == 1)){
        #completely discard the third name. Now generate prefix from just the first two names
        $flag = 2;
        $a = uc(substr $sname[0], 0, 1);
        $b = uc(substr $sname[1], 0, 2);
        $pid = "ENS".$a.$b;
        undef @temp_name;
        push @temp_name, $sname[0];
        push @temp_name, $sname[1];
      }
      else{
        $sname[2] =~ s/$c//gi;
        $c = uc(substr $sname[2], 0, 1);
        $pid = "ENS".$a.$b.$c;
        push @temp_name, $sname[0];
        push @temp_name, $sname[1];
        push @temp_name, $sname[2];
      }
    }
    else{}
    &iterate_prefix($flag, $pid, \@sname, \%pre, $taxid, $accession, $name, $rev_cnt); 
  }  
  else{
    #check that unique prefix meets standard length of 6;
    if (length($pid) <= 5){
      #reverse species name to generate new unique prefix
      if ($rev_cnt < 1){
        @temp_name = split(" ", $name);
        $sname[1] = scalar reverse $temp_name[1];
      }
      elsif ($rev_cnt < 2){
        @temp_name = split(" ", $name);
        $sname[1] = $temp_name[0];
        $sname[0] = $temp_name[1];
      }
      else{
        #generate new prefix at random
        $name =~ s/\s//g; #removing white space from name
        my @chars = split ('', $name);
        my $random_string;
        foreach (1..3){
          # rand @chars will generate a random 
          # number between 0 and scalar @chars
          $random_string.=$chars[rand @chars];
        }
        $pid = "ENS".$random_string;
        if (exists $pre{$pid}){
          $pid = substr($pid,-1,1);
          &iterate_prefix($flag, $pid, \@sname, \%pre, $taxid, $accession, $name, $rev_cnt);
        }
        else{
          return $pid . "\t$taxid\t$accession";
        }
      }
      if ($flag == 2){
        $a = uc(substr $sname[0], 0, 1);
        $b = uc(substr $sname[1], 0, 2);
        $pid = "ENS".$a.$b;
      }
      else{
        $a = uc(substr $sname[0], 0, 1);
        $b = uc(substr $sname[1], 0, 1);
        $c = uc(substr $sname[2], 0, 1);
        $pid = "ENS".$a.$b.$c;
      }
      $rev_cnt++;
      &iterate_prefix($flag, $pid, \@sname, \%pre, $taxid, $accession, $name, $rev_cnt);
    }
    else{
      #final prefix
      $pid = $pid . "\t" . $taxid . "\t" . $accession;
      return $pid;
    }
  }        
}

=head2 register_assembly

 Arg [1]    : None
 Description: Function to perform registration of new assembly(ies)
 Returntype : None
 Exceptions : 

=cut

sub register_assembly{
  my ($self) = @_;
  my $genbank_acc = $_[0] . "." . $_[1];#join chain and version to get complete accession
  my %reg_asm; my $flag = 0; my $report; my $rnaseq_status = ""; my @rep;
  my $start = 0; my $end = 0; my $new_space_id = 0; my $max_version = 0;
  my $assembly_registry = $_[scalar(@_)-1]; 
  #check species space table to see if any record exist for the species
  my $sql = "SELECT current_space from species_space_log WHERE species_id = $_[5]";
  my $sth = $assembly_registry->dbc->prepare($sql);
  $sth->execute();
  my $current_space = $sth->fetchrow_array();
  my $anno_flag = $general_hash->{'import_type'};
  chomp($anno_flag);
  if ($anno_flag !~ /import_refseq|import_flybase|import_community|import_genbank|import_wormbase/){
    $anno_flag = "unannotated";
  }
  
  if ($current_space){
    #get existing versions for assembly chain
    $sql = "SELECT version from assembly WHERE chain = '$_[0]'";
    $sth = $assembly_registry->dbc->prepare($sql);
    $sth->execute();
    while (my $version = $sth->fetchrow_array()) {
      $max_version =  $version if ($max_version < $version);
    }
    if ($max_version){
      if (($max_version < $_[1]) || ($anno_flag =~ /import_refseq|import_flybase|import_community|import_genbank|import_wormbase/)){
        $flag = 2; #Assembly version update required so we re-use same stable id space
        $report = $_[0] . "." . $_[1] . "\t" . $_[4] . "\t" . $_[20] . "\t" . $_[21] . "\t" . $_[7] . "\t" . $_[6] . "\t" . $_[1] . "\t" . $_[25];
        &update_registry_db($flag, $_[0], $_[1], $current_space, $_[3], 0, $_[4], $anno_flag, $_[2], $_[5], $_[1], $_[6], $_[7], $_[8], $_[9], $_[10], $_[11], $_[12], $_[13], $_[14], $_[15], $_[16],
        $_[17], $_[18], $_[19], $_[20], $_[21], $_[22], $_[23], $_[25], $_[26], 'not started',$assembly_registry);
      }
      else{
        die("You should never have a lower version coming as an update");
        return;
      }
    }
    else{#means a different chain exists for species so we increment stable id space id value by 1
      $new_space_id = $current_space + 1;
      #check that new space id does not already exists
      $sql = "SELECT * from stable_id_space WHERE stable_id_space_id = $new_space_id";
      $sth = $assembly_registry->dbc->prepare($sql);
      $sth->execute();
      my $stable_id_space = $sth->fetchrow();
      if(!$stable_id_space) {
        $sql = "SELECT stable_id_space_start, stable_id_space_end from stable_id_space WHERE stable_id_space_id = $current_space";
        $sth = $assembly_registry->dbc->prepare($sql);
        $sth->execute();
        my ($stable_id_space_start,$stable_id_space_end) = $sth->fetchrow();
        $stable_id_space_start = $stable_id_space_start + 1;
        $stable_id_space_end = $stable_id_space_end + 5000000;
        #new stable id is being inserted into the stable id space table
        $flag = 3;
        $report =  $_[0] . "." . $_[1] . "\t" . $_[4] . "\t" . $_[20] . "\t" . $_[21] . "\t" . $_[7] . "\t" . $_[6] . "\t" . $_[26];
        &update_registry_db($flag, $stable_id_space_start, $stable_id_space_end, $_[0], $_[1], $new_space_id, $_[3], 0, $_[4], $anno_flag, $_[2], $_[5], $_[6], $_[7], $_[8], $_[9], $_[10], $_[11], 
        $_[12], $_[13], $_[14], $_[15], $_[16], $_[17], $_[18], $_[19], $_[20], $_[21], $_[22], $_[23],$_[25],$_[26], 'not started',$assembly_registry);
       
      }
      else{
        $flag = 4;
        $report =  $_[0] . "." . $_[1] . "\t" . $_[4] . "\t" . $_[20] . "\t" . $_[21] . "\t" . $_[7] . "\t" . $_[6] . "\t" . $_[26];
        &update_registry_db($flag, $_[0], $_[1], $new_space_id, $_[3], 0, $_[4], $anno_flag, $_[2], $_[5], $_[6], $_[7], $_[8], $_[9], $_[10], $_[11], $_[12], $_[13], $_[14], $_[15], $_[16], $_[17],
        $_[18], $_[19], $_[20], $_[21], $_[22], $_[23], $_[25], $_[26],'not started',$assembly_registry);
      }
    }
  }
  else{#Performing first time registration
    my $stable_id_space = 1;
    $flag = 1;
    &update_registry_db($flag, $_[0], $_[1], $stable_id_space, $_[3], 0, $_[4], $anno_flag, $_[2], $_[5], $_[6], $_[7], $_[8], $_[9], $_[10], $_[11], $_[12], $_[13], $_[14], $_[15], $_[16], $_[17], 
    $_[18], $_[19], $_[20], $_[21], $_[22], $_[23], $_[24],$_[25], $_[26], 'not started',$assembly_registry);
  }
  push @rep, ($flag,$report);

  return @rep;
}

=head2 update_registry_db

 Arg [1]    : None
 Description: Update registry db with new assembly details.
 Returntype : None
 Exceptions : 

=cut

sub update_registry_db{
  my ($self) = @_;
  my $assembly_registry = $_[scalar(@_)-1];
  #Register assembly for the first time
  if ($_[0] == 1){
    $assembly_registry->dbc->do("insert into species_space_log values (?,?,?)", undef, $_[9], $_[3], $_[24]);
    $assembly_registry->dbc->do("insert into assembly (chain, version, stable_id_space_id, species_prefix, is_current, clade, annotated_status, species_id, taxonomy, assembly_path, genome_rep, rnaseq_data) values (?,?,?,?,?,?,?,?,?,?,?,?)", undef, $_[1], $_[2], $_[3], $_[4], $_[5], $_[6], $_[7], $_[8], $_[9], $_[23], $_[26], $_[29]);
    #Store meta data using assembly id value just generated
    my $sql = "SELECT assembly_id from assembly WHERE chain = '$_[1]' and version = $_[2] and annotated_status = '$_[7]'";
    my $sth = $assembly_registry->dbc->prepare($sql);
    $sth->execute();
    my $assembly_id = $sth->fetchrow();
    if($assembly_id) {
      $assembly_registry->dbc->do("insert into meta (assembly_id, assembly_name, assembly_level, wgs_id, scaffold_N50, contig_N50, contig_count, scaffold_count, assembly_type, refseq_accession, assembly_date, top_level_count, total_length, total_gap_length, subspecies_name,common_name, assembly_group, submitter) values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", undef, $assembly_id, $_[10], $_[11], $_[12], $_[13], $_[14], $_[15], $_[16], $_[17], $_[18], $_[19], $_[20], $_[21], $_[22], $_[24], $_[25], $_[27], $_[30]);
    }
  }
  elsif ($_[0] == 2){#this handles the case when we are just updating the version of an assembly
    if ($_[29] =~ /null/i){
      $assembly_registry->dbc->do("insert into assembly (chain, version, stable_id_space_id, species_prefix, is_current, clade, annotated_status, species_id, taxonomy, assembly_path, genome_rep) values (?,?,?,?,?,?,?,?,?,?,?)", undef, $_[1], $_[2], $_[3], $_[4], $_[5], $_[6], $_[7], $_[8], $_[9], $_[24], $_[27]);
    }
    else{
      $assembly_registry->dbc->do("insert into assembly (chain, version, stable_id_space_id, species_prefix, is_current, clade, annotated_status, species_id, taxonomy, assembly_path, genome_rep, rnaseq_data) values (?,?,?,?,?,?,?,?,?,?,?,?)", undef, $_[1], $_[2], $_[3], $_[4], $_[5],$_[6], $_[7], $_[8], $_[9], $_[24], $_[27], $_[29]);
    }
    #Store meta data using assembly id value just generated
    my $sql = "SELECT assembly_id from assembly WHERE chain = '$_[1]' and version = $_[2] and annotated_status = '$_[7]'";
    my $sth = $assembly_registry->dbc->prepare($sql);
    $sth->execute();
    my $assembly_id = $sth->fetchrow();
    if($assembly_id) {
      $assembly_registry->dbc->do("insert into meta (assembly_id, assembly_name, assembly_level, wgs_id, scaffold_N50, contig_N50, contig_count, scaffold_count, assembly_type, refseq_accession, assembly_date, top_level_count, total_length, total_gap_length, subspecies_name, common_name, assembly_group, submitter) values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", undef, $assembly_id, $_[11], $_[12], $_[13], $_[14], $_[15], $_[16], $_[17], $_[18], $_[19], $_[20], $_[21], $_[22], $_[23], $_[25], $_[26], $_[28], $_[30]);
    }
  }
  elsif ($_[0] == 3){#this handles the case when a new assembly for an existing species is being registered
    $assembly_registry->dbc->do("insert into stable_id_space (stable_id_space_id, stable_id_space_start, stable_id_space_end) values (?,?,?)", undef, $_[5], $_[1], $_[2]);
    if ($_[30] =~ /null/i){
      $assembly_registry->dbc->do("insert into assembly (chain, version, stable_id_space_id, species_prefix, is_current, clade, annotated_status, species_id, taxonomy, assembly_path, genome_rep) values (?,?,?,?,?,?,?,?,?,?,?)", undef, $_[3], $_[4], $_[5], $_[6], $_[7], $_[8], $_[9], $_[10], $_[11], $_[25], $_[28]);
    }
    else{
      $assembly_registry->dbc->do("insert into assembly (chain, version, stable_id_space_id, species_prefix, is_current, clade, annotated_status, species_id, taxonomy, assembly_path, genome_rep, rnaseq_data) values (?,?,?,?,?,?,?,?,?,?,?,?)", undef, $_[3], $_[4], $_[5], $_[6], $_[7],$_[8], $_[9], $_[10], $_[11], $_[25], $_[28], $_[30]);
    }
    #Store meta data using assembly id value just generated
    my $sql = "SELECT assembly_id from assembly WHERE chain = '$_[1]' and version = $_[2] and annotated_status = '$_[7]'";
    my $sth = $assembly_registry->dbc->prepare($sql);
    $sth->execute();
    my $assembly_id = $sth->fetchrow();
    if($assembly_id) {
      $assembly_registry->dbc->do("insert into meta (assembly_id, assembly_name, assembly_level, wgs_id, scaffold_N50, contig_N50, contig_count, scaffold_count, assembly_type, refseq_accession, assembly      _date, top_level_count, total_length, total_gap_length, subspecies_name,common_name, assembly_group, submitter) values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", undef, $assembly_id, $_[12], $_[13],       $_[14], $_[15], $_[16], $_[17], $_[18], $_[19], $_[20], $_[21], $_[22], $_[23], $_[24], $_[26], $_[27], $_[29], $_[31]);
      $assembly_registry->dbc->do("update species_space_log set current_space = ? where species_id = ?", undef, $_[5], $_[11]);
    }
  }
  elsif ($_[0] == 4){#this handles the case when a new assembly for an existing species is being registered
    if ($_[28] =~ /null/i){
      $assembly_registry->dbc->do("insert into assembly (chain, version, stable_id_space_id, species_prefix, is_current, clade, annotated_status, species_id, taxonomy, assembly_path, genome_rep) values (?,?,?,?,?,?,?,?,?,?,?)", undef, $_[1], $_[2], $_[3], $_[4], $_[5], $_[6], $_[7],$_[8], $_[9], $_[23], $_[26]);
    }
    else{
      $assembly_registry->dbc->do("insert into assembly (chain, version, stable_id_space_id, species_prefix, is_current, clade, annotated_status, species_id, taxonomy, assembly_path, genome_rep, rnaseq_data) values (?,?,?,?,?,?,?,?,?,?,?,?)", undef, $_[1], $_[2], $_[3], $_[4], $_[5], $_[6], $_[7], $_[8], $_[9], $_[23], $_[26], $_[28]);
    }
    #Store meta data using assembly id value just generated
    my $sql = "SELECT assembly_id from assembly WHERE chain = '$_[1]' and version = $_[2] and annotated_status = '$_[7]'";
    my $sth = $assembly_registry->dbc->prepare($sql);
    $sth->execute();
    my $assembly_id = $sth->fetchrow();
    if($assembly_id) {
      $assembly_registry->dbc->do("insert into meta (assembly_id, assembly_name, assembly_level, wgs_id, scaffold_N50, contig_N50, contig_count, scaffold_count, assembly_type, refseq_accession, assembly_date, top_level_count, total_length, total_gap_length, subspecies_name,common_name, assembly_group, submitter) values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", undef, $assembly_id, $_[10], $_[11], $_[12], $_[13], $_[14], $_[15], $_[16], $_[17], $_[18], $_[19], $_[20], $_[21], $_[22], $_[24], $_[25], $_[27], $_[29]);
      $assembly_registry->dbc->do("update species_space_log set current_space = ? where species_id = ?", undef, $_[3], $_[9]);
    }
  }
  else{}
}
=head2 split_gca

 Arg [1]    : None
 Description: Get chain and version from assembly accession.
 Returntype : None
 Exceptions : 

=cut

sub split_gca {
  my ($self,$chain_version) = @_;

  unless($chain_version =~ /^(GCA\_\d{9})\.(\d+)$/) {
    die("Could not parse versioned GCA. GCA used: ".$chain_version);
  }

  my $chain = $1;
  my $version = $2;

  return($chain,$version);
}

=head2 parse_assembly_report

 Arg [1]    : None
 Description: Download and extract meta data from assembly report.
 Returntype : None
 Exceptions : Throws if no file exists

=cut

sub parse_assembly_report{
  my ($ftp,$accession,$assembly_name,$full_assembly_path) = @_;
  my ($available,$taxon_id,$assembly_level,$refseq_accession,$wgs_id,$bioproject,$genome_rep,$total_length,$scaffold_cnt,$contig_cnt,$contig_N50,$scaffold_N50,$total_gap_length,$top_level_cnt,$assembly_type,$date,$submitter,$common_name);
  $ftp->cwd();
  $ftp->cwd($full_assembly_path);
  my @assembly_report;
  my $report_file_name = $accession."_".$assembly_name."_assembly_stats.txt";
  my $report_file_content;
  my $report_file_handle;
  my $assembly_ftp = 'https://ftp.ncbi.nlm.nih.gov'.$full_assembly_path;
  open($report_file_handle, '>', \$report_file_content) || die("could not open $report_file_name");
  unless($ftp->get($report_file_name, $report_file_handle)) {
    die("Failed to retrieve the assembly report file: ", $ftp->message);
  }
  my @report_file_content = split("\n",$report_file_content);
  $available = 1;
  foreach my $line (@report_file_content) {
    if ($line =~ /^#\s+([^:]+):\s+(.+)/){
      if ($1 eq 'Assembly level') {
        $assembly_level = $2;
      }
      elsif ($1 eq 'Submitter') {
        $submitter = $2;
      }
      elsif ($1 eq 'Taxid') {
        $taxon_id = $2;
      }
      elsif ($1 eq 'WGS project') {
        $wgs_id = $2;
      }
      elsif ($1 eq 'RefSeq assembly accession') {
        $refseq_accession = $2;
      }
      elsif ($1 eq 'BioProject') {
        $bioproject = $2;
      }
      elsif ($1 eq 'Genome representation') {
        $genome_rep = $2;
      }
      elsif ($1 eq 'Date') {
        $date = $2;
      }
      elsif($1 eq 'Organism name') {
        $common_name = $2;
        $common_name = $1 if ($common_name =~ /\(([^)]+)\)/);
        $common_name =~ s/\s+\(.+//;
      }
      elsif($1 eq 'Assembly type') {
        $assembly_type = $2;
      }
     
    }
    else{
      unless($line =~ /^all/) {
       next;
      }
      if($line =~ /total-length\s*(\d+)/) {
       $total_length = $1;
      } elsif($line =~ /scaffold-count\s*(\d+)/) {
        $scaffold_cnt = $1;
      } elsif($line =~ /contig-count\s*(\d+)/) {
        $contig_cnt = $1;
      } elsif($line =~ /contig-N50\s*(\d+)/) {
        $contig_N50 = $1;
      } elsif($line =~ /scaffold-N50\s*(\d+)/) {
        $scaffold_N50 = $1;
      } elsif($line =~ /total-gap-length\s*(\d+)/) {
        $total_gap_length = $1;
      } elsif($line =~ /top-level-count\s*(\d+)/) {
        $top_level_cnt = $1;
      }
    }
  }
  if (!defined(($scaffold_cnt)) || ($scaffold_cnt eq '')){
    $scaffold_cnt = 0;
  }
  if (!defined($contig_cnt) || ($contig_cnt eq '')){
    $contig_cnt = 0;
  }
  if (!defined($scaffold_N50) || ($scaffold_N50 eq '')){
    $scaffold_N50 = 0;
  }
  if (!defined($contig_N50) || ($contig_N50 eq '')){
    $contig_N50 = 0;
  }
  if (!defined($top_level_cnt) || ($top_level_cnt eq '')){
    $top_level_cnt = 0;
  }
  if (!defined($total_length) || ($total_length eq '')){
    $total_length = 0;
  }
  push @assembly_report, ($taxon_id,$assembly_level,$refseq_accession,$wgs_id,$bioproject,$genome_rep,$total_length,$scaffold_cnt,$contig_cnt,$contig_N50,$scaffold_N50,$total_gap_length,$top_level_cnt,$assembly_type,$date,$assembly_ftp,$submitter,$common_name);
  return @assembly_report;
}

=head2 get_bio

 Arg [1]    : bioproject id
 Description: Get project associated with assembly via bioproject id.
 Returntype : None
 Exceptions : 

=cut

sub get_bio{
  my $group = "";
  my $bio_proj = $_[0];
  #Get all Bioproject ids linked to major projects such as DToL, EBP and VGP 
  my @ebp_main_bio_prj = `esearch -query 'PRJNA533106' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc  `;
  my @dtol_main_bio_prj = `esearch -query 'PRJEB40665' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc  `;
  my @vgp_main_bio_prj = `esearch -query 'PRJNA489243' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc  `;
  my @ergapp_main_bio_prj = `esearch -query 'PRJEB47820' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc  `;
  #Get list of projects linked to the same bioproject as assembly being registered
  my @pri_proj = `esearch -query '$bio_proj' -db bioproject | elink -related | efetch -format docsum | xtract -pattern DocumentSummary -element Project_Acc  `;
  if (@pri_proj){
    foreach $bio_proj (@pri_proj){
      chomp($bio_proj);
      if ($bio_proj eq 'PRJNA533106'){#Because DToL is linked to EBP, it is best not to mistakenly assign an EBP assembly to a DToL, hence we first check that the assembly is not under EBP
        return "ebp";
      }
      if (grep (/$bio_proj/,@dtol_main_bio_prj)) {
        return "dtol";
      }
      elsif (grep (/$bio_proj/,@vgp_main_bio_prj)) {
        return "vgp";
      }
      elsif (grep (/$bio_proj/,@ebp_main_bio_prj)) {
         return "ebp";
      }
      elsif (grep (/$bio_proj/,@ergapp_main_bio_prj)) {
         return "ergap";
      }
      else{$group = "ungrouped";}
    }
  }
  else{$group = "ungrouped";}
  return $group; 

}

=head2 get_clade

 Arg [1]    : taxon id
 Description: Get clade for species via taxon id.
 Returntype : None
 Exceptions : 

=cut

sub get_clade{
  my ($self,$taxon_id,$output_path) = @_;
  my ($lineage,$rank,$species_id,$clade,$common_name,$species_name);
  my @result;
  my %clade_setting = ('Rodentia' => 'rodentia', 'Primates' => 'primates', 'Mammalia' => 'mammalia', 'Amphibia' => 'amphibians', 'Teleostei' => 'teleostei', 'Marsupialia' => 'marsupials', 'Aves' => 'aves', 'Sauropsida', => 'reptiles', 'Chondrichthyes' => 'sharks', 'Eukaryota' => 'non_vertebrates', 'Metazoa' => 'metazoa', 'Viral' => 'viral', 'Viruses' => 'viral', 'Viridiplantae' => 'plants', 'Arthropoda' => 'arthropods', 'Lepidoptera' => 'lepidoptera', 'Insecta' => 'insects', 'Hymenoptera' => 'hymenoptera', 'Hemiptera' => 'hemiptera', 'Coleoptera' => 'coleoptera', 'Diptera' => 'diptera', 'Mollusca' => 'mollusca', 'Vertebrata' => 'vertebrates', 'Alveolata' => 'protists', 'Amoebozoa' => 'protists', 'Choanoflagellida' => 'protists', 'Cryptophyta' => 'protists', 'Euglenozoa' => 'protists', 'Fornicata' => 'protists', 'Heterolobosea' => 'protists', 'Parabasalia' => 'protists', 'Rhizaria' => 'protists', 'Stramenopiles' => 'protists', 'Fungi' => 'fungi');

  #Fetching taxonomy info from  NCBI
  my $taxonomy_report = $output_path."/taxonomy_report.txt";
  #call python script to query NCBI taxonomy db for taxonomy report
  system("python3 taxon.py --taxon_id $taxon_id --taxonomy_report $taxonomy_report");
  open(my $dbcore, '<:encoding(UTF-8)', $taxonomy_report) or die "Could not open file '$taxonomy_report' $!";
  my @dblist = <$dbcore>;
  if ($dblist[scalar((@dblist)-1)] =~ /species/){#the species being processed is a subspecies
    my @parent = split(/:/,$dblist[scalar((@dblist)-1)]);
    $species_id = $parent[2];
  }
  else{#We are dealing with a species
    $species_id = $taxon_id;
  }
  if ($dblist[0] =~ /common_name/){#get the common and scientific names of species being processed
    my @c_name = split(/:/,$dblist[0]);
    if (defined($c_name[1])){
      $common_name = $c_name[1];
    }
    $species_name = $c_name[2];
  }
  for my $taxa (reverse @dblist){
    my @line = split(/:/,$taxa);
    $line[1]=~ s/^\s+|\s+$//g;
    if (exists($clade_setting{$line[1]})) {
      if ($species_id == 9606){#specific clade for homo sapiens 
        $clade = "humans";
        last;
      }
      else{
        $clade = $clade_setting{$line[1]};
        last;
      }
    }
  }

  if (!defined($clade)){
    $clade = "non_vertebrates";
  }
  push @result, ($species_id,$clade,$species_name,$common_name,$species_name);
  `rm $taxonomy_report`;
  return @result; 
    
}

=head2 write_output

 Arg [1]    : None
 Description: None
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my ($self) = @_;
}

}
1;

