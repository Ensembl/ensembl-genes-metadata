=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Pipeline::Runnable::FetchTranscriptomicData

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Pipeline::Runnable::FetchTranscriptomicData;

use strict;
use warnings;

use JSON::PP;
use LWP::UserAgent;
use File::Spec::Functions qw(splitpath);
use feature 'say';
use Data::Dumper;
use File::stat;
use File::Spec::Functions;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 param_defaults

 Arg [1]    : None
 Description: Default parameters for the analysis
               ena_base_url => 'https://www.ebi.ac.uk/ena/portal/api/search?display=report',
               files_domain => 'domain=read&result=read_run',
               files_fields => 'run_accession,study_accession,experiment_accession,sample_accession,secondary_sample_accession,instrument_platform,instrument_model,library_layout,library_strategy,nominal_length,read_count,base_count,fastq_ftp,fastq_aspera,fastq_md5,library_source,library_selection,center_name,study_alias,experiment_alias,experiment_title,study_title',
               sample_domain => 'domain=sample&result=sample',
               sample_fields => 'accession,secondary_sample_accession,bio_material,cell_line,cell_type,collected_by,collection_date,country,cultivar,culture_collection,description,dev_stage,ecotype,environmental_sample,first_public,germline,identified_by,isolate,isolation_source,location,mating_type,serotype,serovar,sex,specimen_voucher,strain,sub_species,sub_strain,tissue_lib,tissue_type,variety,tax_id,scientific_name,sample_alias,center_name,protocol_label,project_name,investigation_type,experimental_factor,sample_collection,sequencing_method',
               download_method => 'ftp',
               separator => '\t',
               _read_length => 1, # This is a default that should not exist. Some data do not have the read count and base count
               _centre_name => 'ENA',
               print_all_info => 0, # It will print all the sample information in the description instead of a selection, good to use when checking the CSV file
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
    sample_fields => 'accession,secondary_sample_accession,bio_material,cell_line,cell_type,collected_by,collection_date,country,cultivar,culture_collection,description,dev_stage,ecotype,environmental_sample,first_public,germline,identified_by,isolate,isolation_source,location,mating_type,serotype,serovar,sex,specimen_voucher,strain,sub_species,sub_strain,tissue_lib,tissue_type,variety,tax_id,scientific_name,sample_alias,center_name,protocol_label,project_name,investigation_type,experimental_factor,sample_collection,sequencing_method',
    download_method => 'ftp',
    separator => '\t',
    _read_length => 1,
    _centre_name => 'ENA',
    print_all_info => 0,
    paired_end_only => 1, #by default, module will only add paired-end data to the csv, add "paired_end_only => 0" to pipeline config to include single end data 
    read_type => 'short_read',
    instrument_platform => 'ILLUMINA',  
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Create the query to be based on 'study_accession' or 'taxon_id'.
              An array of values can be given.
              If the 'output_dir' already exists, it will complete early
 Returntype : None
 Exceptions : None

=cut

sub fetch_input {
  my ($self) = @_;
  my $csv_file = "";
  my $csv_full = "";
  $self->param_required('output_dir'); 
  $self->param_required('sp_name');
  my $dir = $self->param('sp_name');
  $dir =~ s/^\s+|\s+$//g;
  $dir =~ s/ /_/g;
  $dir =~ s/'//g;
  $dir =~ s/\./_/g;
  $self->param('species_path', $dir);
  if($self->param_required('read_type') eq 'isoseq') {
    $self->param('paired_end_only',0);
    $self->param('instrument_platform','PACBIO_SMRT'),
    $csv_file = $self->param('output_dir') . $dir . "/" . $dir . "_isoseq.csv";
    #make an alternative csv file that would hold similar formatted data as the short read file
    #this is because the format for the csv file used by minimap is different from that used by bwa
    $csv_full = $self->param('output_dir') . $dir . "/" . $dir . "_isoseq_full.csv";
  }
  else{
  	  $csv_file = $self->param('output_dir') . $dir . "/" . $dir . "_rnaseq.csv";
  }
  if (-e $csv_file) {
    `rm $csv_file $csv_full`;
    $self->_populate_query($self->param('sp_id'), 'tax_tree(%s) AND instrument_platform='.$self->param('instrument_platform').' AND library_source=TRANSCRIPTOMIC');
    #$self->complete_early("'output_dir' exists so I will use that");
  } elsif ($self->param_is_defined('study_accession') and $self->param('study_accession')) {
    $self->_populate_query($self->param('study_accession'), 'study_accession=%s');
  } elsif ($self->param_is_defined('sp_id')){# and $self->param('genus')) {
    $self->_populate_query($self->param('sp_id'), 'tax_tree(%s) AND instrument_platform='.$self->param('instrument_platform').' AND library_source=TRANSCRIPTOMIC');
  } else {
  }
  $self->param('csvfile', $csv_file);
  $self->param('csvfull', $csv_full);
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


=head2 run

 Arg [1]    : None
 Description: Query ENA to find all possibles project that could be used for the RNASeq pipeline
              It will avoid samples which are for non coding RNA analyses or if the individual was
              infected, immunised,...
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;
  my $outdir = $self->param('output_dir') . $self->param('sp_name') . '/';
  my $tax = $outdir . $self->param('sp_id') . ".txt";
  my $csv_path = catdir($outdir, 'csv');
  my $tissue_dict = $self->param('tissue_dict');
  my %csv_data;
  my %samples;
  my $tax_level = $self->param('output_dir') . "taxlevel.txt";
  my $data_level = $self->param('output_dir') . "data_level.txt";
  open DL, (">>$data_level");
  foreach my $query (@{$self->param('query')}) {
    my $ua = LWP::UserAgent->new;
    $ua->env_proxy;
    my $url = join('&', $self->param('ena_base_url'), 'query="'.$query.'"', $self->param('files_domain'), 'fields='.$self->param('files_fields'));
    say "url passed to ftp is $url";
    $self->say_with_header($url);
    my $response = $ua->get($url);
    my $fastq_file = 'fastq_'.$self->param('download_method');
    if ($response->is_success) {
      my $content = $response->decoded_content(ref => 1);
      open M, (">$tax");
      print M Dumper($content);
      if ($content) {
        my %fields_index;
        while ($$content =~ /^(\w+.*)$/mgc) {
          my $line = $1;
          if ($line =~ /^[a-z]/) {
            my $index = 0;
            %fields_index = map { $_ => $index++} split('\t', $line);
          }
          else {
            say "line does not match a-z and line is \t$line";
            # if these two checks below are removed, more time might be needed to prepare the CSV file
            next if ($line =~ / infected | [iIu]mmune| challenge |tomi[zs]ed/); # I do not want to do that but I don't think we have a choice
            if ($line =~ /([Mm]i\w{0,3}RNA)|lncRNA|circRNA|small RNA/){
            }
            if ($line =~ /([Mm]i\w{0,2}RNA)|lncRNA|circRNA|small RNA/){
            }
            if ($line =~ /[Mm]is\w{0,2}RNA|lncRNA|circRNA|small RNA/){ # I do not want to do that but I don't think we have a choice
            }
            else{
            }
     #       say "line does not match a-z and line is \t$line";
            my @row = split("\t", $line);
            my $read_length = $self->param('_read_length');
            my $nominal_length = 0;
            my $calculated_length = 0;

            if ($self->param('paired_end_only')){
               my $third_file = "ftp[^_]*\.fastq\.gz\;ftp";
               $row[$fields_index{$fastq_file}] =~ s/$third_file/ftp/; # sometimes a third combined fastq file exists (it has single end naming format) - discard it

              next if ($row[$fields_index{library_layout}] eq 'SINGLE'); # don't include single end reads
              next if ($row[$fields_index{$fastq_file}] !~ m/.*_1\.fastq\.gz.*_2\.fastq\.gz/);# this will throw out the paired end data that is stored in a single file, i.e. it looks like single end data to a regex
            }

            if ($row[$fields_index{nominal_length}]) {
              $nominal_length = $row[$fields_index{nominal_length}];
              $read_length = $nominal_length;
            }
            if ($nominal_length != $calculated_length and $nominal_length != 0 and $calculated_length != 0) {
              $self->warning("${row[$fields_index{run_accession}]} NOMINAL $nominal_length CALC $calculated_length");
            }
            if ($row[$fields_index{library_layout}] eq 'PAIRED') {
              my @files = split(/;/,$row[$fields_index{$fastq_file}]);
              next if ($files[1] !~ /_2\.fastq\.gz/)
                
              
            }
            $read_length = POSIX::ceil($read_length);
            say "read length = $read_length";
            my %line = (
              run_accession => $row[$fields_index{run_accession}],
              instrument_model => $row[$fields_index{instrument_model}],
              instrument_platform => $row[$fields_index{instrument_platform}],
              library_layout => $row[$fields_index{library_layout}],
              fastq_file => $row[$fields_index{$fastq_file}],
              fastq_md5 => $row[$fields_index{fastq_md5}],
              study_title => $row[$fields_index{study_title}],
              experiment_title => $row[$fields_index{experiment_title}],
              instrument_model => $row[$fields_index{instrument_model}],
              read_length => $read_length,
              center_name => $row[$fields_index{center_name}],
              fastq_ftp => $row[$fields_index{fastq_ftp}],
            );
            $samples{$row[$fields_index{sample_accession}]} = $row[$fields_index{secondary_sample_accession}];
            push(@{$csv_data{$row[$fields_index{study_accession}]}->{$row[$fields_index{sample_accession}]}}, \%line);
          }
        }
        
        my @sample_names = keys %samples;
        my $header;
        unless (-d "$csv_path"){
          `mkdir $csv_path -p`;
        }
        open V, (">>$csv_path/samples_off.txt");
        SAMPLE: foreach my $sample (@sample_names) {
          
          $url = join('&', $self->param('ena_base_url'), 'query="accession='.$sample.'"', $self->param('sample_domain'), 'fields='.$self->param('sample_fields'));
          say "Sample URL is $url";
          $response = $ua->get($url);
          if ($response->is_success) {
            $content = $response->decoded_content();
            open N, (">$csv_path/$sample.csv"); print N Dumper($content);
            if ($content) {
              while ($content =~ /^(\w+.*)$/mgc) {
                my $line = $1;
                if ($line =~ /^[a-z]/) {
                  my $index = 0;
                  %fields_index = map { $_ => $index++} split('\t', $line);
                  $header = $line;
                }
                else {
                  $self->say_with_header($line);
                  my @row = split("\t", $line);
                  my %line = (
                    center_name => $row[$fields_index{center_name}],
                    cell_line => $row[$fields_index{cell_line}],
                    cell_type => $row[$fields_index{cell_type}],
                    dev_stage => $row[$fields_index{dev_stage}],
                    sex => $row[$fields_index{sex}],
                    strain => $row[$fields_index{strain}],
                    sub_species => $row[$fields_index{sub_species}],
                    sub_strain => $row[$fields_index{sub_strain}],
                    tissue_lib => $row[$fields_index{tissue_lib}],
                    tissue_type => $row[$fields_index{tissue_type}],
                    variety => $row[$fields_index{variety}],
                    tax_id => $row[$fields_index{tax_id}],
                    description => $row[$fields_index{description}],
                    sample_collection => $row[$fields_index{sample_collection}],
                    sequencing_method => $row[$fields_index{sequencing_method}],
                    sample_alias => $row[$fields_index{sample_alias}],
                  );
                  my $dh = $ua->default_headers;
                  $ua->default_header('Content-Type' => 'application/json');
                  my $biosd = $ua->get('http://www.ebi.ac.uk/biosamples/samples/'.$sample);
                  if ($biosd->is_success) {
                    $content = $biosd->decoded_content();
                    my $json = JSON::PP->new();
                    my $data = $json->decode($content);
                    open S, (">$csv_path/$sample.txt"); print S Dumper($data);
                    if (exists $data->{characteristics}->{immunization}) {
                    	print V $sample, "\n";
                      delete $samples{$sample};
                      $self->warning("Removed $sample from the set as it has immunization value: ".$data->{characteristics}->{immunization}->[0]->{text});
                      next SAMPLE;
                    }
                    if (exists $data->{characteristics}->{'developmental stage'}){
                      $line{dev_stage} = $data->{characteristics}->{'developmental stage'}->[0]->{text};
                    }
                    elsif (exists $data->{characteristics}->{'development stage'}){
                      $line{dev_stage} = $data->{characteristics}->{'development stage'}->[0]->{text};
                    }
                    elsif (exists $data->{characteristics}->{'stage'}){
                      $line{dev_stage} = $data->{characteristics}->{'stage'}->[0]->{text};
                    }
                    elsif (exists $data->{characteristics}->{'time point'}){
                      $line{dev_stage} = $data->{characteristics}->{'time point'}->[0]->{text};
                    }
                    else{$line{dev_stage} = 'non_existent';}
                    $line{status} = $data->{characteristics}->{healthStatusAtCollection}->[0]->{text}
                      if (exists $data->{characteristics}->{healthStatusAtCollection});
                    if (exists $data->{characteristics}->{age}) {
                      $line{age} = $data->{characteristics}->{age}->[0]->{text};
                      if (exists $data->{characteristics}->{age}->[0]->{unit}) {
                        $line{age} .= ' '.$data->{characteristics}->{age}->[0]->{unit};
                      }
                    }
                    if (exists $data->{characteristics}->{tissue} || exists $data->{characteristics}->{'tissue type'} || exists $data->{characteristics}->{tissue_type}) {
                      $line{organismPart} = $data->{characteristics}->{tissue}->[0]->{text};
                      if (exists $data->{characteristics}->{tissue}->[0]->{ontologyTerms}) {
                        $line{uberon} = $data->{characteristics}->{tissue}->[0]->{ontologyTerms}->[-1];
                      }
                    }
                    elsif (exists $data->{characteristics}->{'organism part'}) {
                      $line{organismPart} = $data->{characteristics}->{'organism part'}->[0]->{text} || $data->{characteristics}->{'tissue type'}->[0]->{text} ||
                                            $data->{characteristics}->{tissue_type}->[0]->{text};
                      if (exists $data->{characteristics}->{'organism part'}->[0]->{ontologyTerms}) {
                        $line{uberon} = $data->{characteristics}->{'organism part'}->[0]->{ontologyTerms}->[-1];
                      }
                      elsif (exists $data->{characteristics}->{tissue}->[0]->{ontologyTerms}) {
                        $line{uberon} = $data->{characteristics}->{tissue}->[0]->{ontologyTerms}->[-1];
                      }
                      elsif (exists $data->{characteristics}->{tissue_type}->[0]->{ontologyTerms}) {
                        $line{uberon} = $data->{characteristics}->{tissue_type}->[0]->{ontologyTerms}->[-1];
                      }
                      elsif (exists $data->{characteristics}->{'tissue type'}->[0]->{ontologyTerms}) {
                        $line{uberon} = $data->{characteristics}->{'tissue type'}->[0]->{ontologyTerms}->[-1];
                      }
                    }
                    elsif (exists $data->{characteristics}->{'source name'}) {
                      
                      $line{organismPart} = $data->{characteristics}->{'source name'}->[0]->{text};
                      
                    }
                    elsif (exists $data->{characteristics}->{'source_name'}) {
                      
                      $line{organismPart} = $data->{characteristics}->{'source_name'}->[0]->{text};
                      
                    }
                    elsif (exists $data->{characteristics}->{'body site'}) {
                      $line{organismPart} = $data->{characteristics}->{'body site'}->[0]->{text};
                      
                    }
                    elsif (exists $data->{characteristics}->{cellType}) {
                      $line{cellType} = $data->{characteristics}->{cellType}->[0]->{text};
                      if (exists $data->{characteristics}->{cellType}->[0]->{ontologyTerms}) {
                        $line{uberon} = $data->{characteristics}->{cellType}->[0]->{ontologyTerms}->[-1];
                      }
                    }
                    elsif (exists $data->{characteristics}->{'cell type'}) {
                      $line{cellType} = $data->{characteristics}->{'cell type'}->[0]->{text};
                      if (exists $data->{characteristics}->{'cell type'}->[0]->{ontologyTerms}) {
                        $line{uberon} = $data->{characteristics}->{'cell type'}->[0]->{ontologyTerms}->[-1];
                      }
                    }
                  }
                  else {
                    $self->warning("Could not connect to BioSample with $sample");
                    my $bioncbi = "$csv_path/$sample.txt";
                    my $sample_name;
                    my @sample;
                    #If biosample info can't be retrieved from ENA, try NCBI
                    `efetch -db biosample -id $sample -format text > $bioncbi`;
                    if (-e $bioncbi) {
                      if ($sample_name = `grep 'tissue' $bioncbi`){
                         @sample = split(/\"/,$sample_name);
                         say "sample is $sample[1]";
                         $line{organismPart} = $sample[1];
                       }
                       elsif ($sample_name = `grep 'source name' $bioncbi`){
                         @sample = split(/\"/,$sample_name);
                         say "sample is $sample[1]";
                         $line{organismPart} = $sample[1];
                       }
                       elsif ($sample_name = `grep 'sample name' $bioncbi`){
                         @sample = split(/\"/,$sample_name);
                         say "sample is $sample[1]";
                         $line{organismPart} = $sample[1];
                       }
                       elsif ($sample_name = `grep 'cell type' $bioncbi`){
                         @sample = split(/\"/,$sample_name);
                         say "sample is $sample[1]";
                         $line{organismPart} = $sample[1];
                       }
                       
                       else{$line{noBio} = 1;}
                    }
                    else{$line{noBio} = 1;}
                  }

                  $ua->default_headers($dh);
                  $samples{$sample} = \%line;
                  print V $row[$fields_index{tax_id}];
                }
              }
            }
          }
        }
      }
      else {
        $self->throw("There was a problem with '$url'");
      }
    }
    else {
      $self->warning('No results with "'.$url);
    }
  }
  if (keys %csv_data) {
    say "found keys";
    
    my $flag = 0;
    foreach my $project (keys %csv_data) {
    	
    	my $v = "";
    	say "Project tested is ", $project, " and it contains the following samples";
      my %dev_stages;
      my %celltypes;
      foreach my $sample (keys %{$csv_data{$project}}) {
        my $flag = 0;
         unless($sample =~ /^SAM/) {
          next;
         }
      	  say "sample in project is ", $sample;
        next unless (exists $samples{$sample});
       if (exists $samples{$sample}->{dev_stage} and ($samples{$sample}->{dev_stage} ne 'non_existent')) {
          
        	say "sample id is at development ", $sample;
          next if ($samples{$sample}->{dev_stage} eq 'sexually immature stage');
          $dev_stages{$samples{$sample}->{dev_stage}} = 1;
       }
      }
      if (scalar(keys(%dev_stages)) > 1) {
        
        foreach my $sample (keys %{$csv_data{$project}}) {
          
          unless($sample =~ /^SAM/) {
            next;
          }
          next unless (exists $samples{$sample});
          #check that the sample query has values set for dev_stage
          if (exists $samples{$sample}->{dev_stage} and $samples{$sample}->{dev_stage} ne '') {
            say "Dev stages passed";
            #check also that if biosample data fails, tissue type is set in sample query
            if ($samples{$sample}->{noBio} and $samples{$sample}->{tissue_type} ne ''){
              $samples{$sample}->{sample_name} =  $samples{$sample}->{tissue_type};
              say "No Biosamples found";
            }
            #this is when dev stage exists and we could download biosample
            elsif (!$samples{$sample}->{noBio} and $samples{$sample}->{organismPart} ne ''){
              $samples{$sample}->{sample_name} =  $samples{$sample}->{organismPart};
              say "Biosamples found with tissue";
            }
            else{
              
              #no tissue type found when biosamples could not be downloaded, use the sample id 
              $samples{$sample}->{sample_name} = $samples{$sample}->{tissue_type} || $samples{$sample}->{cellType} || $sample;
            }
          }
          else {
            #instead of breaking pipeline, show warning where sample within the project has no development stage even though other samples in the project have dev stage
            #use sample id for such cases
            $self->warning('No dev stages for '.$sample.' "'.join('", "', keys %dev_stages).'"');
            
            $samples{$sample}->{sample_name} = $samples{$sample}->{tissue_type} || $samples{$sample}->{cellType} || $samples{$sample}->{organismPart} || $sample;
            
          }

          if ($samples{$sample}->{sex}) {#we don't want sex values in sample name anymore
          }
          open D, ("<$tissue_dict");
          my $tissue_assignment = $self->param('csvfile');
          $tissue_assignment =~ s/.csv/.txt/;
          open BB, (">>$tissue_assignment");
          my @tissues;
          while (<D>){
            push @tissues, $_;
          }
          foreach my $tissue (@tissues){
            chomp($tissue);
            say "checking for tissue match $tissue\t",$samples{$sample}->{sample_name};
            if (index(lc($samples{$sample}->{sample_name}), lc($tissue)) != -1) {
              say "match found";
              $flag = 1;
              last;
            }
            
          }
          
          
          if ($flag == 0){
            say "setting sample name to sample id";
            
            $samples{$sample}->{sample_name} = $sample; 
           
          }
         
          close(D);
          close(BB);
          }
        
      }

      else {#project has no sample with development stage information, then take sample name as any of the below listed options in order
        	
        	foreach my $sample (keys %{$csv_data{$project}}) {
        	  unless($sample =~ /^SAM/) {
        	    next;
        	  }
        	  next unless (exists $samples{$sample});
        	  say "sample without dev is $sample";
          
          $samples{$sample}->{sample_name} = $samples{$sample}->{tissue_type} || $samples{$sample}->{cellType} || $samples{$sample}->{organismPart};# || $samples{$sample}->{sample_alias} || $samples{$sample}->{description};
          say "organism part is ", $samples{$sample}->{organismPart}, " and sample id is $sample" ;
          if (!$samples{$sample}->{sample_name}) {
            #Again instead of failing where no clear sample name has been given, use sample id
            $samples{$sample}->{sample_name} = $sample;
           
          }
          if ($samples{$sample}->{sex}) {#we do not want to use sex as part of sample name anymore
          }
          open D, ("<$tissue_dict");
          my $tissue_assignment = $self->param('csvfile');
          $tissue_assignment =~ s/.csv/.txt/;
          open BB, (">>$tissue_assignment");
          my @tissues;
          while (<D>){
            push @tissues, $_;
          }
          
          foreach my $tissue (@tissues){
            chomp($tissue);
            if (index(lc($samples{$sample}->{sample_name}), lc($tissue)) != -1) {
             
              $flag = 1;
              last;
            }
            else{
              print BB $samples{$sample}->{sample_name},"\n";
            }
            
          }
                    
          if ($flag == 0){
            say "setting sample name to sample id";
            print BB $samples{$sample}->{sample_name},"\t";
            $samples{$sample}->{sample_name} = $sample; 
            print BB $samples{$sample}->{sample_name},"\n";
          }
          
          close(D);
          close(BB);
          }
        
      }
      
    }
    $self->output([\%csv_data, \%samples]);
    $self->param('flow', 1);
    unless (-d $self->param('output_dir') . $self->param('sp_name')){
  	  my $query = "mkdir -p " . $self->param('output_dir') . $self->param('sp_name');
  	  `$query`;
  	}
  	#find species entry in initial candidate csv file and make a copy of species with rnaseq data
  	my $find = "grep -w " . $self->param('sp_name') . " " .  $self->param('species_list') . "_orig_bak >> " . $self->param('species_list');
  	system "$find";
  	if ($self->param('genus') == 1){
  	  #store data level as genus
  	  print DL  $self->param('sp_id') . "\t" . $self->param('sp_name') . "\tgenus\n";
  	  
  	}
  	else{
  	   #store data level as species
  	  print DL  $self->param('sp_id') . "\t" . $self->param('sp_name') . "\tspecies\n";
  	}
  }
  else {
  	$self->param('flow', 0);
    $self->complete_early('Could not find any data for this job');
    
  }
}
 
=head2 write_output

 Arg [1]    : None
 Description: It writes a csv file named 'output_dir' which can be process by the RNA-seq
              pipeline. The header should be
              SM\tID\tis_paired\tfilename\tis_mate_1\tread_length\tis_stranded\tCN\tPL\tDS
              is_mate_1 will always be -1 as we get the data from ENA so the filename will
                be informative of the mates
              is_stranded is 0 as we don't have a correct way of getting this information yet
              DS is made of "study_accession, sample_accession, study_title, experiment_title
                and cell_type if the sample is a cell_type
              It will also return the list of file to download on channel '_branch_to_flow_to',
              usually #2
 Returntype : None
 Exceptions : Throws if it cannot open or close the file 'output_dir'

=cut

sub write_output {
  my ($self) = @_;
  say "csv path is ", $self->param('csvfile');
  	  open(FH, '>>'.$self->param('csvfile')) || $self->throw('Could not open '.$self->param('csvfile'));
  my $data = $self->output;
  my $samples = $data->[1];
  my $download_method = $self->param('download_method');
  my @output_ids; my $source = "";
  foreach my $study_accession (keys %{$data->[0]}) {
    my $study = $data->[0]->{$study_accession};
    foreach my $sample (keys %{$study}) {
    	unless($sample =~ /^SAM/) {
        next;
      }
      next unless (exists $samples->{$sample});
      foreach my $experiment (@{$study->{$sample}}) {
        my @files = split(';', $experiment->{fastq_file});
	    my @fastq_url = split(';', $experiment->{fastq_ftp});
        my @checksums = split(';', $experiment->{fastq_md5});
        my $index = 0;
        foreach my $file (@files) {
          my (undef, undef, $filename) = splitpath($file);
          my $sample_name = $samples->{$sample}->{sample_name};
          $sample_name =~ s/\s+-\s+\w+:\w+$//;
          $sample_name =~ s/[[:space:][:punct:]]+/_/g;
          $sample_name =~ s/_{2,}/_/g;
          $sample_name =~ s/^_|_$//g;
          if ($sample_name !~ /^SAM/i){
            $sample_name =~ s/\d//g;
            $sample_name =~ s/^_//g;
          }
          if (length($sample_name) > 50){
            $sample_name = substr($sample_name,0,50);
          }
          my $description = sprintf("%s, %s%s",
            $experiment->{study_title},
            $experiment->{experiment_title},
            $samples->{$sample}->{cell_type} ? ', '.$samples->{$sample}->{cell_type} : '', );
          if ($self->param('print_all_info')) {
            foreach my $field (values %{$samples->{$sample}}) {
              $description .= ';'.$field if ($field);
            }
          }
          $description =~ tr/:\t/ /;
          if ($samples->{$sample}->{tax_id} eq ''){
             $samples->{$sample}->{tax_id} = $self->param('taxonomy');
          }
          my ($pair) = $filename =~ /\S+_(\d)\.\S+/;
          if (defined($pair) and ($experiment->{library_layout}) eq 'PAIRED'){
            if ($pair == 1) {
                $self->param('is_mate_1',1);
            }
            else {
                $self->param('is_mate_1',0);
            }
          }
          my $taxonomy = "";
          #Set the species id to the taxonomy value when data retrieved comes from genus level
          if ($self->param('genus') == 1){
            $taxonomy = $self->param('taxonomy');
          }
          else{
            $taxonomy = $self->param('sp_id');
          }
          if (($samples->{$sample}->{tax_id} eq '') || !defined($samples->{$sample}->{tax_id})){
            $source = $self->param('sp_id');
          }
          else{
            $source = $samples->{$sample}->{tax_id};
          }
          if($self->param('read_type') eq 'short_read') {
                           print FH sprintf("%s\t%s\t%s\t%s\t%s\t%s, %s, \t%s\t%s\t%s\t%s\n",
				lc($sample_name),
				$experiment->{run_accession},
				$experiment->{library_layout} eq 'PAIRED' ? 1 : 0,
			        $taxonomy,
				$filename,
				$study_accession,
				$sample,
				$file,
				$checksums[$index],
                                $self->param('sp_name'),
                                $source,
				
			  );
		  }
	elsif($self->param('read_type') eq 'isoseq') {
		
		print FH sprintf("%s\t%s\t%s\t%s, %s\t%s\t%s\t%s\t%s\n",
		lc($sample_name),
		$filename,
		$taxonomy,
                $study_accession,
                $sample,
		$file,
		$checksums[$index],
                $self->param('sp_name'),
                $source,
	      );
              #write to a second csv file that would be used during data classification
              #this file would be in similar format to the short read csv file
              open(LR, '>>'.$self->param('csvfull')) || $self->throw('Could not open '.$self->param('csvfull'));
	          print LR sprintf("%s\t%s\t%s\t%s\t%s, %s, %s\t%s\t%s\t%s\n",
				lc($sample_name),
				$experiment->{run_accession},
                                $taxonomy,
				$filename,
				$study_accession,
				$sample,
				$file,
				$checksums[$index],
                                $self->param('sp_name'),
                                $source,
			  );
          
          } else {
            $self->throw('Read type unknown: '.$self->param('read_type'));
          }
          push(@output_ids, {url => $file, download_method => $download_method, checksum => $checksums[$index++]});
        }
      }
    }
  }
  close(FH) || $self->throw('Could not close '.$self->param('csvfile'));
  if ($self->param('flow')){
  	$self->dataflow_output_id(\@output_ids, $self->param('_branch_to_flow_to'));
  }
}

1;

