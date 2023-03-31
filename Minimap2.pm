package Minimap2;

use strict;
use warnings;
use POSIX;
use List::Util qw( min max );
use feature 'say';
use Bio::EnsEMBL::DBSQL::DBConnection;
use File::Basename;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB'); 

sub fetch_input{
	 my ($self) = @_;
	unless (-d $self->param('output_dir')){
		my $dir = $self->param('output_dir');
		`mkdir $dir -p`;
	}
	
}

sub run{
	 my ($self) = @_;
	 my $read = $self->param('input_dir') . $self->param('iid');
	
	 my $genome_index = "";
	 if (-e($self->param('genome_file') . '.fasta.mmi'))
	 {
	   $genome_index = $self->param('genome_file') . '.fasta.mmi';
	  
	 }
	 else{
	 	 say "Indexing genome for minimap operation";
	 	 my $output = $self->param('genome_file') . '.fasta';
	 	 if (-e $output){
	 	 	 $genome_index = $self->param('genome_file') . '.fasta.mmi';
	 	 }
	 	 else{
	 	 	 $output = $self->param('genome_file') . '_toplevel.fasta';
	 	 	 $genome_index = $self->param('genome_file') . '_toplevel.fasta.mmi';
	 	 }
	 	my $minimap_index = $self->param('minimap2_path'). " -d $genome_index $output";
	 	if(system($minimap_index)) {
	  	  $self->throw("Error indexing genome via minimap2\nError code: $?\n");
	  }
	  
	 }
	 
	 my $file1 = "basename " . $self->param('iid') ." .fq";
	 $file1 = `$file1`;
	 my @run = split(/_/,$file1);
	 #use the run id as the output file name
	 my $sam = $run[0];
	 chomp($file1); chomp($sam);
	 ##Now need to separate between single alignment and paired alignment
	 if ($read =~ /\.fq/){
	   
	  # run minimap2
	  my $minimap2_command = $self->param('minimap2_path')." --cs -N 1 -ax splice:hq -u b ".$genome_index." ".$read." > ". $self->param('output_dir') . "/" . $sam. ".sam";
	  $self->warning("Command:\n".$minimap2_command."\n");
	  if(system($minimap2_command)) {
	  	  $self->throw("Error running minimap2\nError code: $?\n");
	  }
	   my $query = $self->param('samtools') . ' view -bS ' . $self->param('output_dir') . "/" . $sam . ".sam > " . $self->param('output_dir') . "/" . $sam . ".bam";
	   say "query is $query";
	   `$query`;
	   $query = $self->param('samtools') . ' sort -O bam ' . $self->param('output_dir') . "/" . $sam . ".bam -o " . $self->param('output_dir') . "/" . $sam . ".sorted.bam";
	   say "query is $query";
	   `$query`;
	   #removing intermediate outputs
	    my $file = $self->param('output_dir') . "/" . $sam;
	   `rm $file.sam $file.bam`;
	 }
	#generating flagstats values
	my $sorted = $self->param('output_dir') . "/" . $sam . ".sorted.bam";
	 say "generating flagstats values for file $sam";
	my $mapped = "samtools flagstat $sorted | awk -F \"[(|\%]\"  'NR == 5 {print \$2}'";
	$mapped = `$mapped`;
	say "mapped is $mapped";
	$sam = "basename " . $sorted ." .sorted.bam";
	$sam = `$sam`;
	chomp($sam);chomp($mapped);
	$sam = $sam . "\t" . $mapped;
	say "Sam is $sam";
	$self->param('alignment', $sam);
	
	 
}

sub write_output{
	 my ($self) = @_;
         my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -user   => $ENV{GBUSER},
        -dbname => $self->param('pipe_db'),
        -host   => $ENV{GBS1},
        -port   => $ENV{GBP1},
        -pass   => $ENV{GBPASS},
        -driver => $ENV{GBDRIVER},
    );
     my @report;
     my $sth = $dba->dbc->prepare("insert into alignment_stats values (?,?,?,?)");
	  my $stats = $self->param('output_dir') . "/flagstats.txt"; 
	  if ($self->param('alignment') =~ m/N\/A/){
	    $self->throw("Check alignment file. Percentage mapped cannot be empty");
	  }
	  else{
	    open ID, (">>$stats");
	    say "file opened for writing flagstat";
	    print ID $self->param('alignment'), "\n";
	    @report = split(/\t/,$self->param('alignment'));
            $sth->bind_param(1,$self->param('accession'));
            $sth->bind_param(2,$report[0]);
            $sth->bind_param(3,$report[1]);
            $sth->bind_param(4,0);
	  }
	  if ($sth->execute){
            say "Mapping stats stored for alignment file $report[0]";
          }
          else{
            $self->throw("Failed to store mapping stats for alignment $report[0]");
          }
	  my $dir = $self->param('output_dir') . "/GCA_*";
	
}
1;

