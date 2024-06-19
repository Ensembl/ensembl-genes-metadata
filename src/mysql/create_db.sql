create database transcriptomic_hackathon;
use transcriptomic_hackathon;
create table data_files (
	file_id int not null auto_increment,
	run_id int not null,
	file_name varchar(255) not null,
	file_url varchar(255) not null,
	md5 varchar(255) not null,
	basic_statistics SET("PASS","FAIL","WARN"),
	per_base_sequence_quality SET("PASS","FAIL","WARN"),
	per_tile_sequence_quality SET("PASS","FAIL","WARN"),
	per_sequence_quality_scores SET("PASS","FAIL","WARN"),
	per_base_sequence_content SET("PASS","FAIL","WARN"),
	per_sequence_gc_content SET("PASS","FAIL","WARN"),
	per_base_n_content SET("PASS","FAIL","WARN"),
	sequence_length_distribution SET("PASS","FAIL","WARN"),
	sequence_duplication_levels SET("PASS","FAIL","WARN"),
	overrepresented_sequences SET("PASS","FAIL","WARN"),
	adapter_content SET("PASS","FAIL","WARN"),
	total_sequences int,
	gc_content int,
	sequence_length int,
	primary	key (file_id)
);
ALTER TABLE data_files ADD UNIQUE (file_name);
create table study (
	study_id int not null auto_increment,
	study_accession varchar(255) not null,
	center_name varchar(255),
	primary key (study_id)
);
ALTER TABLE study ADD UNIQUE (study_accession);

create table align (
	align_id int not null auto_increment,
	run_id int not null,
	assembly_accession varchar(255),
	uniquely_mapped_reads_percentage float,
	percentage_reads_mapped_to_multiple_loci float,
	percentage_reads_unmapped_too_short float,
	primary key (align_id)
);
ALTER TABLE align ADD UNIQUE (run_id);
create table run (
	run_id int not null auto_increment,
	taxon_id int not null,
	run_accession varchar(255) not null,
	qc_status SET("NOT_CHECKED","FILE_ISSUE","QC_PASS","QC_FAIL","ALIGNED"),
	sample_accession varchar(255) not null,
	study_accession varchar(255) not null,
	read_type varchar(255),
	platform varchar(255),
	paired bool not null,
	experiment varchar(255),
	run_description text,
	library_name varchar(255),
	library_selection varchar(255),
	tissue varchar(255),
	cell_line varchar(255),
	cell_type varchar(255),
	strain varchar(255),
	primary key (run_id)
);
ALTER TABLE run ADD UNIQUE (run_accession);

create table meta (
	meta_id int not null auto_increment,
	taxon_id int not null,
	last_check TIMESTAMP,
	primary key (meta_id)
);

