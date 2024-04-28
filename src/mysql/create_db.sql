create database transcriptomic_hackathon;
use transcriptomic_hackathon;
create table meta (
	meta_id int not null auto_increment,
	taxon_id int not null,
	last_check DATE,
	primary key (meta_id)
);
create table data_files (
	file_id int not null auto_increment,
	run_id int not null,
	name varchar(255) not null,
	url varchar(255) not null,
	md5 varchar(255) not null,
	basic_statistics SET("PASS","FAIL","WARN"),
	per_base_sequence_quality SET("PASS","FAIL","WARN"),
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
	read_count int,
	read_length int,
	primary	key (file_id)
);
create table study (
	study_id int not null auto_increment,
	study_accession varchar(255) not null,
	center_name varchar(255),
	primary key (study_id)
);
create table align (
	align_id int not null auto_increment,
	run_id int not null,
	assembly_accession varchar(255),
	percent_mapped int,
	primary key (align_id)
);
create table run (
	run_id int not null auto_increment,
	taxon_id int not null,
	run_accession varchar(255) not null,
	species_taxon_id int not null,
	genus_taxon_id int,
	qc_status SET("NOT_CHECKED","QC_PASS","QC_FAIL","ALIGNED"),
	sample_accession varchar(255) not null,
	study_accession varchar(255) not null,
	read_type varchar(255),
	platform varchar(255),
	paired bool not null,
	experiment varchar(255),
	description text,
	library_name varchar(255),
	library_selection varchar(255),
	tissue varchar(255),
	cell_line varchar(255),
	cell_type varchar(255),
	strain varchar(255),
	primary key (run_id)
);

