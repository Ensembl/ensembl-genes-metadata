USE gb_assembly_metadata ;
DROP TABLE IF EXISTS assembly

CREATE TABLE assembly (
  assembly_id int NOT NULL AUTO_INCREMENT,
  lowest_taxon_id int(15) NOT NULL,
  gca_chain varchar(30) NOT NULL,
  gca_version smallint(3) NOT NULL,
  is_current varchar(30),
  asm_type ENUM('haploid', 'diploid', 'haploid-with-alt-loci', 'unresolved-diploid', 'alternate-pseudohaplotype'),
  asm_level ENUM('Scaffold', 'Contig', 'Chromosome', 'Complete genome'),
  asm_name varchar(225) NOT NULL,
  refseq_accession varchar(30),
  release_date date,
  submitter varchar(225),
  PRIMARY KEY (`assembly_id`),
  CONSTRAINT gca_accession UNIQUE (gca_chain, gca_version)
) ENGINE=MyISAM AUTO_INCREMENT=1 DEFAULT CHARSET=latin1;


DROP TABLE IF EXISTS assembly_metrics

CREATE TABLE assembly_metrics (
  asm_metrics_id int NOT NULL AUTO_INCREMENT,
  assembly_id int NOT NULL,
  metrics_name varchar(50) NOT NULL,
  metrics_value varchar(225) NOT NULL,
  PRIMARY KEY (`asm_metrics_id`),
  FOREIGN KEY (`assembly_id`) REFERENCES assembly(`assembly_id`),
  CONSTRAINT metric_record UNIQUE (assembly_id, metrics_name, metrics_value)
) ENGINE=MyISAM AUTO_INCREMENT=1 DEFAULT CHARSET=latin1;


DROP TABLE IF EXISTS species

CREATE TABLE species (
  species_id int NOT NULL AUTO_INCREMENT,
  lowest_taxon_id int(15) NOT NULL UNIQUE,
  species_taxon_id int(15) NOT NULL,
  scientific_name varchar(225) NOT NULL,
  common_name varchar(225),
  parlance_name varchar(225),
  clade varchar(25),
  PRIMARY KEY (`species_id`),
  FOREIGN KEY (`clade`) REFERENCES clade_settings(`clade`),
  FOREIGN KEY (`lowest_taxon_id`) REFERENCES assembly(`lowest_taxon_id`)
) ENGINE=MyISAM AUTO_INCREMENT=1 DEFAULT CHARSET=latin1;

DROP TABLE IF EXISTS species_prefix

CREATE TABLE species_prefix (
	prefix_id int(10) NOT NULL AUTO_INCREMENT,
	lowest_taxon_id int(15) NOT NULL,
	prefix varchar(7) NOT NULL,
	PRIMARY KEY (`prefix_id`),
	FOREIGN KEY (`lowest_taxon_id`) REFERENCES assembly(`lowest_taxon_id`),
	CONSTRAINT species_prefix UNIQUE (lowest_taxon_id, prefix)
) ENGINE=MyISAM AUTO_INCREMENT=1 DEFAULT CHARSET=latin1;

DROP TABLE IF EXISTS taxonomy

CREATE TABLE taxonomy (
  taxon_clsf_id int NOT NULL AUTO_INCREMENT,
  lowest_taxon_id int(15) NOT NULL,
  taxon_class_id int(15) NOT NULL,
  taxon_class varchar(50) NOT NULL,
  PRIMARY KEY (`taxon_clsf_id`),
  FOREIGN KEY (`lowest_taxon_id`) REFERENCES assembly(`lowest_taxon_id`),
  CONSTRAINT taxon_classification UNIQUE (lowest_taxon_id, taxon_class_id, taxon_class ),
  CONSTRAINT unique_rank_per_taxon UNIQUE (lowest_taxon_id, taxon_class)
) ENGINE=MyISAM AUTO_INCREMENT=1 DEFAULT CHARSET=latin1;


DROP TABLE IF EXISTS taxonomy_names

CREATE TABLE taxonomy_name (
  taxon_class_id int(15) NOT NULL,
  taxon_class_name varchar(50) NOT NULL,
  FOREIGN KEY (`taxon_class_id`) REFERENCES taxonomy(`taxon_class_id`),
  CONSTRAINT taxon_id_name UNIQUE (taxon_class_id, taxon_class_name)
) ENGINE=MyISAM AUTO_INCREMENT=1 DEFAULT CHARSET=latin1;


DROP TABLE IF EXISTS organism

CREATE TABLE organism (
  organism_id int NOT NULL AUTO_INCREMENT,
  assembly_id int NOT NULL UNIQUE,
  biosample_id varchar(225),
  bioproject_id varchar(225),
  dtol_id varchar(30),
  infra_type ENUM ('', 'strain', 'breed', 'cultivar' , 'ecotype', 'isolate' ),
  infra_name varchar(225),
  PRIMARY KEY (`organism_id`),
  FOREIGN KEY (`assembly_id`) REFERENCES assembly(`assembly_id`),
) ENGINE=MyISAM AUTO_INCREMENT=1 DEFAULT CHARSET=latin1;


DROP TABLE IF EXISTS bioproject ;

CREATE TABLE bioproject (
  lineage_id int NOT NULL AUTO_INCREMENT,
  assembly_id int NOT NULL,
  bioproject_id varchar(50) NOT NULL,
  PRIMARY KEY (`lineage_id`),
  FOREIGN KEY (`assembly_id`) REFERENCES assembly(`assembly_id`),
  CONSTRAINT lineage UNIQUE (assembly_id, bioproject_id)
) ENGINE=MyISAM AUTO_INCREMENT=1 DEFAULT CHARSET=latin1;

DROP TABLE IF EXISTS main_bioproject

CREATE TABLE main_bioproject (
  bioproject_id varchar(50) NOT NULL,
  bioproject_name varchar(50) NOT NULL,
  FOREIGN KEY (`bioproject_id`) REFERENCES bioproject(`bioproject_id`),
  CONSTRAINT unique_bioproject UNIQUE (bioproject_id),
  CONSTRAINT unique_id_name UNIQUE (bioproject_id, bioproject_name)
) ENGINE=MyISAM AUTO_INCREMENT=1 DEFAULT CHARSET=latin1;


DROP TABLE IF EXISTS custom_group

CREATE TABLE custom_group (
  group_id int NOT NULL AUTO_INCREMENT,
  group_name varchar(30) NOT NULL,
  group_type ENUM('taxon', 'assembly') NOT NULL,
  item varchar(30) NOT NULL,
  PRIMARY KEY (`group_id`),
  CONSTRAINT group_item UNIQUE (group_name, item)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

DROP TABLE IF EXISTS stable_space

CREATE TABLE stable_space (
  stable_space_id int NOT NULL,
  stable_space_start bigint(20) NOT NULL,
  stable_space_end bigint(20) NOT NULL,
  PRIMARY KEY (`stable_space_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

DROP TABLE IF EXISTS species_spaces

CREATE TABLE species_spaces (
	lowest_taxon_id int(15) NOT NULL,
	assembly_id int(11) NOT NULL,
	stable_space_id int(10) NOT NULl,
	FOREIGN KEY (`lowest_taxon_id`) REFERENCES assembly(`lowest_taxon_id`),
	FOREIGN KEY (`assembly_id`) REFERENCES assembly(`assembly_id`),
	FOREIGN KEY (`stable_space_id`) REFERENCES stable_space(`assembly_id`),
	CONSTRAINT species_space_assign UNIQUE (lowest_taxon_id, assembly_id, stable_space_id)
) ENGINE=MyISAM AUTO_INCREMENT=1 DEFAULT CHARSET=latin1;

DROP TABLE IF EXISTS genebuilder

CREATE TABLE genebuilder (
  genebuilder_id int NOT NULL,
  genebuilder varchar(20),
  PRIMARY KEY (`genebuilder_id`)
) ENGINE=MyISAM AUTO_INCREMENT=1 DEFAULT CHARSET=latin1;

DROP TABLE IF EXISTS genebuild_status

CREATE TABLE genebuild_status (
  genebuild_id int NOT NULL AUTO_INCREMENT,
  assembly_id int NOT NULL,
  gca_accession VARCHAR(20) NOT NULL, 
  gb_status ENUM('in_progress', 'completed', 'handed_over'),
  last_attempt int(2),
  genebuilder varchar(20) NOT NULL,
  annotation_source ENUM('ensembl', 'external','import_refseq', 'import_community', 'import_wormbase', 'import_flybase', 'import_genbank', 'import_noninsdc'),
  annotation_method ENUM('pending','full_genebuild', 'anno', 'braker', 'projection_build', 'mixed_strategy_build','import', 'external_annotation_import'),
  date_started date NOT NULL,
  date_completed date NULL,
  date_completed_beta date NULL,
  release_type ENUM('main', 'beta', 'not_available'),
  release_date date NULL,
  release_date_beta date NULL,
  release_version int(5) NULL,
  release_version_beta int(5) NULL,
  PRIMARY KEY (`genebuild_id`),
  FOREIGN KEY (`assembly_id`) REFERENCES assembly(`assembly_id`),
  FOREIGN KEY (`genebuilder`) REFERENCES genebuilder(`genebuilder`)
  CONSTRAINT species_space_assign UNIQUE (lowest_taxon_id, assembly_id, stable_space_id)
) ENGINE=MyISAM AUTO_INCREMENT=1 DEFAULT CHARSET=latin1;


DROP TABLE IF EXISTS genebuild_metrics

CREATE TABLE genebuild_metrics (
  gb_metrics_id int NOT NULL AUTO_INCREMENT,
  genebuild_id int NOT NULL,
  metrics_name varchar(50),
  metrics_value varchar(225),
  PRIMARY KEY (`gb_metrics_id`),
  FOREIGN KEY (`genebuild_id`) REFERENCES assembly(`genebuild_id`)
) ENGINE=MyISAM AUTO_INCREMENT=1 DEFAULT CHARSET=latin1;


