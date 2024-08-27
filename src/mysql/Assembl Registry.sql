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
  stable_id_prefix smallint(5),
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
  species_prefix varchar(20),
  clade varchar(25),
  PRIMARY KEY (`species_id`),
  FOREIGN KEY (`clade`) REFERENCES clade_settings(`clade`),
  FOREIGN KEY (`lowest_taxon_id`) REFERENCES assembly(`lowest_taxon_id`)
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

DROP TABLE IF EXISTS custom_group

CREATE TABLE custom_group (
  group_id int NOT NULL AUTO_INCREMENT,
  group_name varchar(30) NOT NULL,
  group_type ENUM('taxon', 'assembly') NOT NULL,
  item varchar(30) NOT NULL,
  PRIMARY KEY (`group_id`),
  CONSTRAINT group_item UNIQUE (group_name, item)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


DROP TABLE IF EXISTS group_assembly ;

CREATE TABLE group_assembly (
  group_assembly_id int NOT NULL AUTO_INCREMENT,
  assembly_id int NOT NULL,
  group_name varchar(30) NOT NULL,
  PRIMARY KEY (`group_assembly_id`),
  FOREIGN KEY (`assembly_id`) REFERENCES assembly(`assembly_id`),
  CONSTRAINT group_asm UNIQUE (assembly_id, group_name)


DROP TABLE IF EXISTS stable_space

CREATE TABLE stable_space (
  stable_space_id int NOT NULL,
  stable_space_start bigint(20) NOT NULL,
  stable_space_end bigint(20) NOT NULL,
  PRIMARY KEY (`stable_space_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

DROP TABLE IF EXISTS clade_settings

CREATE TABLE clade_settings (
  clade varchar(55) NOT NULL,
  repbase_lib varchar(55),
  busco_lineage varchar(55),
  species_division varchar(55) NOT NULL,
  uniprot_set varchar(225),
  orthodb_set varchar(225),
  rfam_search_term varchar(55),
  PRIMARY KEY (`clade`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


DROP TABLE IF EXISTS genebuild

CREATE TABLE genebuild (
  genebuild_id int NOT NULL AUTO_INCREMENT,
  assembly_id int NOT NULL,
  status ENUM('in_progress', 'completed', 'handover'),
  genebuilder_id int(5),
  annotation_source ENUM('ensembl', 'refseq', 'community', 'wormbase', 'flybase'),
  annotation_method ENUM('main', 'anno', 'braker', 'imported'),
  date_started datetime,
  date_complete datetime,
  release_type ENUM('main', 'rapid', 'draft'),
  release_date date,
  release_version int(5),
  PRIMARY KEY (`genebuild_id`),
  FOREIGN KEY (`assembly_id`) REFERENCES assembly(`assembly_id`),
  FOREIGN KEY (`genebuilder_id`) REFERENCES genebuilder(`genebuilder_id`)
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


DROP TABLE IF EXISTS genebuilder

CREATE TABLE genebuilder (
  genebuilder_id int NOT NULL,
  genebuilder_name varchar(20) NOT NULL,
  PRIMARY KEY (`genebuilder_id`)
) ENGINE=MyISAM AUTO_INCREMENT=1 DEFAULT CHARSET=latin1;