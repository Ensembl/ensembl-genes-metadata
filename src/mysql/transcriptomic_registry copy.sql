
CREATE DATABASE IF NOT EXISTS transcriptomic_hackathon
    CHARACTER SET utf8mb4
    COLLATE utf8mb4_unicode_ci;

USE transcriptomic_hackathon;

-- ============================================
-- 1. STUDY 
-- ============================================

CREATE TABLE IF NOT EXISTS study (
    study_id INT NOT NULL AUTO_INCREMENT,
    study_accession VARCHAR(255) NOT NULL,
    center_name VARCHAR(255),
    study_title TEXT,

    PRIMARY KEY (study_id),
    UNIQUE KEY uq_study_accession (study_accession)
) ENGINE=InnoDB;


-- ============================================
-- 2. BIOSAMPLE
-- Biological characteristics from BioSamples
-- with ENA fallback when needed
-- ============================================

CREATE TABLE IF NOT EXISTS biosample (
    sample_accession VARCHAR(255) NOT NULL,

    tissue VARCHAR(255),
    age VARCHAR(255),
    dev_stage VARCHAR(255),
    sex VARCHAR(100),
    organ VARCHAR(255),
    sub_tissue VARCHAR(255),

    sample_title TEXT,
    sample_description TEXT,

    experiment_title TEXT,
    experiment_alias VARCHAR(255),
    description TEXT,

    scientific_name VARCHAR(255),
    strain VARCHAR(255),
    isolate VARCHAR(255),
    cell_type VARCHAR(255),
    cell_line VARCHAR(255),
    cultivar VARCHAR(255),
    ecotype VARCHAR(255),

    PRIMARY KEY (sample_accession)
) ENGINE=InnoDB;


-- ============================================
-- 3. SAMPLE CONDITION
-- Biological / experimental conditions
-- ============================================

CREATE TABLE IF NOT EXISTS sample_condition (
    condition_id INT NOT NULL AUTO_INCREMENT,
    sample_accession VARCHAR(255) NOT NULL,

    disease VARCHAR(255),
    infection_status VARCHAR(100),
    treatment TEXT,

    PRIMARY KEY (condition_id),
    UNIQUE KEY uq_condition_sample (sample_accession),

    CONSTRAINT fk_condition_biosample
        FOREIGN KEY (sample_accession)
        REFERENCES biosample(sample_accession)
        ON UPDATE CASCADE
        ON DELETE RESTRICT
) ENGINE=InnoDB;


-- ============================================
-- 4. RUN
-- ENA sequencing metadata and pipeline status
-- ============================================

CREATE TABLE IF NOT EXISTS run (
    run_id INT NOT NULL AUTO_INCREMENT,
    run_accession VARCHAR(255) NOT NULL,

    sample_accession VARCHAR(255) NOT NULL,
    study_accession VARCHAR(255) NOT NULL,

    -- Taxonomy
    taxon_id INT,

    -- ENA metadata
    first_created DATE,
    instrument_platform VARCHAR(255),
    instrument_model VARCHAR(255),

    -- Library metadata
    library_strategy VARCHAR(255),
    library_source VARCHAR(255),
    library_selection VARCHAR(255),
    library_layout VARCHAR(100),
    library_name VARCHAR(255),

    -- FASTQ files (R1 and R2)
    fastq_1 TEXT,
    fastq_2 TEXT,
    md5_1 CHAR(32),
    md5_2 CHAR(32),


    -- Pipeline progress
    processing_status ENUM(
        'METADATA_RETRIEVAL_FAILED',
        'STAR_ISSUE',
        'NOT_CHECKED',
        'FILE_ISSUE',
        'QC_PASS',
        'QC_PROCESSED',
        'QC_FAIL',
        'ALIGNED',
        'DOWNLOAD_FAILED'
    ) NOT NULL DEFAULT 'PENDING',

    -- Overall quality-control result
    qc_result ENUM(
        'NOT_CHECKED',
        'PASS',
        'WARN',
        'FAIL'
    ) NOT NULL DEFAULT 'NOT_CHECKED',

    PRIMARY KEY (run_id),
    UNIQUE KEY uq_run_accession (run_accession),

    KEY idx_run_sample (sample_accession),
    KEY idx_run_study (study_accession),
    KEY idx_run_taxon (taxon_id),
    KEY idx_run_processing (processing_status),

    CONSTRAINT fk_run_biosample
        FOREIGN KEY (sample_accession)
        REFERENCES biosample(sample_accession)
        ON UPDATE CASCADE
        ON DELETE RESTRICT,

    CONSTRAINT fk_run_study
        FOREIGN KEY (study_accession)
        REFERENCES study(study_accession)
        ON UPDATE CASCADE
        ON DELETE RESTRICT
) ENGINE=InnoDB;


-- ============================================
-- 5. DATA FILES
-- Individual FASTQ files and FastQC results
-- ============================================

CREATE TABLE IF NOT EXISTS data_files (
    file_id INT NOT NULL AUTO_INCREMENT,
    run_id INT NOT NULL,

    file_name VARCHAR(255) NOT NULL,
    file_url TEXT NOT NULL,
    md5 VARCHAR(255) NOT NULL,

    -- FastQC module results
    basic_statistics ENUM('PASS','FAIL','WARN'),
    per_base_sequence_quality ENUM('PASS','FAIL','WARN'),
    per_tile_sequence_quality ENUM('PASS','FAIL','WARN'),
    per_sequence_quality_scores ENUM('PASS','FAIL','WARN'),
    per_base_sequence_content ENUM('PASS','FAIL','WARN'),
    per_sequence_gc_content ENUM('PASS','FAIL','WARN'),
    per_base_n_content ENUM('PASS','FAIL','WARN'),
    sequence_length_distribution ENUM('PASS','FAIL','WARN'),
    sequence_duplication_levels ENUM('PASS','FAIL','WARN'),
    overrepresented_sequences ENUM('PASS','FAIL','WARN'),
    adapter_content ENUM('PASS','FAIL','WARN'),

    -- FastQC statistics
    total_sequences BIGINT UNSIGNED,
    gc_content INT,
    sequence_length INT,

    PRIMARY KEY (file_id),
    UNIQUE KEY uq_file_name (file_name),
    KEY idx_file_run (run_id),

    CONSTRAINT fk_data_files_run
        FOREIGN KEY (run_id)
        REFERENCES run(run_id)
        ON UPDATE CASCADE
        ON DELETE RESTRICT
) ENGINE=InnoDB;


-- ============================================
-- 6. ALIGN
-- STAR alignment results
-- ============================================

CREATE TABLE IF NOT EXISTS align (
    align_id INT NOT NULL AUTO_INCREMENT,
    run_id INT NOT NULL,

    assembly_accession VARCHAR(255),

    uniquely_mapped_reads_percentage FLOAT,
    percentage_reads_mapped_to_multiple_loci FLOAT,
    percentage_reads_unmapped_too_short FLOAT,

    PRIMARY KEY (align_id),
    UNIQUE KEY uq_align_run (run_id),

    CONSTRAINT fk_align_run
        FOREIGN KEY (run_id)
        REFERENCES run(run_id)
        ON UPDATE CASCADE
        ON DELETE RESTRICT
) ENGINE=InnoDB;


-- ============================================
-- 7. META
-- Taxon-level metadata refresh tracking
-- ============================================

CREATE TABLE IF NOT EXISTS meta (
    meta_id INT NOT NULL AUTO_INCREMENT,
    taxon_id INT NOT NULL,
    last_check DATE,

    PRIMARY KEY (meta_id),
    KEY idx_meta_taxon (taxon_id)
) ENGINE=InnoDB;

-- ============================================
-- 8. TISSUE PREDICTION
-- Tissue prediction results for each sample, 
-- including original tissue, predicted tissue, and associated scores.
-- ============================================
CREATE TABLE tissue_prediction (
    prediction_id BIGINT AUTO_INCREMENT PRIMARY KEY,
    sample_accession VARCHAR(255) NOT NULL,

    original_tissue TEXT,
    predicted_tissue VARCHAR(255),
    ontology_id VARCHAR(100),
    ontology_source VARCHAR(50),

    semantic_similarity FLOAT,
    lexical_similarity FLOAT,
    context_score FLOAT,
    taxonomy_score FLOAT,
    combined_score FLOAT,

    prediction_method VARCHAR(100),
    model_version VARCHAR(255),

    validation_status ENUM(
        'AUTO_ACCEPTED',
        'REVIEW',
        'REJECTED',
        'UNMATCHED'
    ) NOT NULL DEFAULT 'REVIEW',

    explanation TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,

    INDEX idx_prediction_sample (sample_accession),
    INDEX idx_validation_status (validation_status)
) ENGINE=InnoDB
  DEFAULT CHARSET=utf8mb4;