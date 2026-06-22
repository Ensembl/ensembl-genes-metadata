"""
Tests for genome_report.py and busco_utils.py
Run with: pytest tests/gsoc/module1/test_genome_report.py -v
"""

# pylint: disable=redefined-outer-name,missing-function-docstring

import pytest
import pandas as pd
from metadata_app.backend.app.services.gsoc.module1.busco_utils import (
    parse_busco_string,
    get_busco_complete,
    busco_quality_label,
)
from metadata_app.backend.app.services.gsoc.module1.genome_report import (
    extract_genome_report,
    GenomeReport,
)


@pytest.fixture
def sample_anno_wide() -> pd.DataFrame:
    """
    Synthetic anno_wide DataFrame mirroring annotations_service.generate_tables() output.
    No database connection required.
    """
    return pd.DataFrame(
        [
            {
                "gca": "GCA_000001405.29",
                "scientific_name": "Homo sapiens",
                "common_name": "human",
                "lowest_taxon_id": 9606,
                "internal_clade": "Primates",
                "annotation_method": "full_genebuild",
                "genebuilder": "genebuild_team",
                "gb_status": "live",
                "annotation_source": "ensembl",
                "bioproject_id": "PRJNA31257",
                "associated_project": "Human Genome Project",
                "release_date": "2023-01-15",
                "date_status_update": "2023-01-15",
                "last_genebuild_update": "2023-01-10",
                "protein_busco": "C:94.3%[S:91.2%,D:3.1%],F:2.1%,M:3.6%,n:255",
                "protein_busco_lineage": "primates_odb10",
                "protein_busco_version": "5.4.3",
                "assembly_busco": "C:98.1%[S:97.0%,D:1.1%],F:0.9%,M:1.0%,n:255",
                "assembly_busco_lineage": "primates_odb10",
                "assembly_busco_version": "5.4.3",
                "coding_genes": 20442,
                "latest_annotated": "Yes",
                "annotated_version": 29.0,
                "assembly_version": 29.0,
                "ftp": "https://ftp.ebi.ac.uk/pub/ensemblorganisms/Homo_sapiens/GCA_000001405.29/",
            },
            {
                "gca": "GCA_000002285.2",
                "scientific_name": "Gallus gallus",
                "common_name": "chicken",
                "lowest_taxon_id": 9031,
                "internal_clade": "Aves",
                "annotation_method": "full_genebuild",
                "genebuilder": "genebuild_team",
                "gb_status": "live",
                "annotation_source": "ensembl",
                "bioproject_id": "PRJNA10808",
                "associated_project": "Chicken Genome Project",
                "release_date": "2022-06-01",
                "date_status_update": "2022-06-01",
                "last_genebuild_update": "2022-05-28",
                "protein_busco": "C:85.0%[S:83.0%,D:2.0%],F:5.0%,M:10.0%,n:255",
                "protein_busco_lineage": "aves_odb10",
                "protein_busco_version": "5.4.3",
                "assembly_busco": "C:90.0%[S:89.0%,D:1.0%],F:4.0%,M:6.0%,n:255",
                "assembly_busco_lineage": "aves_odb10",
                "assembly_busco_version": "5.4.3",
                "coding_genes": 16736,
                "latest_annotated": "Yes",
                "annotated_version": 2.0,
                "assembly_version": 2.0,
                "ftp": "https://ftp.ebi.ac.uk/pub/ensemblorganisms/Gallus_gallus/GCA_000002285.2/",
            },
            {
                "gca": "GCA_000003625.1",
                "scientific_name": "Danio rerio",
                "common_name": "zebrafish",
                "lowest_taxon_id": 7955,
                "internal_clade": "Teleostei",
                "annotation_method": "full_genebuild",
                "genebuilder": "genebuild_team",
                "gb_status": "live",
                "annotation_source": "ensembl",
                "bioproject_id": "PRJNA13922",
                "associated_project": "Zebrafish Genome Project",
                "release_date": "",
                "date_status_update": "",
                "last_genebuild_update": "",
                "protein_busco": "",
                "protein_busco_lineage": "",
                "protein_busco_version": "",
                "assembly_busco": "",
                "assembly_busco_lineage": "",
                "assembly_busco_version": "",
                "coding_genes": "",
                "latest_annotated": "No",
                "annotated_version": 1.0,
                "assembly_version": 2.0,
                "ftp": None,
            },
        ]
    )


class TestParseBuscoString:
    """Tests for parse_busco_string()"""

    def test_full_valid_string(self):
        result = parse_busco_string("C:94.3%[S:91.2%,D:3.1%],F:2.1%,M:3.6%,n:255")
        assert result["complete"] == 94.3
        assert result["single_copy"] == 91.2
        assert result["duplicated"] == 3.1
        assert result["fragmented"] == 2.1
        assert result["missing"] == 3.6
        assert result["n_genes"] == 255

    def test_empty_string_returns_none_values(self):
        result = parse_busco_string("")
        core_keys = {k: v for k, v in result.items() if k != "extra"}
        assert all(v is None for v in core_keys.values())
        assert not result["extra"]

    def test_none_input_returns_none_values(self):
        result = parse_busco_string(None)
        core_keys = {k: v for k, v in result.items() if k != "extra"}
        assert all(v is None for v in core_keys.values())
        assert not result["extra"]

    def test_malformed_string_returns_none_values(self):
        result = parse_busco_string("not_a_busco_string")
        core_keys = {k: v for k, v in result.items() if k != "extra"}
        assert all(v is None for v in core_keys.values())
        assert not result["extra"]

    def test_integer_complete_value(self):
        result = parse_busco_string("C:100%[S:98%,D:2%],F:0%,M:0%,n:300")
        assert result["complete"] == 100.0
        assert result["n_genes"] == 300

    def test_extra_field_captured(self):
        result = parse_busco_string(
            "C:97.8%[S:93.6%,D:4.2%],F:0.3%,M:1.9%,n:7207,E:4.5%"
        )
        assert result["complete"] == 97.8
        assert result["n_genes"] == 7207
        assert result["extra"] == {"E": 4.5}

    def test_no_extra_field_gives_empty_dict(self):
        result = parse_busco_string("C:94.3%[S:91.2%,D:3.1%],F:2.1%,M:3.6%,n:255")
        assert not result["extra"]

    def test_whitespace_is_stripped(self):
        result = parse_busco_string("  C:94.3%[S:91.2%,D:3.1%],F:2.1%,M:3.6%,n:255  ")
        assert result["complete"] == 94.3
        assert result["n_genes"] == 255

    def test_lowercase_letters_are_accepted(self):
        result = parse_busco_string("c:94.3%[s:91.2%,d:3.1%],f:2.1%,m:3.6%,n:255")
        assert result["complete"] == 94.3
        assert result["single_copy"] == 91.2
        assert result["n_genes"] == 255


class TestGetBuscoComplete:
    """Tests for get_busco_complete()"""

    def test_returns_float(self):
        result = get_busco_complete("C:94.3%[S:91.2%,D:3.1%],F:2.1%,M:3.6%,n:255")
        assert result == 94.3

    def test_returns_none_for_empty(self):
        result = get_busco_complete("")
        assert result is None

    def test_returns_none_for_none(self):
        result = get_busco_complete(None)
        assert result is None


class TestBuscoQualityLabel:
    """Tests for busco_quality_label()"""

    def test_excellent(self):
        assert busco_quality_label(96.0) == "Excellent"
        assert busco_quality_label(95.0) == "Excellent"

    def test_good(self):
        assert busco_quality_label(90.0) == "Good"
        assert busco_quality_label(85.0) == "Good"

    def test_moderate(self):
        assert busco_quality_label(75.0) == "Moderate"
        assert busco_quality_label(70.0) == "Moderate"

    def test_poor(self):
        assert busco_quality_label(50.0) == "Poor"
        assert busco_quality_label(0.0) == "Poor"

    def test_none_returns_unknown(self):
        assert busco_quality_label(None) == "Unknown"


class TestExtractGenomeReport:
    """Tests for extract_genome_report()"""

    def test_returns_genome_report_instance(self, sample_anno_wide):
        report = extract_genome_report("GCA_000001405.29", sample_anno_wide)
        assert isinstance(report, GenomeReport)

    def test_correct_gca_extracted(self, sample_anno_wide):
        report = extract_genome_report("GCA_000001405.29", sample_anno_wide)
        assert report.gca == "GCA_000001405.29"

    def test_scientific_name(self, sample_anno_wide):
        report = extract_genome_report("GCA_000001405.29", sample_anno_wide)
        assert report.scientific_name == "Homo sapiens"

    def test_protein_busco_complete_parsed(self, sample_anno_wide):
        report = extract_genome_report("GCA_000001405.29", sample_anno_wide)
        assert report.protein_busco_complete == 94.3

    def test_protein_busco_quality_label(self, sample_anno_wide):
        report = extract_genome_report("GCA_000001405.29", sample_anno_wide)
        assert report.protein_busco_quality == "Good"

    def test_assembly_busco_complete_parsed(self, sample_anno_wide):
        report = extract_genome_report("GCA_000001405.29", sample_anno_wide)
        assert report.assembly_busco_complete == 98.1

    def test_coding_genes(self, sample_anno_wide):
        report = extract_genome_report("GCA_000001405.29", sample_anno_wide)
        assert report.coding_genes == 20442

    def test_internal_clade(self, sample_anno_wide):
        report = extract_genome_report("GCA_000001405.29", sample_anno_wide)
        assert report.internal_clade == "Primates"

    def test_ftp_link(self, sample_anno_wide):
        report = extract_genome_report("GCA_000001405.29", sample_anno_wide)
        assert report.ftp.startswith("https://ftp.ebi.ac.uk")

    def test_missing_gca_raises_value_error(self, sample_anno_wide):
        with pytest.raises(ValueError, match="GCA_999999999.1"):
            extract_genome_report("GCA_999999999.1", sample_anno_wide)

    def test_empty_busco_fields_return_none(self, sample_anno_wide):
        report = extract_genome_report("GCA_000003625.1", sample_anno_wide)
        assert report.protein_busco_complete is None
        assert report.protein_busco_quality == "Unknown"

    def test_empty_coding_genes_returns_none(self, sample_anno_wide):
        report = extract_genome_report("GCA_000003625.1", sample_anno_wide)
        assert report.coding_genes is None

    def test_second_genome_correct(self, sample_anno_wide):
        report = extract_genome_report("GCA_000002285.2", sample_anno_wide)
        assert report.scientific_name == "Gallus gallus"
        assert report.protein_busco_complete == 85.0
        assert report.protein_busco_quality == "Good"
        assert report.coding_genes == 16736

    def test_clade_context_defaults_to_none(self, sample_anno_wide):
        report = extract_genome_report("GCA_000001405.29", sample_anno_wide)
        assert report.clade_median_busco is None
        assert report.clade_median_coding_genes is None
        assert report.clade_sample_size is None

    def test_protein_busco_extra_populated(self, sample_anno_wide):
        report = extract_genome_report("GCA_000001405.29", sample_anno_wide)
        assert not report.protein_busco_extra

    def test_assembly_busco_extra_populated(self, sample_anno_wide):
        report = extract_genome_report("GCA_000001405.29", sample_anno_wide)
        assert not report.assembly_busco_extra


class TestMultipleAnnotationRows:
    """
    Tests for the 'current annotation' selection rule when a single GCA
    has multiple genebuild_status rows (see code review 2026-06-18:
    extract_genome_report() previously used genome_rows.iloc[0] without
    any ordering rule, which could silently pick the wrong annotation).
    """

    @pytest.fixture
    def duplicate_gca_anno_wide(self) -> pd.DataFrame:
        """Two rows for the same GCA with different date_status_update values."""
        return pd.DataFrame(
            [
                {
                    "gca": "GCA_000001515.5",
                    "scientific_name": "Pan troglodytes",
                    "common_name": "chimpanzee",
                    "lowest_taxon_id": 9598,
                    "internal_clade": "Primates",
                    "annotation_method": "full_genebuild",
                    "genebuilder": "genebuild_team",
                    "gb_status": "superseded",
                    "annotation_source": "ensembl",
                    "bioproject_id": "PRJNA13184",
                    "associated_project": "",
                    "release_date": "2020-05-01",
                    "date_status_update": "2020-05-01",
                    "last_genebuild_update": "2020-04-20",
                    "protein_busco": "C:88.0%[S:85.0%,D:3.0%],F:5.0%,M:7.0%,n:255",
                    "protein_busco_lineage": "primates_odb10",
                    "protein_busco_version": "5.4.3",
                    "assembly_busco": "",
                    "assembly_busco_lineage": "",
                    "assembly_busco_version": "",
                    "coding_genes": 19000,
                    "latest_annotated": "No",
                    "annotated_version": 1.0,
                    "assembly_version": 5.0,
                    "ftp": "",
                },
                {
                    "gca": "GCA_000001515.5",
                    "scientific_name": "Pan troglodytes",
                    "common_name": "chimpanzee",
                    "lowest_taxon_id": 9598,
                    "internal_clade": "Primates",
                    "annotation_method": "full_genebuild",
                    "genebuilder": "genebuild_team",
                    "gb_status": "live",
                    "annotation_source": "ensembl",
                    "bioproject_id": "PRJNA13184",
                    "associated_project": "",
                    "release_date": "2024-09-01",
                    "date_status_update": "2024-09-01",
                    "last_genebuild_update": "2024-08-15",
                    "protein_busco": "C:96.5%[S:95.0%,D:1.5%],F:1.0%,M:2.5%,n:255",
                    "protein_busco_lineage": "primates_odb10",
                    "protein_busco_version": "5.4.3",
                    "assembly_busco": "",
                    "assembly_busco_lineage": "",
                    "assembly_busco_version": "",
                    "coding_genes": 20100,
                    "latest_annotated": "Yes",
                    "annotated_version": 2.0,
                    "assembly_version": 5.0,
                    "ftp": (
                        "https://ftp.ebi.ac.uk/pub/ensemblorganisms/"
                        "Pan_troglodytes/GCA_000001515.5/"
                    ),
                },
            ]
        )

    def test_most_recent_row_is_selected(self, duplicate_gca_anno_wide):
        report = extract_genome_report("GCA_000001515.5", duplicate_gca_anno_wide)
        assert report.gb_status == "live"
        assert report.protein_busco_complete == 96.5
        assert report.coding_genes == 20100

    def test_warning_logged_for_multiple_rows(self, duplicate_gca_anno_wide, caplog):
        with caplog.at_level("WARNING"):
            extract_genome_report("GCA_000001515.5", duplicate_gca_anno_wide)
        assert any(
            "Multiple annotation rows found" in record.message
            for record in caplog.records
        )

    def test_single_row_does_not_warn(self, sample_anno_wide, caplog):
        with caplog.at_level("WARNING"):
            extract_genome_report("GCA_000001405.29", sample_anno_wide)
        assert not any(
            "Multiple annotation rows found" in record.message
            for record in caplog.records
        )
