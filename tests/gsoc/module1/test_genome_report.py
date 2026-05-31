"""
Tests for busco_utils.py
Run with: pytest tests/gsoc/module1/test_genome_report.py -v
"""

import pytest
from metadata_app.backend.app.services.gsoc.module1.busco_utils import (
    parse_busco_string,
    get_busco_complete,
    busco_quality_label,
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
        assert all(v is None for v in result.values())

    def test_none_input_returns_none_values(self):
        result = parse_busco_string(None)
        assert all(v is None for v in result.values())

    def test_malformed_string_returns_none_values(self):
        result = parse_busco_string("not_a_busco_string")
        assert all(v is None for v in result.values())

    def test_integer_complete_value(self):
        result = parse_busco_string("C:100%[S:98%,D:2%],F:0%,M:0%,n:300")
        assert result["complete"] == 100.0
        assert result["n_genes"] == 300


class TestGetBuscoComplete:
    """Tests for get_busco_complete()"""

    def test_returns_float(self):
        result = get_busco_complete("C:94.3%[S:91.2%,D:3.1%],F:2.1%,M:3.6%,n:255")
        assert result == 94.3
        assert isinstance(result, float)

    def test_returns_none_for_empty(self):
        assert get_busco_complete("") is None

    def test_returns_none_for_none(self):
        assert get_busco_complete(None) is None


class TestBuscoQualityLabel:
    """Tests for busco_quality_label()"""

    def test_excellent(self):
        assert busco_quality_label(95.0) == "Excellent"
        assert busco_quality_label(100.0) == "Excellent"

    def test_good(self):
        assert busco_quality_label(85.0) == "Good"
        assert busco_quality_label(94.9) == "Good"

    def test_moderate(self):
        assert busco_quality_label(70.0) == "Moderate"
        assert busco_quality_label(84.9) == "Moderate"

    def test_poor(self):
        assert busco_quality_label(69.9) == "Poor"
        assert busco_quality_label(0.0) == "Poor"

    def test_none_returns_unknown(self):
        assert busco_quality_label(None) == "Unknown"
