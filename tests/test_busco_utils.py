# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Tests for BUSCO summary parsing."""

from metadata_app.backend.app.services.busco_utils import parse_busco_string


def test_parse_busco_string_extracts_all_standard_components():
    result = parse_busco_string("C:98.2%[S:97.5%,D:0.7%],F:0.8%,M:1.0%,n:255")

    assert result == {
        "complete": 98.2,
        "single_copy": 97.5,
        "duplicated": 0.7,
        "fragmented": 0.8,
        "missing": 1.0,
        "n_genes": 255,
    }


def test_parse_busco_string_accepts_integer_percentages():
    result = parse_busco_string("C:100%[S:98%,D:2%],F:0%,M:0%,n:300")

    assert result["complete"] == 100.0
    assert result["single_copy"] == 98.0
    assert result["duplicated"] == 2.0
    assert result["n_genes"] == 300


def test_parse_busco_string_keeps_missing_components_as_none():
    result = parse_busco_string("C:95.0%,F:2.0%,n:100")

    assert result["complete"] == 95.0
    assert result["fragmented"] == 2.0
    assert result["single_copy"] is None
    assert result["duplicated"] is None
    assert result["missing"] is None
    assert result["n_genes"] == 100


def test_parse_busco_string_handles_empty_and_malformed_values():
    assert all(value is None for value in parse_busco_string("").values())
    assert all(value is None for value in parse_busco_string("not a BUSCO string").values())
    assert all(value is None for value in parse_busco_string(None).values())
