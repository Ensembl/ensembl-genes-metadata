from metadata_app.backend.app.services.busco_parser import (
    BuscoScores,
    parse_busco_string,
)


def test_parse_full_busco_string():
    value = "C:98.2%[S:97.5%,D:0.7%],F:0.8%,M:1.0%,n:255"

    scores = parse_busco_string(value)

    assert scores == BuscoScores(
        complete=98.2,
        single=97.5,
        duplicated=0.7,
        fragmented=0.8,
        missing=1.0,
        searched=255,
    )


def test_parse_completeness_only_busco_string():
    scores = parse_busco_string("C:91.4%")

    assert scores == BuscoScores(complete=91.4)


def test_parse_returns_empty_scores_for_none():
    assert parse_busco_string(None) == BuscoScores()


def test_parse_returns_empty_scores_for_blank_string():
    assert parse_busco_string("   ") == BuscoScores()


def test_parse_returns_empty_scores_for_malformed_string():
    assert parse_busco_string("not a busco value") == BuscoScores()
