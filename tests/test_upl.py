"""Tests for upl_finder_fast.upl module."""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pytest

from upl_finder_fast.upl import (
    _find_all,
    _load_upl_probes_tsv,
    _revcomp,
    find_upl_matches,
    load_upl_probes,
)


# ---------------------------------------------------------------------------
# _revcomp
# ---------------------------------------------------------------------------

def test_revcomp_basic():
    assert _revcomp("ACTG") == "CAGT"


def test_revcomp_n():
    assert _revcomp("ACNGT") == "ACNGT"


def test_revcomp_palindrome():
    seq = "ACGT"
    assert _revcomp(seq) == "ACGT"


def test_revcomp_lowercase():
    assert _revcomp("actg") == "CAGT"


# ---------------------------------------------------------------------------
# _find_all
# ---------------------------------------------------------------------------

def test_find_all_single_hit():
    assert _find_all("AACGTAA", "CGT") == [2]


def test_find_all_multiple_hits():
    assert _find_all("AAAA", "AA") == [0, 1, 2]


def test_find_all_no_hit():
    assert _find_all("AAAA", "CC") == []


def test_find_all_empty_sub():
    assert _find_all("AAAA", "") == []


# ---------------------------------------------------------------------------
# _load_upl_probes_tsv
# ---------------------------------------------------------------------------

def test_load_tsv_basic():
    text = "p1\tACTG\np2\tCCCC\n"
    result = _load_upl_probes_tsv(text)
    assert result == {"p1": "ACTG", "p2": "CCCC"}


def test_load_tsv_skip_comments_and_blank():
    text = "# comment\n\np1\tACTG\n"
    result = _load_upl_probes_tsv(text)
    assert result == {"p1": "ACTG"}


def test_load_tsv_whitespace_separated():
    text = "p1 ACTG\n"
    result = _load_upl_probes_tsv(text)
    assert result == {"p1": "ACTG"}


def test_load_tsv_normalizes_case():
    # _load_upl_probes_tsv returns the sequence as-is; normalization happens in load_upl_probes.
    text = "p1\tactg\n"
    result = _load_upl_probes_tsv(text)
    assert result["p1"] == "actg"


def test_load_upl_probes_normalizes_case(tmp_path):
    f = tmp_path / "probes.tsv"
    f.write_text("p1\tactg\n", encoding="utf-8")
    result = load_upl_probes(str(f))
    assert result["p1"] == "ACTG"


def test_load_tsv_invalid_line():
    text = "only_one_field\n"
    with pytest.raises(ValueError, match="expected"):
        _load_upl_probes_tsv(text)


# ---------------------------------------------------------------------------
# load_upl_probes
# ---------------------------------------------------------------------------

def test_load_upl_probes_tsv_file(tmp_path):
    f = tmp_path / "probes.tsv"
    f.write_text("p1\tACTG\np2\tGGGG\n", encoding="utf-8")
    result = load_upl_probes(str(f))
    assert result == {"p1": "ACTG", "p2": "GGGG"}


def test_load_upl_probes_txt_file(tmp_path):
    f = tmp_path / "probes.txt"
    f.write_text("p1\tACTG\n", encoding="utf-8")
    result = load_upl_probes(str(f))
    assert result == {"p1": "ACTG"}


def test_load_upl_probes_json_file(tmp_path):
    f = tmp_path / "probes.json"
    f.write_text(json.dumps({"p1": "ACTG", "p2": "TTTT"}), encoding="utf-8")
    result = load_upl_probes(str(f))
    assert result == {"p1": "ACTG", "p2": "TTTT"}


def test_load_upl_probes_file_not_found():
    with pytest.raises(FileNotFoundError):
        load_upl_probes("/nonexistent/path/probes.tsv")


def test_load_upl_probes_pkl_raises():
    with tempfile.NamedTemporaryFile(suffix=".pkl", delete=False) as tf:
        tf.write(b"dummy")
        path = tf.name
    with pytest.raises(ValueError, match=r"\.pkl"):
        load_upl_probes(path)


def test_load_upl_probes_unsupported_extension(tmp_path):
    f = tmp_path / "probes.csv"
    f.write_text("p1,ACTG\n", encoding="utf-8")
    with pytest.raises(ValueError):
        load_upl_probes(str(f))


def test_load_upl_probes_empty_file_raises(tmp_path):
    f = tmp_path / "probes.tsv"
    f.write_text("# only comments\n", encoding="utf-8")
    with pytest.raises(ValueError, match="no valid probe"):
        load_upl_probes(str(f))


# ---------------------------------------------------------------------------
# find_upl_matches
# ---------------------------------------------------------------------------

_PROBES = {"p1": "ACTG", "p2": "CCCC"}


def test_find_upl_matches_forward():
    hits = find_upl_matches("NNACTGNN", _PROBES)
    match = next((h for h in hits if h.probe_id == "p1" and h.strand == "+"), None)
    assert match is not None
    assert match.start == 2
    assert match.end == 5


def test_find_upl_matches_reverse_complement():
    # ACTG revcomp = CAGT
    hits = find_upl_matches("NNCAGTNN", _PROBES)
    match = next((h for h in hits if h.probe_id == "p1" and h.strand == "-"), None)
    assert match is not None
    assert match.start == 2
    assert match.end == 5


def test_find_upl_matches_no_hit():
    hits = find_upl_matches("NNGGGGNN", _PROBES)
    assert not any(h.probe_id == "p1" for h in hits)


def test_find_upl_matches_sorted_by_start():
    probes = {"p1": "ACTG"}
    seq = "ACTGNNACTG"
    hits = find_upl_matches(seq, probes)
    starts = [h.start for h in hits if h.probe_id == "p1" and h.strand == "+"]
    assert starts == [0, 6]


def test_find_upl_matches_palindrome_not_duplicated():
    # ACGT is its own reverse complement
    probes = {"pal": "ACGT"}
    hits = find_upl_matches("ACGT", probes)
    # palindrome: only + strand hit should appear (no duplicate -)
    pal_hits = [h for h in hits if h.probe_id == "pal"]
    assert len(pal_hits) == 1
    assert pal_hits[0].strand == "+"


def test_find_upl_matches_empty_probes():
    hits = find_upl_matches("ACTG", {})
    assert hits == []


def test_find_upl_matches_lowercase_sequence():
    hits = find_upl_matches("nnactgnn", _PROBES)
    match = next((h for h in hits if h.probe_id == "p1" and h.strand == "+"), None)
    assert match is not None
