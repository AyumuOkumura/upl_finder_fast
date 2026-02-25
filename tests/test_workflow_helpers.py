"""Tests for upl_finder_fast.workflow module helpers."""
from __future__ import annotations

import pytest

from upl_finder_fast.workflow import (
    _has_3p_gc_run,
    _has_poly_run,
    _junction_flags,
    _normalize_seq,
    _probe_preference,
    _score_candidate,
    _MAX_SPECIFICITY_RETRIES,
)
from upl_finder_fast.design import PrimerPairCandidate


# ---------------------------------------------------------------------------
# _normalize_seq
# ---------------------------------------------------------------------------

def test_normalize_seq_uppercase():
    assert _normalize_seq("actg") == "ACTG"


def test_normalize_seq_u_to_t():
    assert _normalize_seq("AUCG") == "ATCG"


def test_normalize_seq_strips_non_acgtn():
    assert _normalize_seq("ACTG-N ") == "ACTGN"


def test_normalize_seq_empty():
    assert _normalize_seq("") == ""


# ---------------------------------------------------------------------------
# _has_3p_gc_run
# ---------------------------------------------------------------------------

def test_has_3p_gc_run_true():
    assert _has_3p_gc_run("ATGCCC") is True


def test_has_3p_gc_run_false():
    assert _has_3p_gc_run("ATGCCA") is False


def test_has_3p_gc_run_exactly_at_limit():
    # run_len=3: exactly 3 GC at end → True
    assert _has_3p_gc_run("ATGCGG", run_len=3) is True
    # 2 GC at end with run_len=3 → False
    assert _has_3p_gc_run("ATGCAG", run_len=3) is False


def test_has_3p_gc_run_empty():
    assert _has_3p_gc_run("") is False


# ---------------------------------------------------------------------------
# _has_poly_run
# ---------------------------------------------------------------------------

def test_has_poly_run_true():
    assert _has_poly_run("ACTAAAAT") is True


def test_has_poly_run_false():
    assert _has_poly_run("ACTGACTG") is False


def test_has_poly_run_empty():
    assert _has_poly_run("") is False


def test_has_poly_run_single_char():
    assert _has_poly_run("A") is False


def test_has_poly_run_run_len_4_exact():
    assert _has_poly_run("AAAA") is True
    assert _has_poly_run("AAA") is False


# ---------------------------------------------------------------------------
# _junction_flags
# ---------------------------------------------------------------------------

def test_junction_flags_no_exon_info():
    spanning, detail = _junction_flags(boundaries=[], product_start=0, left_len=5, right_len=5, product_size=20)
    assert spanning is False
    assert detail == "no_exon_info"


def test_junction_flags_left_primer_spans():
    # boundary at position 3 means last base of exon 1 is at index 3
    # left primer covers positions [10..14], product_start=10
    # boundary at 12 is inside left primer → spans
    spanning, detail = _junction_flags(boundaries=[12], product_start=10, left_len=5, right_len=5, product_size=20)
    assert spanning is True
    assert detail == "left_primer_spans_junction"


def test_junction_flags_amplicon_spans():
    # boundary at 15 is inside amplicon [10..29], but not inside left primer [10..14]
    spanning, detail = _junction_flags(boundaries=[17], product_start=10, left_len=5, right_len=5, product_size=20)
    assert spanning is True
    assert detail == "amplicon_spans_junction"


def test_junction_flags_no_junction():
    # boundary is outside the amplicon
    spanning, detail = _junction_flags(boundaries=[5], product_start=10, left_len=5, right_len=5, product_size=20)
    assert spanning is False
    assert detail == "no_junction"


def test_junction_flags_right_primer_spans():
    # amplicon [0..19], right primer [14..19], right_start_guess=14
    # boundary at 15: 14 <= 15 < 19 → right spans
    # With the fix (right checked before amplicon), this returns the more precise label.
    spanning, detail = _junction_flags(boundaries=[15], product_start=0, left_len=5, right_len=6, product_size=20)
    assert spanning is True
    assert detail == "right_primer_spans_junction (approx)"


def test_junction_flags_amplicon_spans_not_in_primers():
    # junction strictly between left primer [0..4] and right primer [14..19]
    # boundary at 10 is in the amplicon interior (not in either primer)
    spanning, detail = _junction_flags(boundaries=[10], product_start=0, left_len=5, right_len=6, product_size=20)
    assert spanning is True
    assert detail == "amplicon_spans_junction"


# ---------------------------------------------------------------------------
# _probe_preference
# ---------------------------------------------------------------------------

def test_probe_preference_symmetric_center():
    # probe exactly in center with equal distances should score higher than off-center
    score_center = _probe_preference(product_size=100, probe_start=46, probe_end=53, dist_l=40, dist_r=40)
    score_off = _probe_preference(product_size=100, probe_start=10, probe_end=17, dist_l=4, dist_r=76)
    assert score_center > score_off


# ---------------------------------------------------------------------------
# _score_candidate
# ---------------------------------------------------------------------------

def _make_candidate(**kwargs) -> PrimerPairCandidate:
    defaults = dict(
        left_seq="ATATATAT",
        right_seq="TATATATA",
        product_size=80,
        left_len=8,
        right_len=8,
        tm_left=59.5,
        tm_right=59.5,
        gc_left=50.0,
        gc_right=50.0,
        pair_penalty=1.0,
        left_start=0,
        right_start=72,
    )
    defaults.update(kwargs)
    return PrimerPairCandidate(**defaults)


def test_score_candidate_junction_adds_bonus():
    cand = _make_candidate()
    score_no_jx = _score_candidate(cand, junction_spanning=False, dist_l=20, dist_r=20, probe_mid=40.0, product_size=80)
    score_jx = _score_candidate(cand, junction_spanning=True, dist_l=20, dist_r=20, probe_mid=40.0, product_size=80)
    assert score_jx > score_no_jx


def test_score_candidate_no_penalty_field():
    cand = _make_candidate(pair_penalty=None)
    score = _score_candidate(cand, junction_spanning=False, dist_l=10, dist_r=10, probe_mid=40.0, product_size=80)
    # No crash, score should still be a number
    assert isinstance(score, float)


# ---------------------------------------------------------------------------
# _MAX_SPECIFICITY_RETRIES constant
# ---------------------------------------------------------------------------

def test_max_specificity_retries_positive():
    assert isinstance(_MAX_SPECIFICITY_RETRIES, int)
    assert _MAX_SPECIFICITY_RETRIES > 0
