"""Tests for upl_finder_fast.specificity module."""
from __future__ import annotations

import pytest

from upl_finder_fast.specificity import BlastHit, filter_hits_in_silico


def _make_hit(qseq: str, sseq: str) -> BlastHit:
    """Create a BlastHit with the given query/subject sequences."""
    return BlastHit(
        sseqid="chr1",
        sstart=100,
        send=107,
        sstrand="plus",
        pident=87.5,
        length=len(qseq.replace("-", "")),
        mismatch=sum(1 for q, s in zip(qseq, sseq) if q != s and q != "-" and s != "-"),
        qstart=1,
        qend=len(qseq),
        qseq=qseq,
        sseq=sseq,
    )


# ---------------------------------------------------------------------------
# require_terminal_mismatch=False (default behavior)
# ---------------------------------------------------------------------------

def test_terminal_match_filtered_normally_when_flag_false():
    """
    Default behavior: a hit that passes mismatch filter is discarded as safe,
    even if it has a perfect terminal base match.
    """
    # 8-base primer; 2 mismatches in 3' window (last 5 bases), last base matches
    # qseq: ACGTACGT
    # sseq: ACGTXCGT  (mismatch at pos 4, terminal G matches)
    # mismatch_total=1, mismatch_3p=1 → both >= thresholds (min=1,min3p=1) → filtered (safe)
    hit = _make_hit("ACGTACGT", "ACGTXCGT")
    result = filter_hits_in_silico(
        hits=[hit],
        primer_len=8,
        min_mismatches_total=1,
        min_mismatches_3p=1,
        three_prime_window=5,
        ignore_mismatches_total_ge=6,
        require_terminal_mismatch=False,
    )
    assert len(result) == 0  # safe → filtered out


def test_terminal_match_retained_when_flag_set():
    """
    When require_terminal_mismatch=True, a hit with a perfect 3'-terminal base
    match is retained as a potential off-target even if mismatch counts would
    otherwise mark it safe.
    """
    # 8-base primer; qseq=ACGTACGT sseq=ACGTXCGT
    # mismatch_total=1, mismatch_3p=1 → normally safe (both >= threshold of 1)
    # But terminal base (T vs T) matches → retained when flag is True
    hit = _make_hit("ACGTACGT", "ACGTXCGT")
    result = filter_hits_in_silico(
        hits=[hit],
        primer_len=8,
        min_mismatches_total=1,
        min_mismatches_3p=1,
        three_prime_window=5,
        ignore_mismatches_total_ge=6,
        require_terminal_mismatch=True,
    )
    assert len(result) == 1  # retained as off-target


def test_terminal_mismatch_always_filtered_when_flag_set():
    """
    When require_terminal_mismatch=True, a hit with a terminal base mismatch
    follows normal filter logic (if mismatch counts exceed thresholds, it's safe).
    """
    # qseq=ACGTACGA sseq=ACGTXCGX
    # terminal: A vs X → mismatch at terminal → terminal_match=False
    # mismatch_total=2, mismatch_3p=2 → both >= threshold of 2 → normally safe
    # With flag=True but terminal_match=False → still safe (filtered)
    hit = _make_hit("ACGTACGA", "ACGTXCGX")
    result = filter_hits_in_silico(
        hits=[hit],
        primer_len=8,
        min_mismatches_total=2,
        min_mismatches_3p=2,
        three_prime_window=5,
        ignore_mismatches_total_ge=6,
        require_terminal_mismatch=True,
    )
    assert len(result) == 0  # safe → filtered out


def test_excessive_mismatches_always_filtered():
    """Hits with excessive mismatches are always filtered, even with require_terminal_mismatch=True."""
    # qseq vs sseq with 6 mismatches (ignore_mismatches_total_ge=6 → filtered)
    hit = _make_hit("ACGTACGT", "X" * 8)
    result = filter_hits_in_silico(
        hits=[hit],
        primer_len=8,
        min_mismatches_total=1,
        min_mismatches_3p=1,
        three_prime_window=5,
        ignore_mismatches_total_ge=6,
        require_terminal_mismatch=True,
    )
    assert len(result) == 0  # too many mismatches → always filtered


def test_filter_hits_default_no_require_terminal_mismatch():
    """Default require_terminal_mismatch=False preserves existing behavior."""
    # Perfect match hit: mismatch_total=0, mismatch_3p=0 → retained as off-target
    hit = _make_hit("ACGTACGT", "ACGTACGT")
    result = filter_hits_in_silico(
        hits=[hit],
        primer_len=8,
        min_mismatches_total=2,
        min_mismatches_3p=2,
        three_prime_window=5,
        ignore_mismatches_total_ge=6,
    )
    assert len(result) == 1  # perfect match retained as off-target
