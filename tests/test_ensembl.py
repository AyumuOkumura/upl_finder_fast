"""Tests for upl_finder_fast.ensembl module cache release tracking."""
from __future__ import annotations

import pytest

from upl_finder_fast.ensembl import EnsemblClient


def test_cache_release_mismatch_warns(tmp_path):
    """A stale cache (different release) should emit a UserWarning."""
    release_file = tmp_path / "ensembl_release.txt"
    release_file.write_text("109", encoding="utf-8")

    client = EnsemblClient.__new__(EnsemblClient)
    client.base_url = "https://rest.ensembl.org"
    client.cache_dir = tmp_path
    client._ensembl_release = "112"

    with pytest.warns(UserWarning, match="release mismatch"):
        client._check_cache_release()


def test_cache_release_no_warn_when_same(tmp_path):
    """No warning when cached release matches current release."""
    release_file = tmp_path / "ensembl_release.txt"
    release_file.write_text("112", encoding="utf-8")

    client = EnsemblClient.__new__(EnsemblClient)
    client.base_url = "https://rest.ensembl.org"
    client.cache_dir = tmp_path
    client._ensembl_release = "112"

    # Should not warn
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        client._check_cache_release()  # must not raise


def test_cache_release_no_warn_when_unknown(tmp_path):
    """No warning when current release is unknown (API unreachable)."""
    release_file = tmp_path / "ensembl_release.txt"
    release_file.write_text("109", encoding="utf-8")

    client = EnsemblClient.__new__(EnsemblClient)
    client.base_url = "https://rest.ensembl.org"
    client.cache_dir = tmp_path
    client._ensembl_release = "unknown"

    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        client._check_cache_release()  # must not raise even though mismatch


def test_cache_release_written_when_known(tmp_path):
    """Release number should be written to ensembl_release.txt when known."""
    client = EnsemblClient.__new__(EnsemblClient)
    client.base_url = "https://rest.ensembl.org"
    client.cache_dir = tmp_path
    client._ensembl_release = "112"

    client._check_cache_release()

    release_file = tmp_path / "ensembl_release.txt"
    assert release_file.exists()
    assert release_file.read_text(encoding="utf-8").strip() == "112"


def test_cache_release_not_written_when_unknown(tmp_path):
    """Release file should not be written when release is unknown."""
    client = EnsemblClient.__new__(EnsemblClient)
    client.base_url = "https://rest.ensembl.org"
    client.cache_dir = tmp_path
    client._ensembl_release = "unknown"

    client._check_cache_release()

    release_file = tmp_path / "ensembl_release.txt"
    assert not release_file.exists()
