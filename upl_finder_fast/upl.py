from __future__ import annotations

import json
import pickle
from dataclasses import dataclass
from typing import Any
from pathlib import Path

try:  # optional Rust acceleration
    import upl_finder_rust as _upl_rust  # type: ignore
except Exception:  # pragma: no cover
    _upl_rust = None

_RUST_FINDER: Any | None = None
_RUST_PROBES_ID: int | None = None


def load_upl_probes(path: str) -> dict[str, str]:
    p = Path(path)
    if p.suffix == ".pkl" and p.exists():
        with p.open("rb") as f:
            obj = pickle.load(f)
    else:
        cache = p.with_suffix(".pkl")
        if cache.exists() and (not p.exists() or cache.stat().st_mtime >= p.stat().st_mtime):
            try:
                with cache.open("rb") as f:
                    obj = pickle.load(f)
            except Exception:
                obj = None
        else:
            obj = None

    if obj is None:
        obj = json.loads(p.read_text(encoding="utf-8"))
    if not isinstance(obj, dict):
        raise ValueError("UPL JSON must be a dict: {probe_id: sequence}")
    probes: dict[str, str] = {}
    for k, v in obj.items():
        if not isinstance(k, str) or not isinstance(v, str):
            continue
        probes[k] = v.upper().replace(" ", "").replace("\n", "")
    if p.exists() and p.suffix != ".pkl":
        try:
            with cache.open("wb") as f:
                pickle.dump(probes, f, protocol=pickle.HIGHEST_PROTOCOL)
        except Exception:
            pass
    return probes


@dataclass(frozen=True)
class ProbeMatch:
    probe_id: str
    probe_seq: str
    start: int
    end: int


def find_upl_matches(sequence: str, upl_probes: dict[str, str]) -> list[ProbeMatch]:
    if _upl_rust is not None:
        global _RUST_FINDER, _RUST_PROBES_ID
        if _RUST_FINDER is None or _RUST_PROBES_ID != id(upl_probes):
            try:
                _RUST_FINDER = _upl_rust.UplFinder(upl_probes)
                _RUST_PROBES_ID = id(upl_probes)
            except Exception:
                _RUST_FINDER = None
                _RUST_PROBES_ID = None
        if _RUST_FINDER is not None:
            try:
                hits = _RUST_FINDER.find(sequence)
                hits.sort(key=lambda x: (x.start, x.probe_id))
                return hits
            except Exception:
                pass
    seq = sequence.upper()
    candidates = list(upl_probes.items())
    # Safe prefilter: if probe fully matches, its prefix k-mer must appear in seq.
    kmer_len = 8
    if len(seq) >= 150 and len(candidates) >= 20 and len(seq) >= kmer_len:
        kmers = {seq[i : i + kmer_len] for i in range(len(seq) - kmer_len + 1)}
        filtered: list[tuple[str, str]] = []
        for probe_id, probe_seq in candidates:
            if len(probe_seq) < kmer_len:
                filtered.append((probe_id, probe_seq))
            elif probe_seq[:kmer_len] in kmers:
                filtered.append((probe_id, probe_seq))
        candidates = filtered
    hits: list[ProbeMatch] = []
    for probe_id, probe_seq in candidates:
        start = seq.find(probe_seq)
        while start != -1:
            hits.append(ProbeMatch(probe_id=probe_id, probe_seq=probe_seq, start=start, end=start + len(probe_seq) - 1))
            start = seq.find(probe_seq, start + 1)
    hits.sort(key=lambda x: (x.start, x.probe_id))
    return hits
