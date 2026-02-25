from __future__ import annotations

import os
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

try:  # optional Rust acceleration
    import upl_finder_rust as _upl_rust  # type: ignore
except Exception:  # pragma: no cover
    _upl_rust = None

_RUST_FINDER: Any | None = None
_RUST_PROBES_ID: int | None = None


def rust_extension_available() -> bool:
    return _upl_rust is not None


def set_rust_mode(mode: str) -> None:
    v = mode.strip().lower()
    if v not in {"auto", "on", "off"}:
        raise ValueError(f"Unsupported rust mode: {mode}")
    if v == "on":
        os.environ["UPL_FINDER_DISABLE_RUST"] = "0"
    elif v == "off":
        os.environ["UPL_FINDER_DISABLE_RUST"] = "1"


def _rust_enabled() -> bool:
    return _upl_rust is not None and os.getenv("UPL_FINDER_DISABLE_RUST", "0") != "1"


def _get_rust_finder(upl_probes: dict[str, str]) -> Any | None:
    global _RUST_FINDER, _RUST_PROBES_ID
    if not _rust_enabled():
        return None
    if _RUST_FINDER is None or _RUST_PROBES_ID != id(upl_probes):
        try:
            _RUST_FINDER = _upl_rust.UplFinder(upl_probes)
            _RUST_PROBES_ID = id(upl_probes)
        except Exception:
            _RUST_FINDER = None
            _RUST_PROBES_ID = None
    return _RUST_FINDER


def load_upl_probes(path: str) -> dict[str, str]:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    suffix = p.suffix.lower()
    obj: object
    if suffix == ".pkl":
        raise ValueError("UPL probe data in .pkl is no longer supported; use .tsv/.txt or .json instead")
    if suffix == ".json":
        obj = json.loads(p.read_text(encoding="utf-8"))
    elif suffix in {".tsv", ".txt"}:
        obj = _load_upl_probes_tsv(p.read_text(encoding="utf-8"))
    else:
        raise ValueError("UPL probe data must be one of: .json, .tsv, .txt")
    if not isinstance(obj, dict):
        raise ValueError("UPL probe data must be a dict: {probe_id: sequence}")
    probes: dict[str, str] = {}
    for k, v in obj.items():
        if not isinstance(k, str) or not isinstance(v, str):
            continue
        probes[k] = v.upper().replace(" ", "").replace("\n", "")
    if not probes:
        raise ValueError("UPL probe file contained no valid probe sequences")
    return probes


def _load_upl_probes_tsv(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for i, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if "\t" in line:
            parts = line.split("\t")
        else:
            parts = line.split()
        if len(parts) < 2:
            raise ValueError(f"Invalid TSV at line {i}: expected <probe_id>\\t<sequence>")
        probe_id = str(parts[0]).strip()
        seq = str(parts[1]).strip()
        if not probe_id or not seq:
            raise ValueError(f"Invalid TSV at line {i}: expected <probe_id>\\t<sequence>")
        out[probe_id] = seq
    return out


@dataclass(frozen=True)
class ProbeMatch:
    probe_id: str
    probe_seq: str
    start: int
    end: int
    strand: str  # "+" when matched as-is, "-" when matched as reverse-complement


def rank_candidates_with_rust(
    *,
    template: str,
    upl_probes: dict[str, str],
    candidates: list[dict[str, Any]],
    exon_boundaries: list[int],
    min_probe_offset_bp: int,
) -> list[dict[str, Any]] | None:
    if not _rust_enabled():
        return None
    try:
        rows = _upl_rust.rank_candidates(
            template,
            upl_probes,
            candidates,
            exon_boundaries,
            int(min_probe_offset_bp),
        )
        return [dict(r) for r in rows]
    except Exception:
        return None


def find_upl_matches(sequence: str, upl_probes: dict[str, str]) -> list[ProbeMatch]:
    seq = sequence.upper()
    candidates = list(upl_probes.items())
    kmer_len = 8
    if len(seq) >= 150 and len(candidates) >= 20 and len(seq) >= kmer_len:
        kmers = {seq[i : i + kmer_len] for i in range(len(seq) - kmer_len + 1)}
        filtered: list[tuple[str, str]] = []
        for probe_id, probe_seq in candidates:
            if len(probe_seq) < kmer_len:
                filtered.append((probe_id, probe_seq))
            else:
                rc = _revcomp(probe_seq)
                if probe_seq[:kmer_len] in kmers or rc[:kmer_len] in kmers:
                    filtered.append((probe_id, probe_seq))
        candidates = filtered

    hits: list[ProbeMatch] = []
    seen: set[tuple[str, int, int, str]] = set()
    allowed_probe_ids = {pid for pid, _ in candidates}

    rust_finder = _get_rust_finder(upl_probes)
    if rust_finder is not None:
        try:
            rust_hits = rust_finder.find(seq)
            for h in rust_hits:
                pid = getattr(h, "probe_id", None)
                s = getattr(h, "start", None)
                e = getattr(h, "end", None)
                strand = str(getattr(h, "strand", "+"))
                if not isinstance(pid, str) or not isinstance(s, int) or not isinstance(e, int):
                    continue
                if pid not in allowed_probe_ids:
                    continue
                probe_seq = upl_probes.get(pid)
                if probe_seq is None:
                    continue
                key = (pid, s, e, strand)
                if key in seen:
                    continue
                seen.add(key)
                hits.append(ProbeMatch(probe_id=pid, probe_seq=probe_seq, start=s, end=e, strand=strand))
            hits.sort(key=lambda x: (x.start, x.probe_id, x.strand))
            return hits
        except Exception:
            pass

    for probe_id, probe_seq in candidates:
        probe_seq = probe_seq.upper()
        for s in _find_all(seq, probe_seq):
            e = s + len(probe_seq) - 1
            key = (probe_id, s, e, "+")
            if key in seen:
                continue
            seen.add(key)
            hits.append(ProbeMatch(probe_id=probe_id, probe_seq=probe_seq, start=s, end=e, strand="+"))

        rc = _revcomp(probe_seq)
        if rc != probe_seq:
            for s in _find_all(seq, rc):
                e = s + len(rc) - 1
                key = (probe_id, s, e, "-")
                if key in seen:
                    continue
                seen.add(key)
                hits.append(ProbeMatch(probe_id=probe_id, probe_seq=probe_seq, start=s, end=e, strand="-"))

    hits.sort(key=lambda x: (x.start, x.probe_id, x.strand))
    return hits


_RC_TABLE = str.maketrans({"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"})


def _revcomp(seq: str) -> str:
    s = seq.upper().translate(_RC_TABLE)
    return s[::-1]


def _find_all(seq: str, sub: str) -> list[int]:
    if not sub:
        return []
    out: list[int] = []
    start = seq.find(sub)
    while start != -1:
        out.append(start)
        start = seq.find(sub, start + 1)
    return out
