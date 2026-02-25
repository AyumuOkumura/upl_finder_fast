from __future__ import annotations

import hashlib
import json
import os
import time
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any

import requests


class Species(str, Enum):
    HUMAN = "homo_sapiens"
    MOUSE = "mus_musculus"


@dataclass(frozen=True)
class TranscriptRecord:
    transcript_id: str
    display_name: str | None
    biotype: str | None
    is_canonical: bool | None
    length: int | None


@dataclass(frozen=True)
class TranscriptDetails:
    transcript_id: str
    transcript_display_name: str | None
    gene_id: str | None
    gene_symbol: str | None
    strand: int | None
    exons: list[dict[str, Any]]
    cdna_sequence: str


class EnsemblClient:
    def __init__(self, base_url: str = "https://rest.ensembl.org", cache_dir: str = "data/cache/ensembl") -> None:
        self.base_url = base_url.rstrip("/")
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def resolve_refseq_mrna_to_ensembl_transcript(self, species: Species, refseq_id: str) -> tuple[str, list[str]]:
        """
        Resolve RefSeq mRNA transcript ID (e.g., NM_000546) to Ensembl transcript ID.

        Uses Ensembl REST xrefs endpoint:
          /xrefs/symbol/{species}/{refseq_id}?external_db=RefSeq_mRNA
        """
        rid = refseq_id.strip()
        if not rid:
            raise ValueError("Empty RefSeq ID")
        payload = self._get_json(
            f"/xrefs/symbol/{species.value}/{rid}",
            params={"external_db": "RefSeq_mRNA", "content-type": "application/json"},
        )
        tx_ids: list[str] = []
        seen: set[str] = set()
        for item in payload or []:
            if not isinstance(item, dict):
                continue
            if item.get("type") != "transcript":
                continue
            tid = item.get("id")
            if not isinstance(tid, str) or not tid:
                continue
            if tid in seen:
                continue
            seen.add(tid)
            tx_ids.append(tid)
        if not tx_ids:
            raise RuntimeError(f"RefSeq ID could not be resolved to Ensembl transcript: {rid} (species={species.value})")
        return (tx_ids[0], tx_ids)

    def lookup_gene_transcripts(self, species: Species, gene_symbol: str) -> list[TranscriptRecord]:
        payload = self._get_json(f"/lookup/symbol/{species.value}/{gene_symbol}", params={"expand": "1"})
        transcripts = payload.get("Transcript") or []

        out: list[TranscriptRecord] = []
        for t in transcripts:
            tid = t.get("id")
            if not isinstance(tid, str) or not tid:
                continue
            is_canonical = _to_bool_or_none(t.get("is_canonical"))
            out.append(
                TranscriptRecord(
                    transcript_id=tid,
                    display_name=t.get("display_name"),
                    biotype=t.get("biotype"),
                    is_canonical=is_canonical,
                    length=t.get("length"),
                )
            )
        out.sort(
            key=lambda x: (
                (x.biotype or "") != "protein_coding",
                x.is_canonical is not True,
                -(x.length or 0),
                x.transcript_id,
            )
        )
        return out

    def get_transcript_details(self, transcript_id: str) -> TranscriptDetails:
        meta = self._get_json(f"/lookup/id/{transcript_id}", params={"expand": "1"})
        seq = self._get_text(f"/sequence/id/{transcript_id}", params={"type": "cdna"}).strip()
        exons = meta.get("Exon") or []
        gene_id = meta.get("Parent") if isinstance(meta.get("Parent"), str) else meta.get("Parent", {}).get("id")
        gene_symbol = None
        if isinstance(gene_id, str) and gene_id:
            gene_meta = self._get_json(f"/lookup/id/{gene_id}")
            gene_symbol = gene_meta.get("display_name")
        return TranscriptDetails(
            transcript_id=transcript_id,
            transcript_display_name=meta.get("display_name"),
            gene_id=gene_id,
            gene_symbol=gene_symbol,
            strand=meta.get("strand"),
            exons=list(exons),
            cdna_sequence=seq,
        )

    def get_transcript_gene_symbol(self, transcript_id: str) -> str | None:
        meta = self._get_json(f"/lookup/id/{transcript_id}")
        gene_id = meta.get("Parent") if isinstance(meta.get("Parent"), str) else meta.get("Parent", {}).get("id")
        if isinstance(gene_id, str) and gene_id:
            gene_meta = self._get_json(f"/lookup/id/{gene_id}")
            return gene_meta.get("display_name")
        return None

    def _cache_key(self, path: str, params: dict[str, str] | None) -> str:
        raw = json.dumps({"path": path, "params": params or {}}, sort_keys=True).encode("utf-8")
        return hashlib.sha256(raw).hexdigest()

    def _get_json(self, path: str, params: dict[str, str] | None = None) -> dict[str, Any]:
        key = self._cache_key(path, params)
        fp = self.cache_dir / f"{key}.json"
        if fp.exists():
            return json.loads(fp.read_text(encoding="utf-8"))

        url = f"{self.base_url}{path}"
        headers = {"Content-Type": "application/json"}
        r = requests.get(url, params=params, headers=headers, timeout=30)
        if r.status_code != 200:
            raise RuntimeError(f"Ensembl error {r.status_code} for {path}: {r.text[:200]}")
        obj = r.json()
        fp.write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding="utf-8")
        time.sleep(0.05)
        return obj

    def _get_text(self, path: str, params: dict[str, str] | None = None) -> str:
        key = self._cache_key(path, params)
        fp = self.cache_dir / f"{key}.txt"
        if fp.exists():
            return fp.read_text(encoding="utf-8")

        url = f"{self.base_url}{path}"
        headers = {"Content-Type": "text/plain"}
        r = requests.get(url, params=params, headers=headers, timeout=30)
        if r.status_code != 200:
            raise RuntimeError(f"Ensembl error {r.status_code} for {path}: {r.text[:200]}")
        text = r.text
        fp.write_text(text, encoding="utf-8")
        time.sleep(0.05)
        return text


def _to_bool_or_none(v: object) -> bool | None:
    if isinstance(v, bool):
        return v
    if isinstance(v, int):
        return v != 0
    return None
