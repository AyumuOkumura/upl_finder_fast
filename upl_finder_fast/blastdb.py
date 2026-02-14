from __future__ import annotations

import gzip
import shutil
import subprocess
import urllib.request
from pathlib import Path
from typing import Callable


ENSEMBL_BASE = "https://ftp.ensembl.org/pub/current_fasta"

BLAST_SOURCES = [
    (
        "human",
        "cdna",
        f"{ENSEMBL_BASE}/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
    ),
    (
        "human",
        "genome",
        f"{ENSEMBL_BASE}/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
    ),
    (
        "mouse",
        "cdna",
        f"{ENSEMBL_BASE}/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz",
    ),
    (
        "mouse",
        "genome",
        f"{ENSEMBL_BASE}/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz",
    ),
]


def ensure_default_blastdbs(
    base_dir: Path, *, progress_cb: Callable[[float, str], None] | None = None
) -> None:
    db_dir = base_dir / "data" / "blastdb"
    db_dir.mkdir(parents=True, exist_ok=True)
    total = len(BLAST_SOURCES)
    for i, (species, kind, url) in enumerate(BLAST_SOURCES, start=1):
        prefix = db_dir / f"{species}_{kind}"
        if _blast_db_exists(prefix):
            continue
        _ensure_makeblastdb()
        if progress_cb is not None:
            pct = 2.0 + (8.0 * i / max(1, total))
            progress_cb(pct, f"Downloading BLAST DB: {species} {kind}...")
        _download_and_makeblastdb(url=url, prefix=prefix)


def _download_and_makeblastdb(*, url: str, prefix: Path) -> None:
    gz_path = prefix.with_suffix(".fa.gz")
    fa_path = prefix.with_suffix(".fa")
    try:
        _download(url, gz_path)
        _gunzip(gz_path, fa_path)
        _makeblastdb(fa_path, prefix)
    finally:
        if gz_path.exists():
            gz_path.unlink()
        if fa_path.exists():
            fa_path.unlink()


def _download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url) as r, open(dest, "wb") as f:
        shutil.copyfileobj(r, f)


def _gunzip(src: Path, dest: Path) -> None:
    with gzip.open(src, "rb") as f_in, open(dest, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)


def _makeblastdb(fasta: Path, prefix: Path) -> None:
    cmd = ["makeblastdb", "-in", str(fasta), "-dbtype", "nucl", "-out", str(prefix)]
    p = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if p.returncode != 0:
        raise RuntimeError(f"makeblastdb failed (code {p.returncode}): {p.stderr[:200]}")


def _ensure_makeblastdb() -> None:
    if shutil.which("makeblastdb") is None:
        raise RuntimeError("makeblastdb が見つかりません。BLAST+ をインストールしてください。")


def _blast_db_exists(prefix: Path) -> bool:
    if prefix.exists():
        return True
    candidates = [f"{prefix}.nin", f"{prefix}.nsq", f"{prefix}.nhr", f"{prefix}.nal"]
    return any(Path(c).exists() for c in candidates)
