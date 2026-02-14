from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class PrimerSpecificity:
    exact_hit_count: int
    example_subjects: list[str]


@dataclass(frozen=True)
class BlastHit:
    sseqid: str
    sstart: int
    send: int
    sstrand: str
    pident: float
    length: int
    mismatch: int
    qstart: int | None = None
    qend: int | None = None
    qseq: str | None = None
    sseq: str | None = None
    mismatch_total: int | None = None
    mismatch_3p: int | None = None

    def three_prime_pos(self) -> int:
        # For plus strand primers, 3' end is at the higher coordinate.
        # For minus strand primers, 3' end is at the lower coordinate.
        # In BLAST output: plus has sstart < send (returns send via max()),
        # minus has sstart > send (returns send via min()).
        return max(self.sstart, self.send) if self.sstrand == "plus" else min(self.sstart, self.send)


def blastn_available(blastn_path: str = "blastn") -> bool:
    return shutil.which(blastn_path) is not None


def blastn_exact_hit_count(*, primer_seq: str, db: str, blastn_path: str = "blastn", max_targets: int = 500) -> PrimerSpecificity:
    """
    BLASTN (blastn-short) で primer の完全一致ヒット数を概算します。
    - outfmt6で返し、pident=100 & length=len(primer) & mismatch=0 のみをカウントします。
    - DBは makeblastdb 済みのものを想定します。
    """
    hits = blastn_hits(primer_seq=primer_seq, db=db, blastn_path=blastn_path, max_targets=max_targets)
    primer = primer_seq.upper().replace(" ", "").replace("\n", "")
    exact_hits = [h for h in hits if _is_exact_hit(h, primer_len=len(primer))]
    subjects: list[str] = []
    for h in exact_hits:
        if len(subjects) < 5 and h.sseqid not in subjects:
            subjects.append(h.sseqid)
    return PrimerSpecificity(exact_hit_count=len(exact_hits), example_subjects=subjects)


def blastn_hits(*, primer_seq: str, db: str, blastn_path: str = "blastn", max_targets: int = 500) -> list[BlastHit]:
    """
    BLASTN (blastn-short) で primer ヒットを取得します。
    outfmt6 + sstrand を使い、方向と座標を保持します。
    """
    if not blastn_available(blastn_path):
        raise RuntimeError(f"blastn が見つかりません: {blastn_path}")

    if not _blast_db_exists(db):
        raise RuntimeError(f"BLAST DB が見つかりません: {db}")

    primer = primer_seq.upper().replace(" ", "").replace("\n", "")
    if not primer:
        return []

    with tempfile.TemporaryDirectory() as td:
        q = Path(td) / "q.fa"
        q.write_text(f">q\n{primer}\n", encoding="utf-8")

        cmd = [
            blastn_path,
            "-task",
            "blastn-short",
            "-db",
            db,
            "-query",
            str(q),
            "-max_target_seqs",
            str(int(max_targets)),
            "-dust",
            "no",
            "-soft_masking",
            "false",
            "-outfmt",
            "6 sseqid sstart send sstrand pident length mismatch",
        ]
        p = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if p.returncode != 0:
            raise RuntimeError(f"blastn failed (code {p.returncode}): {p.stderr[:200]}")

        hits: list[BlastHit] = []
        for line in p.stdout.splitlines():
            parts = line.split("\t")
            if len(parts) < 7:
                continue
            sseqid = parts[0]
            try:
                sstart = int(parts[1])
                send = int(parts[2])
                sstrand = parts[3]
                pident = float(parts[4])
                length = int(parts[5])
                mismatch = int(parts[6])
            except Exception:
                continue
            hits.append(
                BlastHit(
                    sseqid=sseqid,
                    sstart=sstart,
                    send=send,
                    sstrand=sstrand,
                    pident=pident,
                    length=length,
                    mismatch=mismatch,
                    qstart=None,
                    qend=None,
                    qseq=None,
                    sseq=None,
                    mismatch_total=None,
                    mismatch_3p=None,
                )
            )
        return hits


def blastn_hits_detailed(
    *, primer_seq: str, db: str, blastn_path: str = "blastn", max_targets: int = 500
) -> list[BlastHit]:
    """
    BLASTN (blastn-short) で primer ヒットを取得します（ミスマッチ位置推定用）。
    qseq/sseq を含めて返します。
    """
    if not blastn_available(blastn_path):
        raise RuntimeError(f"blastn が見つかりません: {blastn_path}")

    if not _blast_db_exists(db):
        raise RuntimeError(f"BLAST DB が見つかりません: {db}")

    primer = primer_seq.upper().replace(" ", "").replace("\n", "")
    if not primer:
        return []

    with tempfile.TemporaryDirectory() as td:
        q = Path(td) / "q.fa"
        q.write_text(f">q\n{primer}\n", encoding="utf-8")

        cmd = [
            blastn_path,
            "-task",
            "blastn-short",
            "-db",
            db,
            "-query",
            str(q),
            "-max_target_seqs",
            str(int(max_targets)),
            "-dust",
            "no",
            "-soft_masking",
            "false",
            "-outfmt",
            "6 sseqid sstart send sstrand pident length mismatch qstart qend qseq sseq",
        ]
        p = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if p.returncode != 0:
            raise RuntimeError(f"blastn failed (code {p.returncode}): {p.stderr[:200]}")

        hits: list[BlastHit] = []
        for line in p.stdout.splitlines():
            parts = line.split("\t")
            if len(parts) < 11:
                continue
            sseqid = parts[0]
            try:
                sstart = int(parts[1])
                send = int(parts[2])
                sstrand = parts[3]
                pident = float(parts[4])
                length = int(parts[5])
                mismatch = int(parts[6])
                qstart = int(parts[7])
                qend = int(parts[8])
                qseq = parts[9]
                sseq = parts[10]
            except Exception:
                continue
            hits.append(
                BlastHit(
                    sseqid=sseqid,
                    sstart=sstart,
                    send=send,
                    sstrand=sstrand,
                    pident=pident,
                    length=length,
                    mismatch=mismatch,
                    qstart=qstart,
                    qend=qend,
                    qseq=qseq,
                    sseq=sseq,
                    mismatch_total=None,
                    mismatch_3p=None,
                )
            )
        return hits


def count_amplicon_products(
    *,
    left_hits: Iterable[BlastHit],
    right_hits: Iterable[BlastHit],
    min_size: int,
    max_size: int,
) -> int:
    """
    left は plus strand、right は minus strand の完全一致ヒットを前提に、
    産物サイズ内のPCR産物候補数を数えます。
    """
    left_by_contig: dict[str, list[BlastHit]] = {}
    right_by_contig: dict[str, list[BlastHit]] = {}
    for h in left_hits:
        left_by_contig.setdefault(h.sseqid, []).append(h)
    for h in right_hits:
        right_by_contig.setdefault(h.sseqid, []).append(h)

    products: set[tuple[str, int, int]] = set()
    for contig, lhits in left_by_contig.items():
        rhits = right_by_contig.get(contig)
        if not rhits:
            continue
        for l in lhits:
            if l.sstrand != "plus":
                continue
            l3 = l.three_prime_pos()
            for r in rhits:
                if r.sstrand != "minus":
                    continue
                r3 = r.three_prime_pos()
                if r3 <= l3:
                    continue
                size = r3 - l3 + 1
                if min_size <= size <= max_size:
                    products.add((contig, l3, r3))
    return len(products)


def _is_exact_hit(h: BlastHit, *, primer_len: int) -> bool:
    return h.mismatch == 0 and h.length == primer_len and h.pident == 100.0


def filter_hits_in_silico(
    *,
    hits: list[BlastHit],
    primer_len: int,
    min_mismatches_total: int,
    min_mismatches_3p: int,
    three_prime_window: int,
    ignore_mismatches_total_ge: int,
) -> list[BlastHit]:
    """
    Filter BLAST hits to identify primers that could cause off-target amplification.
    
    Returns hits that are potential off-targets (not filtered out as "safe").
    
    A hit is considered "safe" (filtered out) if BOTH conditions are met:
    - Total mismatches >= min_mismatches_total AND
    - 3' window mismatches >= min_mismatches_3p
    
    This reflects PCR biology: off-target amplification requires BOTH
    overall sequence similarity AND good 3' end matching.
    
    Hits with excessive mismatches (>= ignore_mismatches_total_ge) are also
    filtered as they're unlikely to amplify under any conditions.
    """
    out: list[BlastHit] = []
    for h in hits:
        stats = _mismatch_stats(
            qseq=h.qseq,
            sseq=h.sseq,
            primer_len=primer_len,
            three_prime_window=three_prime_window,
        )
        if stats is None:
            continue
        mismatch_total, mismatch_3p, terminal_match = stats
        if mismatch_total >= ignore_mismatches_total_ge:
            continue
        if (mismatch_total >= min_mismatches_total) and (mismatch_3p >= min_mismatches_3p):
            continue
        out.append(
            BlastHit(
                sseqid=h.sseqid,
                sstart=h.sstart,
                send=h.send,
                sstrand=h.sstrand,
                pident=h.pident,
                length=h.length,
                mismatch=h.mismatch,
                qstart=h.qstart,
                qend=h.qend,
                qseq=h.qseq,
                sseq=h.sseq,
                mismatch_total=mismatch_total,
                mismatch_3p=mismatch_3p,
            )
        )
    return out


def _mismatch_stats(
    *, qseq: str | None, sseq: str | None, primer_len: int, three_prime_window: int
) -> tuple[int, int, bool] | None:
    """
    Calculate mismatch statistics for primer-template alignment.
    
    Returns (total_mismatches, mismatches_in_3p_window, terminal_base_matches) or None.
    
    The 3' end of primers is critical for PCR specificity because:
    - DNA polymerase extends from 3' end
    - Mismatches at 3' end strongly inhibit primer extension
    - Terminal mismatch is especially important
    
    Filters out alignments with gaps (indels) as they indicate poor binding.
    """
    if not qseq or not sseq:
        return None
    pairs: list[tuple[str, str]] = []
    for q, s in zip(qseq, sseq):
        if q == "-" or s == "-":
            continue
        pairs.append((q, s))
    if len(pairs) != primer_len:
        # Alignment doesn't cover full primer length (has gaps) - reject
        return None
    mismatch_total = sum(1 for q, s in pairs if q != s)
    window = pairs[-max(1, int(three_prime_window)) :]
    mismatch_3p = sum(1 for q, s in window if q != s)
    terminal_match = window[-1][0] == window[-1][1] if window else False
    return (mismatch_total, mismatch_3p, terminal_match)


def _blast_db_exists(prefix: str) -> bool:
    """
    BLAST DB のprefixは単体ファイルではなく、.nin/.nsq/.nhr などのセット。
    そのため prefix 自体の存在ではなく、想定拡張子の存在で判定する。
    """
    p = Path(prefix)
    if p.exists():
        return True
    candidates = [f"{prefix}.nin", f"{prefix}.nsq", f"{prefix}.nhr", f"{prefix}.nal"]
    return any(Path(c).exists() for c in candidates)
