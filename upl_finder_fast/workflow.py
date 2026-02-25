from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
from datetime import datetime
from typing import Any, Callable
from pathlib import Path
import os
from concurrent.futures import ThreadPoolExecutor
from time import perf_counter

from upl_finder_fast.blastdb import ensure_default_blastdbs
from upl_finder_fast.design import PrimerPairCandidate, design_primers_primer3
from upl_finder_fast.ensembl import EnsemblClient, Species, TranscriptDetails
from upl_finder_fast.specificity import (
    blastn_available,
    blastn_hits,
    blastn_hits_detailed,
    filter_hits_in_silico,
    BlastHit,
    count_amplicon_products,
)
from upl_finder_fast.upl import find_upl_matches, rank_candidates_with_rust


@dataclass(frozen=True)
class DesignInputs:
    species: Species
    input_type: str
    raw_input: str
    primer_tm_opt: float
    primer_tm_min: float
    primer_tm_max: float
    primer_tm_diff_max: float
    product_size_min: int
    product_size_max: int
    min_probe_offset_bp: int
    max_pairs: int
    selected_transcript_id: str | None = None
    specificity_mode: str = "none"
    blastn_path: str = "blastn"
    transcriptome_blast_db: str | None = None
    genome_blast_db: str | None = None
    specificity_top_n: int = 20
    max_offtarget_amplicons: int = 1
    min_mismatches_total: int = 2
    min_mismatches_3p: int = 2
    three_prime_window: int = 5
    ignore_mismatches_total_ge: int = 6
    max_target_amplicon_size: int = 1000
    blastn_parallel_jobs: int = 0


@dataclass(frozen=True)
class RankedPrimerPair:
    left_seq: str
    right_seq: str
    product_size: int
    tm_left: float | None
    tm_right: float | None
    gc_left: float | None
    gc_right: float | None
    pair_penalty: float | None

    upl_probe_id: str
    upl_probe_seq: str
    upl_probe_strand: str
    probe_start_in_amplicon: int
    probe_end_in_amplicon: int
    dist_left_3p_to_probe: int
    dist_probe_to_right_3p: int
    left_start: int
    left_end: int
    right_start: int
    right_end: int
    left_len: int
    right_len: int
    amplicon_start: int
    amplicon_end: int

    junction_spanning: bool
    junction_spanning_detail: str

    tx_left_exact_hits: int | None
    tx_right_exact_hits: int | None
    genome_left_exact_hits: int | None
    genome_right_exact_hits: int | None
    tx_product_offtargets: int | None
    genome_product_offtargets: int | None
    score: float

    tx_amplicons: list[dict[str, Any]] = field(default_factory=list)
    genome_amplicons: list[dict[str, Any]] = field(default_factory=list)
    specificity_passed: bool | None = None

    def to_row(self) -> dict[str, Any]:
        return {
            "score": self.score,
            "left_seq": self.left_seq,
            "right_seq": self.right_seq,
            "product_size": self.product_size,
            "tm_left": self.tm_left,
            "tm_right": self.tm_right,
            "gc_left": self.gc_left,
            "gc_right": self.gc_right,
            "pair_penalty": self.pair_penalty,
            "upl_probe_id": self.upl_probe_id,
            "upl_probe_seq": self.upl_probe_seq,
            "upl_probe_strand": self.upl_probe_strand,
            "probe_start": self.probe_start_in_amplicon,
            "probe_end": self.probe_end_in_amplicon,
            "dist_left3p_probe": self.dist_left_3p_to_probe,
            "dist_probe_right3p": self.dist_probe_to_right_3p,
            "left_start": self.left_start,
            "left_end": self.left_end,
            "right_start": self.right_start,
            "right_end": self.right_end,
            "left_len": self.left_len,
            "right_len": self.right_len,
            "amplicon_start": self.amplicon_start,
            "amplicon_end": self.amplicon_end,
            "junction_spanning": self.junction_spanning,
            "junction_detail": self.junction_spanning_detail,
            "tx_left_exact_hits": self.tx_left_exact_hits,
            "tx_right_exact_hits": self.tx_right_exact_hits,
            "genome_left_exact_hits": self.genome_left_exact_hits,
            "genome_right_exact_hits": self.genome_right_exact_hits,
            "tx_product_offtargets": self.tx_product_offtargets,
            "genome_product_offtargets": self.genome_product_offtargets,
            "tx_amplicons": self.tx_amplicons,
            "genome_amplicons": self.genome_amplicons,
            "specificity_passed": self.specificity_passed,
        }


@dataclass
class DesignResult:
    species: Species
    transcript_info: dict[str, Any] | None
    pairs: list[RankedPrimerPair]
    warnings: list[str]
    specificity_params: dict[str, Any] | None = None
    design_started_at: str | None = None
    design_elapsed_sec: float | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "species": self.species.value,
            "transcript_info": self.transcript_info,
            "pairs": [asdict(p) for p in self.pairs],
            "warnings": self.warnings,
            "specificity_params": self.specificity_params,
            "design_started_at": self.design_started_at,
            "design_elapsed_sec": self.design_elapsed_sec,
        }


def run_design_workflow(
    *,
    inputs: DesignInputs,
    ensembl: EnsemblClient,
    upl_probes: dict[str, str],
    progress_cb: Callable[[float, str], None] | None = None,
) -> DesignResult:
    warnings: list[str] = []
    if progress_cb is not None:
        progress_cb(0.0, "Preparing input...")

    transcript_details: TranscriptDetails | None = None
    template: str
    transcript_info: dict[str, Any] | None = None
    target_transcript_ids: set[str] | None = None

    if inputs.input_type == "Paste cDNA sequence":
        template = _normalize_seq(inputs.raw_input)
        transcript_info = {"transcript_id": None, "gene_symbol": None, "note": "pasted_sequence"}
    elif inputs.input_type == "Ensembl Transcript ID":
        transcript_details = ensembl.get_transcript_details(inputs.raw_input)
        template = _normalize_seq(transcript_details.cdna_sequence)
        transcript_info = _transcript_info_dict(transcript_details)
    elif inputs.input_type == "RefSeq Transcript ID":
        refseq_id = inputs.raw_input
        chosen_tid, candidates = ensembl.resolve_refseq_mrna_to_ensembl_transcript(inputs.species, refseq_id)
        if len(candidates) > 1:
            warnings.append(
                f"RefSeq IDから複数のEnsembl transcriptが解決されました。先頭を採用します。/ "
                f"Multiple Ensembl transcripts resolved from RefSeq {refseq_id}; using {chosen_tid}; "
                f"candidates={candidates}"
            )
        transcript_details = ensembl.get_transcript_details(chosen_tid)
        template = _normalize_seq(transcript_details.cdna_sequence)
        transcript_info = _transcript_info_dict(transcript_details)
        transcript_info["refseq_transcript_id_input"] = refseq_id
        transcript_info["selected_transcript_id"] = chosen_tid
    else:
        gene_symbol = inputs.raw_input
        records = ensembl.lookup_gene_transcripts(inputs.species, gene_symbol)
        if not records:
            raise RuntimeError(f"No transcripts found for gene symbol: {gene_symbol}")
        target_transcript_ids = {r.transcript_id for r in records}

        chosen = records[0]
        if inputs.selected_transcript_id:
            m = [r for r in records if r.transcript_id == inputs.selected_transcript_id]
            if m:
                chosen = m[0]
            else:
                warnings.append(
                    "指定transcriptが候補に見つからないため自動選択にフォールバックしました。"
                    f" / Fallback to auto-selected transcript: {inputs.selected_transcript_id}"
                )

        transcript_details = ensembl.get_transcript_details(chosen.transcript_id)
        template = _normalize_seq(transcript_details.cdna_sequence)
        transcript_info = _transcript_info_dict(transcript_details)
        transcript_info["gene_symbol_input"] = gene_symbol
        transcript_info["selected_transcript_id"] = chosen.transcript_id
        transcript_info["selected_transcript_biotype"] = chosen.biotype

        if len(records) > 1:
            warnings.append(
                "複数転写産物があります。/ Multiple transcripts found. "
                f"Selected: {chosen.transcript_id} (biotype={chosen.biotype}, canonical={chosen.is_canonical})"
            )

    target_gene_symbol: str | None = transcript_details.gene_symbol if transcript_details is not None else None
    if target_transcript_ids is None and target_gene_symbol:
        try:
            target_transcript_ids = {r.transcript_id for r in ensembl.lookup_gene_transcripts(inputs.species, target_gene_symbol)}
        except Exception:
            target_transcript_ids = None

    if progress_cb is not None:
        progress_cb(15.0, "Designing primer candidates with Primer3...")

    candidates = design_primers_primer3(
        template,
        product_size_min=inputs.product_size_min,
        product_size_max=inputs.product_size_max,
        primer_tm_opt=inputs.primer_tm_opt,
        primer_tm_min=inputs.primer_tm_min,
        primer_tm_max=inputs.primer_tm_max,
        primer_tm_diff_max=inputs.primer_tm_diff_max,
        max_pairs=inputs.max_pairs,
    )

    boundaries = _cdna_exon_boundaries(transcript_details) if transcript_details is not None else []

    ranked: list[RankedPrimerPair] = []
    rust_rows = rank_candidates_with_rust(
        template=template,
        upl_probes=upl_probes,
        candidates=[asdict(c) for c in candidates],
        exon_boundaries=boundaries,
        min_probe_offset_bp=inputs.min_probe_offset_bp,
    )
    if rust_rows is not None:
        ranked = [_ranked_pair_from_rust_row(row) for row in rust_rows]
        if progress_cb is not None:
            progress_cb(60.0, "Matching UPL probes...")
    else:
        total_cands = len(candidates)
        for i, cand in enumerate(candidates):
            product_start, product_end = cand.product_range()
            amplicon = template[product_start:product_end]
            if len(amplicon) != cand.product_size:
                continue

            if _has_3p_gc_run(cand.left_seq) or _has_3p_gc_run(cand.right_seq):
                continue
            if _has_poly_run(cand.left_seq) or _has_poly_run(cand.right_seq):
                continue

            hits = find_upl_matches(amplicon, upl_probes)
            if not hits:
                continue

            junction_spanning, detail = _junction_flags(
                boundaries=boundaries,
                product_start=product_start,
                left_len=cand.left_len,
                right_len=cand.right_len,
                product_size=cand.product_size,
            )

            best_hit: tuple[float, object, int, int] | None = None
            for hit in hits:
                # Amplicon-relative coordinates (0-based, inclusive)
                left_3p = cand.left_len - 1
                right_3p = cand.product_size - 1
                right_primer_start = cand.product_size - cand.right_len

                # Always forbid any overlap between primers and probe, even when min_probe_offset_bp=0.
                # left primer occupies [0, left_3p], right primer occupies [right_primer_start, right_3p].
                if int(hit.start) <= left_3p:
                    continue
                if int(hit.end) >= right_primer_start:
                    continue

                # Distance in bp between primer 3' end and nearest probe base (0 when adjacent)
                dist_l = int(hit.start) - left_3p - 1
                dist_r = right_3p - int(hit.end) - 1
                if dist_l < inputs.min_probe_offset_bp or dist_r < inputs.min_probe_offset_bp:
                    continue
                pref = _probe_preference(
                    product_size=cand.product_size,
                    probe_start=int(hit.start),
                    probe_end=int(hit.end),
                    dist_l=dist_l,
                    dist_r=dist_r,
                )
                if best_hit is None or pref > best_hit[0]:
                    best_hit = (pref, hit, dist_l, dist_r)

            if best_hit is None:
                continue

            hit = best_hit[1]
            dist_l = best_hit[2]
            dist_r = best_hit[3]
            probe_mid = (int(hit.start) + int(hit.end)) / 2.0
            score = _score_candidate(
                cand,
                junction_spanning=junction_spanning,
                dist_l=dist_l,
                dist_r=dist_r,
                probe_mid=probe_mid,
                product_size=cand.product_size,
            )
            left_start = cand.left_start
            left_end = cand.left_start + cand.left_len - 1
            right_end = product_start + cand.product_size - 1
            right_start = right_end - cand.right_len + 1
            ranked.append(
                RankedPrimerPair(
                    left_seq=cand.left_seq,
                    right_seq=cand.right_seq,
                    product_size=cand.product_size,
                    tm_left=cand.tm_left,
                    tm_right=cand.tm_right,
                    gc_left=cand.gc_left,
                    gc_right=cand.gc_right,
                    pair_penalty=cand.pair_penalty,
                    upl_probe_id=str(getattr(hit, "probe_id")),
                    upl_probe_seq=str(getattr(hit, "probe_seq")),
                    upl_probe_strand=str(getattr(hit, "strand", "+")),
                    probe_start_in_amplicon=int(hit.start),
                    probe_end_in_amplicon=int(hit.end),
                    dist_left_3p_to_probe=dist_l,
                    dist_probe_to_right_3p=dist_r,
                    left_start=left_start,
                    left_end=left_end,
                    right_start=right_start,
                    right_end=right_end,
                    left_len=cand.left_len,
                    right_len=cand.right_len,
                    amplicon_start=product_start,
                    amplicon_end=product_end - 1,
                    junction_spanning=junction_spanning,
                    junction_spanning_detail=detail,
                    tx_left_exact_hits=None,
                    tx_right_exact_hits=None,
                    genome_left_exact_hits=None,
                    genome_right_exact_hits=None,
                    tx_product_offtargets=None,
                    genome_product_offtargets=None,
                    tx_amplicons=[],
                    genome_amplicons=[],
                    score=score,
                )
            )
            if progress_cb is not None and total_cands > 0:
                if i == total_cands - 1 or i % 10 == 0:
                    pct = 15.0 + (45.0 * (i + 1) / total_cands)
                    progress_cb(pct, "Matching UPL probes...")

    ranked = _dedup_primer_pairs(ranked)
    ranked = _apply_exon_constraint_with_fallback(ranked, warnings=warnings)
    ranked = _maybe_add_specificity(
        ranked,
        inputs=inputs,
        ensembl=ensembl,
        warnings=warnings,
        target_gene_symbol=target_gene_symbol,
        target_transcript_ids=target_transcript_ids,
        progress_cb=progress_cb,
    )
    ranked.sort(key=lambda x: (-x.score, x.pair_penalty is None, x.pair_penalty or 1e9))
    if progress_cb is not None:
        progress_cb(100.0, "Done.")
    spec_params = None
    if inputs.specificity_mode != "none":
        spec_params = {
            "mode": inputs.specificity_mode,
            "max_offtarget_amplicons": inputs.max_offtarget_amplicons,
            "min_mismatches_total": inputs.min_mismatches_total,
            "min_mismatches_3p": inputs.min_mismatches_3p,
            "three_prime_window": inputs.three_prime_window,
            "ignore_mismatches_total_ge": inputs.ignore_mismatches_total_ge,
            "max_target_amplicon_size": inputs.max_target_amplicon_size,
            "top_n": inputs.specificity_top_n,
            "blastn_parallel_jobs": inputs.blastn_parallel_jobs,
        }
    return DesignResult(
        species=inputs.species,
        transcript_info=transcript_info,
        pairs=ranked,
        warnings=warnings,
        specificity_params=spec_params,
    )


def run_design_with_specificity_retries(
    *,
    inputs: DesignInputs,
    ensembl: EnsemblClient,
    upl_probes: dict[str, str],
    progress_cb: Callable[[float, str], None] | None = None,
) -> tuple[DesignResult, list[DesignResult]]:
    attempts: list[DesignResult] = []
    current_inputs = inputs
    if current_inputs.specificity_mode != "none":
        ensure_default_blastdbs(Path.cwd(), progress_cb=progress_cb)
    while True:
        if progress_cb is not None and attempts:
            progress_cb(
                55.0,
                f"Retrying specificity with max_offtarget_amplicons={current_inputs.max_offtarget_amplicons}...",
            )
        start_dt = datetime.now().isoformat(timespec="seconds")
        t0 = perf_counter()
        result = run_design_workflow(
            inputs=current_inputs,
            ensembl=ensembl,
            upl_probes=upl_probes,
            progress_cb=progress_cb,
        )
        result.design_started_at = start_dt
        result.design_elapsed_sec = perf_counter() - t0
        attempts.append(result)
        if current_inputs.specificity_mode == "none":
            break
        confirmed = [p for p in result.pairs if getattr(p, "specificity_passed", None) is True]
        if confirmed:
            break
        if not result.pairs:
            break
        if all(getattr(p, "specificity_passed", None) is None for p in result.pairs):
            result.warnings.append(
                "Specificity を再評価できません（BLAST未実行/失敗の可能性）。/ "
                "Specificity re-check skipped due to BLAST issues."
            )
            break
        result.warnings.append(
            "Specificity confirmed が0件のため max_offtarget_amplicons を増やして再評価します。/ "
            "Retrying with higher max_offtarget_amplicons."
        )
        current_inputs = replace(current_inputs, max_offtarget_amplicons=current_inputs.max_offtarget_amplicons + 2)
    return (attempts[-1], attempts)


def _normalize_seq(seq: str) -> str:
    s = seq.upper().replace("U", "T")
    s = "".join([c for c in s if c in {"A", "C", "G", "T", "N"}])
    return s


def _transcript_info_dict(t: TranscriptDetails) -> dict[str, Any]:
    return {
        "transcript_id": t.transcript_id,
        "transcript_display_name": t.transcript_display_name,
        "gene_id": t.gene_id,
        "gene_symbol": t.gene_symbol,
        "strand": t.strand,
        "exon_count": len(t.exons),
        "cdna_length": len(t.cdna_sequence),
    }


def _cdna_exon_boundaries(t: TranscriptDetails) -> list[int]:
    if not t.exons:
        return []
    strand = int(t.strand or 1)
    exons = sorted(t.exons, key=lambda e: int(e.get("start", 0)), reverse=(strand == -1))
    boundaries: list[int] = []
    acc = 0
    for i, e in enumerate(exons):
        exon_len = int(e.get("end")) - int(e.get("start")) + 1
        acc += exon_len
        if i < len(exons) - 1:
            boundaries.append(acc - 1)  # 0-based index of last base of exon i
    return boundaries


def _ranked_pair_from_rust_row(row: dict[str, Any]) -> RankedPrimerPair:
    tm_left_raw = row.get("tm_left")
    tm_right_raw = row.get("tm_right")
    gc_left_raw = row.get("gc_left")
    gc_right_raw = row.get("gc_right")
    pair_penalty_raw = row.get("pair_penalty")
    return RankedPrimerPair(
        left_seq=str(row["left_seq"]),
        right_seq=str(row["right_seq"]),
        product_size=int(row["product_size"]),
        tm_left=None if tm_left_raw is None else float(tm_left_raw),
        tm_right=None if tm_right_raw is None else float(tm_right_raw),
        gc_left=None if gc_left_raw is None else float(gc_left_raw),
        gc_right=None if gc_right_raw is None else float(gc_right_raw),
        pair_penalty=None if pair_penalty_raw is None else float(pair_penalty_raw),
        upl_probe_id=str(row["upl_probe_id"]),
        upl_probe_seq=str(row["upl_probe_seq"]),
        upl_probe_strand=str(row.get("upl_probe_strand", "+")),
        probe_start_in_amplicon=int(row["probe_start_in_amplicon"]),
        probe_end_in_amplicon=int(row["probe_end_in_amplicon"]),
        dist_left_3p_to_probe=int(row["dist_left_3p_to_probe"]),
        dist_probe_to_right_3p=int(row["dist_probe_to_right_3p"]),
        left_start=int(row["left_start"]),
        left_end=int(row["left_end"]),
        right_start=int(row["right_start"]),
        right_end=int(row["right_end"]),
        left_len=int(row["left_len"]),
        right_len=int(row["right_len"]),
        amplicon_start=int(row["amplicon_start"]),
        amplicon_end=int(row["amplicon_end"]),
        junction_spanning=bool(row.get("junction_spanning", False)),
        junction_spanning_detail=str(row.get("junction_spanning_detail", "no_junction")),
        tx_left_exact_hits=None,
        tx_right_exact_hits=None,
        genome_left_exact_hits=None,
        genome_right_exact_hits=None,
        tx_product_offtargets=None,
        genome_product_offtargets=None,
        tx_amplicons=[],
        genome_amplicons=[],
        score=float(row["score"]),
    )


def _junction_flags(
    *, boundaries: list[int], product_start: int, left_len: int, right_len: int, product_size: int
) -> tuple[bool, str]:
    if not boundaries:
        return (False, "no_exon_info")

    # Coordinates in transcript cDNA (0-based, inclusive)
    left_0 = product_start
    left_1 = product_start + left_len - 1
    right_3p = product_start + product_size - 1

    spans_left = any((left_0 <= b < left_1) for b in boundaries)
    spans_amplicon = any((product_start <= b < product_start + product_size - 1) for b in boundaries)

    if spans_left:
        return (True, "left_primer_spans_junction")
    if spans_amplicon:
        # amplicon spans junction even if left primer doesn't (right primer may be in another exon)
        return (True, "amplicon_spans_junction")
    # right primer junction-spanning check
    right_start_guess = right_3p - (right_len - 1)
    spans_right = any((right_start_guess <= b < right_3p) for b in boundaries)
    if spans_right:
        return (True, "right_primer_spans_junction (approx)")
    return (False, "no_junction")


def _probe_preference(*, product_size: int, probe_start: int, probe_end: int, dist_l: int, dist_r: int) -> float:
    """
    同一primer pairで複数probe候補がある場合の選好スコア。
    - 3'端から十分離れていることは上流で担保（min_probe_offset_bp）
    - できるだけ中央寄り & 左右の距離がバランスしているものを優先
    """
    min_dist = min(dist_l, dist_r)
    imbalance = abs(dist_l - dist_r)
    probe_mid = (probe_start + probe_end) / 2.0
    center = (product_size - 1) / 2.0
    center_offset = abs(probe_mid - center)
    return (min_dist * 0.2) - (imbalance * 0.5) - (center_offset * 0.2)


def _score_candidate(
    cand: PrimerPairCandidate, *, junction_spanning: bool, dist_l: int, dist_r: int, probe_mid: float, product_size: int
) -> float:
    score = 0.0
    if cand.pair_penalty is not None:
        score += max(0.0, 10.0 - cand.pair_penalty)
    if junction_spanning:
        score += 5.0
    # Prefer probe not too close to primer 3' ends (handled by min_probe_offset_bp),
    # and prefer roughly central / balanced placement.
    min_dist = min(dist_l, dist_r)
    imbalance = abs(dist_l - dist_r)
    center = (product_size - 1) / 2.0
    center_offset = abs(probe_mid - center)
    score += min(30.0, float(min_dist)) * 0.05
    score -= min(60.0, float(imbalance)) * 0.02
    score -= min(60.0, float(center_offset)) * 0.01
    return score


def _dedup_primer_pairs(pairs: list[RankedPrimerPair]) -> list[RankedPrimerPair]:
    if not pairs:
        return pairs
    best: dict[tuple[str, str, int, int, int], RankedPrimerPair] = {}
    for p in pairs:
        key = (p.left_seq, p.right_seq, int(p.left_start), int(p.right_start), int(p.product_size))
        prev = best.get(key)
        if prev is None:
            best[key] = p
            continue
        if p.score > prev.score:
            best[key] = p
            continue
        if p.score == prev.score:
            if min(p.dist_left_3p_to_probe, p.dist_probe_to_right_3p) > min(
                prev.dist_left_3p_to_probe, prev.dist_probe_to_right_3p
            ):
                best[key] = p
                continue
    out = list(best.values())
    out.sort(key=lambda x: (-x.score, x.pair_penalty is None, x.pair_penalty or 1e9))
    return out


def _apply_exon_constraint_with_fallback(
    pairs: list[RankedPrimerPair], *, warnings: list[str]
) -> list[RankedPrimerPair]:
    if not pairs:
        return pairs
    if all(p.junction_spanning_detail == "no_exon_info" for p in pairs):
        warnings.append(
            "Exon情報が無いため、junction跨ぎ条件を検証できませんでした。/ "
            "No exon info; junction constraint not enforced."
        )
        return pairs
    primer_span = [
        p
        for p in pairs
        if p.junction_spanning_detail
        in {
            "left_primer_spans_junction",
            "right_primer_spans_junction (approx)",
        }
    ]
    if primer_span:
        return primer_span
    amplicon_span = [p for p in pairs if p.junction_spanning_detail == "amplicon_spans_junction"]
    if amplicon_span:
        warnings.append(
            "exon-exon境界を跨ぐプライマーが無いため、amplicon内junctionを許容しました。/ "
            "No exon-spanning primers; allowing amplicon-spanning junctions."
        )
        return amplicon_span
    warnings.append(
        "exon-exon境界/イントロン跨ぎ条件を満たす候補が無いため、注釈付きで制約を緩和しました。/ "
        "No exon/intron-spanning candidates; constraint relaxed with annotation."
    )
    return pairs


def _has_3p_gc_run(seq: str, run_len: int = 3) -> bool:
    """
    Check if primer has consecutive G/C bases at 3' end.
    Having 3+ consecutive G/C at 3' end can cause:
    - Strong secondary structures
    - Non-specific binding
    This is stricter than Primer3's PRIMER_GC_CLAMP which just requires
    a minimum number of G/C in the last few bases.
    """
    count = 0
    for base in reversed(seq.upper()):
        if base in {"G", "C"}:
            count += 1
            if count >= run_len:
                return True
        else:
            break
    return False


def _has_poly_run(seq: str, run_len: int = 4) -> bool:
    """
    Check for runs of identical bases (poly-A, poly-T, poly-G, poly-C).
    Poly-runs can cause:
    - Mispriming
    - Secondary structure formation
    Note: Primer3's PRIMER_MAX_POLY_X=3 should already filter these,
    but this provides an additional safety check.
    """
    if not seq:
        return False
    last = seq[0]
    count = 1
    for ch in seq[1:]:
        if ch == last:
            count += 1
            if count >= run_len:
                return True
        else:
            last = ch
            count = 1
    return False


def _resolve_parallel_jobs(n: int) -> int:
    if n <= 0:
        return max(1, (os.cpu_count() or 1) - 1)
    return max(1, n)


def _list_amplicon_products(
    *,
    left_hits: list[BlastHit],
    right_hits: list[BlastHit],
    min_size: int,
    max_size: int,
) -> list[dict[str, Any]]:
    left_by_contig: dict[str, list[BlastHit]] = {}
    right_by_contig: dict[str, list[BlastHit]] = {}
    for h in left_hits:
        left_by_contig.setdefault(h.sseqid, []).append(h)
    for h in right_hits:
        right_by_contig.setdefault(h.sseqid, []).append(h)

    products: dict[tuple[str, int, int], dict[str, Any]] = {}
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
                    key = (contig, l3, r3)
                    if key not in products:
                        products[key] = {
                            "contig": contig,
                            "left_3p": l3,
                            "right_3p": r3,
                            "size": size,
                            "left_mismatch_total": l.mismatch_total,
                            "right_mismatch_total": r.mismatch_total,
                            "left_mismatch_3p": l.mismatch_3p,
                            "right_mismatch_3p": r.mismatch_3p,
                        }
    out = list(products.values())
    out.sort(key=lambda x: (x["contig"], x["left_3p"], x["right_3p"]))
    return out


def _maybe_add_specificity(
    pairs: list[RankedPrimerPair],
    *,
    inputs: DesignInputs,
    ensembl: EnsemblClient | None,
    warnings: list[str],
    target_gene_symbol: str | None = None,
    target_transcript_ids: set[str] | None = None,
    progress_cb: Callable[[float, str], None] | None = None,
) -> list[RankedPrimerPair]:
    if inputs.specificity_mode == "none":
        return pairs
    if inputs.transcriptome_blast_db is None or inputs.genome_blast_db is None:
        warnings.append(
            "Specificity: local_blast を選択しましたが、Transcriptome/Genome のDBパスが未設定です。"
            " / Transcriptome/Genome DB path is missing."
        )
        return pairs
    if not blastn_available(inputs.blastn_path):
        warnings.append(f"Specificity: blastn が見つかりません: {inputs.blastn_path} / blastn not found")
        return pairs

    target_tx_ids = target_transcript_ids if target_transcript_ids else None
    gene_cache: dict[str, str | None] = {}

    def _gene_symbol_for_transcript(tid: str) -> str | None:
        if target_tx_ids is not None and target_gene_symbol and tid in target_tx_ids:
            return target_gene_symbol
        if tid in gene_cache:
            return gene_cache[tid]
        if ensembl is None:
            return None
        if not tid.startswith("ENS"):
            gene_cache[tid] = None
            return None
        try:
            sym = ensembl.get_transcript_gene_symbol(tid)
        except Exception:
            sym = None
        gene_cache[tid] = sym
        return sym

    def _annotate_tx_amplicons(items: list[dict[str, Any]]) -> None:
        for a in items:
            tid = a.get("contig")
            if not isinstance(tid, str) or not tid:
                continue
            if target_tx_ids is not None:
                if target_gene_symbol and tid in target_tx_ids:
                    a["gene_symbol"] = target_gene_symbol
                else:
                    a["gene_symbol"] = None
                continue
            a["gene_symbol"] = _gene_symbol_for_transcript(tid)

    def _effective_tx_product_count(items: list[dict[str, Any]]) -> int:
        """
        Transcriptome側のアンプリコン数の数え方:
        - target_gene_symbol が分かる場合、同一遺伝子由来の複数variant(転写産物)は「1つの目的産物」として扱う
        - それ以外（別遺伝子 or gene不明）は非特異的産物として個別にカウントする
        """
        if target_tx_ids is None and not target_gene_symbol:
            return len(items)
        target_seen = False
        non_target = 0
        for a in items:
            tid = a.get("contig")
            gene = a.get("gene_symbol")
            is_target = False
            if target_tx_ids is not None and isinstance(tid, str) and tid:
                is_target = tid in target_tx_ids
            elif target_gene_symbol:
                is_target = bool(gene) and (gene == target_gene_symbol)
            a["is_target_gene"] = is_target
            if is_target:
                target_seen = True
            else:
                non_target += 1
        return (1 if target_seen else 0) + non_target

    def _evaluate_pair(p: RankedPrimerPair) -> RankedPrimerPair:
        if inputs.specificity_mode == "local_blast (in_silico_pcr)":
            tx_l_hits_raw = blastn_hits_detailed(
                primer_seq=p.left_seq, db=inputs.transcriptome_blast_db, blastn_path=inputs.blastn_path
            )
            tx_r_hits_raw = blastn_hits_detailed(
                primer_seq=p.right_seq, db=inputs.transcriptome_blast_db, blastn_path=inputs.blastn_path
            )
            g_l_hits_raw = blastn_hits_detailed(
                primer_seq=p.left_seq, db=inputs.genome_blast_db, blastn_path=inputs.blastn_path
            )
            g_r_hits_raw = blastn_hits_detailed(
                primer_seq=p.right_seq, db=inputs.genome_blast_db, blastn_path=inputs.blastn_path
            )

            tx_l_hits = filter_hits_in_silico(
                hits=tx_l_hits_raw,
                primer_len=len(p.left_seq),
                min_mismatches_total=inputs.min_mismatches_total,
                min_mismatches_3p=inputs.min_mismatches_3p,
                three_prime_window=inputs.three_prime_window,
                ignore_mismatches_total_ge=inputs.ignore_mismatches_total_ge,
            )
            tx_r_hits = filter_hits_in_silico(
                hits=tx_r_hits_raw,
                primer_len=len(p.right_seq),
                min_mismatches_total=inputs.min_mismatches_total,
                min_mismatches_3p=inputs.min_mismatches_3p,
                three_prime_window=inputs.three_prime_window,
                ignore_mismatches_total_ge=inputs.ignore_mismatches_total_ge,
            )
            g_l_hits = filter_hits_in_silico(
                hits=g_l_hits_raw,
                primer_len=len(p.left_seq),
                min_mismatches_total=inputs.min_mismatches_total,
                min_mismatches_3p=inputs.min_mismatches_3p,
                three_prime_window=inputs.three_prime_window,
                ignore_mismatches_total_ge=inputs.ignore_mismatches_total_ge,
            )
            g_r_hits = filter_hits_in_silico(
                hits=g_r_hits_raw,
                primer_len=len(p.right_seq),
                min_mismatches_total=inputs.min_mismatches_total,
                min_mismatches_3p=inputs.min_mismatches_3p,
                three_prime_window=inputs.three_prime_window,
                ignore_mismatches_total_ge=inputs.ignore_mismatches_total_ge,
            )
        elif inputs.specificity_mode == "local_blast (exact match count)":
            tx_l_hits = _exact_hits(
                blastn_hits(primer_seq=p.left_seq, db=inputs.transcriptome_blast_db, blastn_path=inputs.blastn_path),
                primer_len=len(p.left_seq),
            )
            tx_r_hits = _exact_hits(
                blastn_hits(primer_seq=p.right_seq, db=inputs.transcriptome_blast_db, blastn_path=inputs.blastn_path),
                primer_len=len(p.right_seq),
            )
            g_l_hits = _exact_hits(
                blastn_hits(primer_seq=p.left_seq, db=inputs.genome_blast_db, blastn_path=inputs.blastn_path),
                primer_len=len(p.left_seq),
            )
            g_r_hits = _exact_hits(
                blastn_hits(primer_seq=p.right_seq, db=inputs.genome_blast_db, blastn_path=inputs.blastn_path),
                primer_len=len(p.right_seq),
            )
        else:
            return p

        tx_amplicons = _list_amplicon_products(
            left_hits=tx_l_hits,
            right_hits=tx_r_hits,
            min_size=inputs.product_size_min,
            max_size=inputs.max_target_amplicon_size,
        )
        _annotate_tx_amplicons(tx_amplicons)
        tx_products = _effective_tx_product_count(tx_amplicons)
        g_amplicons = _list_amplicon_products(
            left_hits=g_l_hits,
            right_hits=g_r_hits,
            min_size=inputs.product_size_min,
            max_size=inputs.max_target_amplicon_size,
        )
        _annotate_tx_amplicons(g_amplicons)
        g_products = len(g_amplicons)

        penalty = 0.0
        for c in (len(tx_l_hits), len(tx_r_hits), len(g_l_hits), len(g_r_hits)):
            if c <= 1:
                continue
            penalty += min(5.0, (c - 1) * 1.0)
        if tx_products > inputs.max_offtarget_amplicons:
            penalty += 10.0 + (tx_products - inputs.max_offtarget_amplicons) * 2.0
        if g_products > inputs.max_offtarget_amplicons:
            penalty += 10.0 + (g_products - inputs.max_offtarget_amplicons) * 2.0

        specificity_passed = (tx_products <= inputs.max_offtarget_amplicons) and (
            g_products <= inputs.max_offtarget_amplicons
        )

        return RankedPrimerPair(
            **_with_specificity(
                p,
                len(tx_l_hits),
                len(tx_r_hits),
                len(g_l_hits),
                len(g_r_hits),
                tx_products,
                g_products,
                tx_amplicons,
                g_amplicons,
                p.score - penalty,
                specificity_passed,
            )
        )

    out = list(pairs)
    top_n = max(1, int(inputs.specificity_top_n or 20))
    target_pairs = pairs[:top_n]
    jobs = _resolve_parallel_jobs(inputs.blastn_parallel_jobs)
    if jobs <= 1 or len(target_pairs) <= 1:
        completed = 0
        for i, p in enumerate(target_pairs):
            try:
                out[i] = _evaluate_pair(p)
            except Exception as e:
                warnings.append(f"Specificity: BLAST 実行でエラー（以降スキップ）: {e} / BLAST error")
                break
            completed += 1
            if progress_cb is not None:
                pct = 60.0 + (35.0 * completed / max(1, len(target_pairs)))
                progress_cb(pct, "Checking specificity with BLAST...")
        return out

    futures: list[tuple[int, object]] = []
    completed = 0
    with ThreadPoolExecutor(max_workers=jobs) as executor:
        for idx, p in enumerate(target_pairs):
            futures.append((idx, executor.submit(_evaluate_pair, p)))
        for idx, fut in futures:
            try:
                res = fut.result()
            except Exception as e:
                warnings.append(f"Specificity: BLAST 実行でエラー（以降スキップ）: {e} / BLAST error")
                break
            out[idx] = res
            completed += 1
            if progress_cb is not None:
                pct = 60.0 + (35.0 * completed / max(1, len(target_pairs)))
                progress_cb(pct, "Checking specificity with BLAST...")
    return out


def _with_specificity(
    p: RankedPrimerPair,
    tx_left: int,
    tx_right: int,
    genome_left: int,
    genome_right: int,
    tx_products: int,
    genome_products: int,
    tx_amplicons: list[dict[str, Any]],
    genome_amplicons: list[dict[str, Any]],
    new_score: float,
    specificity_passed: bool,
) -> dict[str, Any]:
    d: dict[str, Any] = dict(p.__dict__)
    d["tx_left_exact_hits"] = tx_left
    d["tx_right_exact_hits"] = tx_right
    d["genome_left_exact_hits"] = genome_left
    d["genome_right_exact_hits"] = genome_right
    d["tx_product_offtargets"] = tx_products
    d["genome_product_offtargets"] = genome_products
    d["tx_amplicons"] = tx_amplicons
    d["genome_amplicons"] = genome_amplicons
    d["score"] = new_score
    d["specificity_passed"] = specificity_passed
    return d


def _exact_hits(hits: list, *, primer_len: int) -> list:
    return [h for h in hits if getattr(h, "mismatch", None) == 0 and getattr(h, "length", None) == primer_len and getattr(h, "pident", None) == 100.0]
