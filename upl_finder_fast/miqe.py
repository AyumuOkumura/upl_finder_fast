from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass
class MIQEReport:
    species: str | None = None
    target_gene: str | None = None
    transcript_id: str | None = None

    primer_left: str | None = None
    primer_right: str | None = None
    upl_probe_id: str | None = None
    upl_probe_seq: str | None = None
    amplicon_length_bp: int | None = None

    primer_tm_left: float | None = None
    primer_tm_right: float | None = None
    primer_gc_left: float | None = None
    primer_gc_right: float | None = None

    exon_junction_spanning: bool | None = None

    # 実験側（フォーム入力を想定 / MVPでは空のまま）
    rt_priming_method: str | None = None
    qpcr_mastermix: str | None = None
    qpcr_instrument: str | None = None
    efficiency_percent: float | None = None
    efficiency_method: str | None = None
    reference_gene: str | None = None
    normalization_method: str | None = None

    @classmethod
    def from_design_result(cls, result: Any) -> "MIQEReport":
        best = result.pairs[0] if result.pairs else None
        return cls(
            species=getattr(result, "species", None).value if getattr(result, "species", None) is not None else None,
            target_gene=(result.transcript_info or {}).get("gene_symbol") if result.transcript_info else None,
            transcript_id=(result.transcript_info or {}).get("transcript_id") if result.transcript_info else None,
            primer_left=(best.left_seq if best else None),
            primer_right=(best.right_seq if best else None),
            upl_probe_id=(best.upl_probe_id if best else None),
            upl_probe_seq=(best.upl_probe_seq if best else None),
            amplicon_length_bp=(best.product_size if best else None),
            primer_tm_left=(best.tm_left if best else None),
            primer_tm_right=(best.tm_right if best else None),
            primer_gc_left=(best.gc_left if best else None),
            primer_gc_right=(best.gc_right if best else None),
            exon_junction_spanning=(best.junction_spanning if best else None),
        )


def miqe_missing_fields(report: MIQEReport) -> list[str]:
    missing: list[str] = []
    required = [
        ("species", report.species),
        ("target_gene", report.target_gene),
        ("transcript_id", report.transcript_id),
        ("primer_left", report.primer_left),
        ("primer_right", report.primer_right),
        ("upl_probe_id", report.upl_probe_id),
        ("amplicon_length_bp", report.amplicon_length_bp),
    ]
    for k, v in required:
        if v is None or (isinstance(v, str) and not v.strip()):
            missing.append(k)
    return missing


def miqe_experimental_fields_missing(report: MIQEReport) -> list[str]:
    """
    アプリ管轄外（実験で得る/記録することが多い）項目の未入力リスト。
    ここは「警告」ではなく、任意のガイドとして扱う想定。
    """
    missing: list[str] = []
    required = [
        ("qpcr_mastermix", report.qpcr_mastermix),
        ("qpcr_instrument", report.qpcr_instrument),
        ("efficiency_percent", report.efficiency_percent),
        ("efficiency_method", report.efficiency_method),
        ("normalization_method", report.normalization_method),
    ]
    for k, v in required:
        if v is None or (isinstance(v, str) and not v.strip()):
            missing.append(k)
    return missing
