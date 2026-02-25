from __future__ import annotations

from decimal import Decimal, ROUND_HALF_UP
from typing import Any


def _probe_tm_simple(seq: str) -> float:
    s = seq.upper()
    a = s.count("A")
    t = s.count("T")
    g = s.count("G")
    c = s.count("C")
    return float(2 * (a + t) + 4 * (g + c))


def _round_1dp_half_up(v: float) -> str:
    return str(Decimal(str(v)).quantize(Decimal("0.0"), rounding=ROUND_HALF_UP))


def _amplicon_diagram_lines(p: Any, *, width: int = 60) -> list[str]:
    length = int(getattr(p, "product_size", 0) or 0)
    if length <= 0:
        return ["(no amplicon)"]
    w = max(10, min(width, length))
    scale = (length - 1) / (w - 1) if w > 1 else 1.0

    def _to_idx(pos: int) -> int:
        if length <= 1:
            return 0
        idx = int(round(pos / scale))
        return max(0, min(w - 1, idx))

    primers = ["-"] * w
    probe = ["-"] * w
    a_start = int(getattr(p, "amplicon_start", 0) or 0)
    left_start = int(getattr(p, "left_start", 0) or 0) - a_start
    left_len = int(getattr(p, "left_len", 0) or 0)
    right_start = int(getattr(p, "right_start", 0) or 0) - a_start
    right_len = int(getattr(p, "right_len", 0) or 0)
    probe_start = int(getattr(p, "probe_start_in_amplicon", 0) or 0)
    probe_end = int(getattr(p, "probe_end_in_amplicon", 0) or 0)

    def _mark(arr: list[str], start: int, length_bp: int, ch: str) -> None:
        if length_bp <= 0:
            return
        s = max(0, min(length - 1, start))
        e = max(0, min(length - 1, start + length_bp - 1))
        if e < s:
            return
        i0 = _to_idx(s)
        i1 = _to_idx(e)
        for i in range(i0, i1 + 1):
            arr[i] = ch if arr[i] == "-" else "X"

    _mark(primers, left_start, left_len, "L")
    _mark(primers, right_start, right_len, "R")
    _mark(probe, probe_start, max(1, probe_end - probe_start + 1), "P")

    legend_primers = "Primers: L=left  R=right  X=overlap(scaled)"
    legend_probe = "Probe: P=probe"
    coord = f"0{' ' * (w - 2)}{length - 1}" if w >= 2 else "0"

    # Exact (non-scaled) overlap in bp, for sanity-checking.
    left_end = left_start + left_len - 1
    right_end = right_start + right_len - 1
    ov_left = max(0, min(probe_end, left_end) - max(probe_start, left_start) + 1)
    ov_right = max(0, min(probe_end, right_end) - max(probe_start, right_start) + 1)
    return [
        f"amplicon_len={length}bp",
        coord,
        "".join(primers),
        legend_primers,
        "".join(probe),
        legend_probe,
        f"overlap_bp(left,probe)={ov_left} overlap_bp(probe,right)={ov_right}",
        f"left=[{left_start},{left_start + left_len - 1}] right=[{right_start},{right_start + right_len - 1}] "
        f"probe=[{probe_start},{probe_end}]",
    ]


def _field_descriptions_lines() -> list[str]:
    return [
        "**項目の意味 / Field Descriptions**",
        "- `species`: 設計対象の生物種 / Species used for design",
        "- `design_started_at`: 設計開始日時 / Design start datetime",
        "- `design_elapsed_min`: 設計時間（分） / Total design time (minutes)",
        "- `transcript_info`: 選択された転写産物の情報 / Selected transcript metadata",
        "- `pairs`: プライマーペア候補一覧 / List of primer pair candidates",
        "- `pairs[].left_seq`: 左プライマー配列 / Left primer sequence",
        "- `pairs[].right_seq`: 右プライマー配列 / Right primer sequence",
        "- `pairs[].tm_left`: 左プライマーTm（℃） / Left primer Tm (°C)",
        "- `pairs[].tm_right`: 右プライマーTm（℃） / Right primer Tm (°C)",
        "- `pairs[].gc_left`: 左プライマーGC% / Left primer GC%",
        "- `pairs[].gc_right`: 右プライマーGC% / Right primer GC%",
        "- `pairs[].left_start/left_end`: 左プライマー位置 / Left primer position",
        "- `pairs[].right_start/right_end`: 右プライマー位置 / Right primer position",
        "- `pairs[].left_len/right_len`: プライマー長 / Primer length",
        "- `pairs[].amplicon_start/amplicon_end`: 産物位置 / Amplicon position",
        "- `pairs[].tx_*`/`genome_*`: 特異性情報 / Specificity info",
        "- `pairs[].upl_probe_id/seq/tm`: UPL probe情報 / UPL probe info",
        "- `pairs[].upl_probe_strand`: probeの向き（amplicon内の一致方向） / Probe match strand within amplicon",
        "- `specificity_params.min_mismatches_total`: ミスマッチ下限(全体) / Min mismatches (total)",
        "- `specificity_params.min_mismatches_3p`: 3'窓ミスマッチ下限 / Min mismatches (3' window)",
        "- `specificity_params.ignore_mismatches_total_ge`: 無視するミスマッチ下限 / Ignore mismatches >= N",
        "- `specificity_params.max_target_amplicon_size`: 最大アンプリコン長 / Max target amplicon size",
        "",
    ]


def design_result_markdown(
    result: Any, *, include_title: bool = True, include_field_descriptions: bool = True
) -> str:
    transcript_info = result.transcript_info or {}
    species_obj = getattr(result, "species", None)
    species_value = getattr(species_obj, "value", species_obj)
    lines: list[str] = []
    if include_title:
        lines.append("# Primer UPL Design Result / 結果")
        lines.append("")
    if include_field_descriptions:
        lines.extend(_field_descriptions_lines())
    lines.append("**結果 / Result**")
    lines.append("")
    lines.append(f"Species: `{species_value}`")
    started_at = getattr(result, "design_started_at", None)
    elapsed_sec = getattr(result, "design_elapsed_sec", None)
    if started_at:
        lines.append(f"Design started at: `{started_at}`")
    if elapsed_sec is not None:
        elapsed_min = float(elapsed_sec) / 60.0
        lines.append(f"Design elapsed (min): `{elapsed_min:.2f}`")
    lines.append("")
    lines.append("**Transcript Info**")
    lines.append("| Key | Value |")
    lines.append("| --- | --- |")
    for key in [
        "transcript_id",
        "transcript_display_name",
        "gene_id",
        "gene_symbol",
        "strand",
        "exon_count",
        "cdna_length",
        "gene_symbol_input",
        "selected_transcript_id",
        "selected_transcript_biotype",
        "note",
    ]:
        value = transcript_info.get(key)
        lines.append(f"| `{key}` | `{value}` |")
    lines.append("")
    lines.append("**Warnings**")
    lines.append("| Warning |")
    lines.append("| --- |")
    warnings = getattr(result, "warnings", []) or []
    if warnings:
        for w in warnings:
            lines.append(f"| {w} |")
    else:
        lines.append("| None |")
    lines.append("")
    spec_params = getattr(result, "specificity_params", None)
    if spec_params:
        lines.append("**Specificity Params**")
        lines.append("| Key | Value |")
        lines.append("| --- | --- |")
        for k, v in spec_params.items():
            lines.append(f"| `{k}` | `{v}` |")
        lines.append("")
    all_pairs = getattr(result, "pairs", []) or []
    specific_pairs = [p for p in all_pairs if getattr(p, "specificity_passed", None) is True]
    lines.append(f"**Pairs Count (All)**: `{len(all_pairs)}`")
    lines.append(f"**Pairs Count (Specificity confirmed)**: `{len(specific_pairs)}`")
    lines.append("")
    pairs = specific_pairs[:5]
    if not pairs:
        lines.append("No primer pairs with confirmed specificity were found.")
        lines.append("")
    for i, p in enumerate(pairs, start=1):
        lines.append(f"**Pair {i}**")
        lines.append("| Key | Value |")
        lines.append("| --- | --- |")
        upl_tm = _probe_tm_simple(getattr(p, "upl_probe_seq", "") or "")
        for key in [
            "left_seq",
            "right_seq",
            "product_size",
            "tm_left",
            "tm_right",
            "gc_left",
            "gc_right",
            "left_start",
            "left_end",
            "right_start",
            "right_end",
            "left_len",
            "right_len",
            "amplicon_start",
            "amplicon_end",
            "pair_penalty",
            "upl_probe_id",
            "upl_probe_seq",
            "upl_probe_strand",
            "probe_start_in_amplicon",
            "probe_end_in_amplicon",
            "dist_left_3p_to_probe",
            "dist_probe_to_right_3p",
            "junction_spanning",
            "junction_spanning_detail",
            "tx_left_exact_hits",
            "tx_right_exact_hits",
            "genome_left_exact_hits",
            "genome_right_exact_hits",
            "tx_product_offtargets",
            "genome_product_offtargets",
            "specificity_passed",
            "score",
        ]:
            value = getattr(p, key)
            if key in {"tm_left", "tm_right"} and isinstance(value, (int, float)):
                value = _round_1dp_half_up(float(value))
            lines.append(f"| `{key}` | `{value}` |")
        lines.append(f"| `upl_probe_tm` | `{upl_tm:.1f}` |")
        lines.append("")
        lines.append("**Amplicon diagram (relative) / 位置の簡易図**")
        lines.append("```text")
        lines.extend(_amplicon_diagram_lines(p))
        lines.append("```")
        lines.append("")
        lines.append("Predicted Amplicons (Transcriptome)")
        tx_amplicons = getattr(p, "tx_amplicons", []) or []
        if tx_amplicons:
            for a in tx_amplicons:
                mm_t = a.get("left_mismatch_total")
                mm_r = a.get("right_mismatch_total")
                mm3_t = a.get("left_mismatch_3p")
                mm3_r = a.get("right_mismatch_3p")
                gene = a.get("gene_symbol")
                extra = ""
                if mm_t is not None or mm_r is not None:
                    extra = f", mismatches L/R={mm_t}/{mm_r}, 3p={mm3_t}/{mm3_r}"
                if gene:
                    extra = f"{extra}, gene={gene}"
                is_target = a.get("is_target_gene")
                if is_target is True:
                    extra = f"{extra}, target_gene=1"
                elif is_target is False:
                    extra = f"{extra}, target_gene=0"
                contig = a.get("contig")
                if gene and contig:
                    contig = f"{contig} ({gene})"
                lines.append(
                    f"- {contig}:{a.get('left_3p')}-{a.get('right_3p')} (size={a.get('size')}{extra})"
                )
        else:
            lines.append("- None")
        lines.append("")
        lines.append("Predicted Amplicons (Genome)")
        g_amplicons = getattr(p, "genome_amplicons", []) or []
        if g_amplicons:
            for a in g_amplicons:
                mm_t = a.get("left_mismatch_total")
                mm_r = a.get("right_mismatch_total")
                mm3_t = a.get("left_mismatch_3p")
                mm3_r = a.get("right_mismatch_3p")
                gene = a.get("gene_symbol")
                extra = ""
                if mm_t is not None or mm_r is not None:
                    extra = f", mismatches L/R={mm_t}/{mm_r}, 3p={mm3_t}/{mm3_r}"
                gene_label = gene if gene else "NA"
                extra = f"{extra}, gene={gene_label}"
                contig = a.get("contig")
                contig_disp = f"{contig} ({gene_label})" if contig else contig
                lines.append(
                    f"- {contig_disp}:{a.get('left_3p')}-{a.get('right_3p')} (size={a.get('size')}{extra})"
                )
        else:
            lines.append("- None")
        lines.append("")
    return "\n".join(lines)


def design_results_markdown(results: list[Any]) -> str:
    if not results:
        return ""
    if len(results) == 1:
        return design_result_markdown(results[0])
    lines: list[str] = []
    lines.append("# Primer UPL Design Results (All Attempts) / 全試行結果")
    lines.append("")
    lines.extend(_field_descriptions_lines())
    for i, res in enumerate(results, start=1):
        spec_params = getattr(res, "specificity_params", None) or {}
        max_off = spec_params.get("max_offtarget_amplicons")
        lines.append(f"## Attempt {i} (max_offtarget_amplicons={max_off})")
        lines.append("")
        lines.append(design_result_markdown(res, include_field_descriptions=False))
        lines.append("")
    return "\n".join(lines)
