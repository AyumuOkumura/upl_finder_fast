from __future__ import annotations

from dataclasses import dataclass

try:
    import primer3  # type: ignore
except Exception:  # pragma: no cover
    primer3 = None


@dataclass(frozen=True)
class PrimerPairCandidate:
    left_seq: str
    right_seq: str
    product_size: int
    left_len: int
    right_len: int
    tm_left: float | None
    tm_right: float | None
    gc_left: float | None
    gc_right: float | None
    pair_penalty: float | None
    left_start: int
    right_start: int

    def product_range(self) -> tuple[int, int]:
        return (self.left_start, self.left_start + self.product_size)


def design_primers_primer3(
    template: str,
    *,
    product_size_min: int,
    product_size_max: int,
    primer_tm_opt: float,
    primer_tm_min: float,
    primer_tm_max: float,
    primer_tm_diff_max: float,
    max_pairs: int,
) -> list[PrimerPairCandidate]:
    if primer3 is None:
        raise RuntimeError(
            "primer3-py が見つかりません。`pip install -r requirements.txt` を実行して依存関係を導入してください。"
        )
    seq_args = {"SEQUENCE_TEMPLATE": template}
    global_args = {
        "PRIMER_NUM_RETURN": int(max_pairs),
        "PRIMER_PRODUCT_SIZE_RANGE": [[int(product_size_min), int(product_size_max)]],
        "PRIMER_OPT_SIZE": 20,
        "PRIMER_MIN_SIZE": 18,
        "PRIMER_MAX_SIZE": 24,
        "PRIMER_OPT_TM": float(primer_tm_opt),
        "PRIMER_MIN_TM": float(primer_tm_min),
        "PRIMER_MAX_TM": float(primer_tm_max),
        "PRIMER_PAIR_MAX_DIFF_TM": float(primer_tm_diff_max),
        "PRIMER_MIN_GC": 30.0,
        "PRIMER_MAX_GC": 80.0,
        "PRIMER_GC_CLAMP": 2,
        "PRIMER_MAX_POLY_X": 3,
        "PRIMER_MAX_SELF_ANY": 8.0,
        "PRIMER_MAX_SELF_END": 3.0,
        "PRIMER_PAIR_MAX_COMPL_ANY": 8.0,
        "PRIMER_PAIR_MAX_COMPL_END": 3.0,
        "PRIMER_MAX_HAIRPIN_TH": 24.0,
        "PRIMER_SALT_MONOVALENT": 50.0,
        "PRIMER_DNA_CONC": 50.0,
        "PRIMER_MAX_NS_ACCEPTED": 0,
    }

    # primer3-py deprecated designPrimers; use design_primers when available
    try:
        res = primer3.bindings.design_primers(seq_args, global_args)
    except AttributeError:  # pragma: no cover
        res = primer3.bindings.designPrimers(seq_args, global_args)
    n = int(res.get("PRIMER_PAIR_NUM_RETURNED", 0) or 0)
    out: list[PrimerPairCandidate] = []
    for i in range(min(n, max_pairs)):
        left = res.get(f"PRIMER_LEFT_{i}_SEQUENCE")
        right = res.get(f"PRIMER_RIGHT_{i}_SEQUENCE")
        if not left or not right:
            continue
        product_size = int(res.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE"))
        left_pos = res.get(f"PRIMER_LEFT_{i}")
        right_pos = res.get(f"PRIMER_RIGHT_{i}")
        left_start = int(left_pos[0]) if isinstance(left_pos, (list, tuple)) and len(left_pos) >= 2 else 0
        left_len = int(left_pos[1]) if isinstance(left_pos, (list, tuple)) and len(left_pos) >= 2 else len(left)
        right_len = int(right_pos[1]) if isinstance(right_pos, (list, tuple)) and len(right_pos) >= 2 else len(right)
        right_start = int(right_pos[0]) if isinstance(right_pos, (list, tuple)) and len(right_pos) >= 2 else 0

        out.append(
            PrimerPairCandidate(
                left_seq=left,
                right_seq=right,
                product_size=product_size,
                left_len=left_len,
                right_len=right_len,
                tm_left=_to_float(res.get(f"PRIMER_LEFT_{i}_TM")),
                tm_right=_to_float(res.get(f"PRIMER_RIGHT_{i}_TM")),
                gc_left=_to_float(res.get(f"PRIMER_LEFT_{i}_GC_PERCENT")),
                gc_right=_to_float(res.get(f"PRIMER_RIGHT_{i}_GC_PERCENT")),
                pair_penalty=_to_float(res.get(f"PRIMER_PAIR_{i}_PENALTY")),
                left_start=left_start,
                right_start=right_start,
            )
        )
    return out


def _to_float(x: object) -> float | None:
    try:
        return float(x)  # type: ignore[arg-type]
    except Exception:
        return None
