from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
from typing import Any

import tomllib
import sys
import os
import importlib.resources as importlib_resources

from upl_finder_fast.ensembl import EnsemblClient, Species
from upl_finder_fast.report import design_results_markdown
from upl_finder_fast.upl import load_upl_probes, rust_extension_available, set_rust_mode
from upl_finder_fast.workflow import DesignInputs, run_design_with_specificity_retries


def _load_logo_text() -> str:
    # Priority:
    # 1) explicit env var
    # 2) local package file (upl_finder_fast/UPL_Finder.txt)
    # 3) packaged resources fallback
    candidates: list[Path] = []
    env = os.getenv("UPL_FINDER_LOGO_PATH")
    if env:
        candidates.append(Path(env))
    candidates.append(Path(__file__).with_name("UPL_Finder.txt"))

    for p in candidates:
        try:
            if p.is_file():
                return p.read_text(encoding="utf-8").rstrip("\n")
        except Exception:
            pass

    try:
        return importlib_resources.files("upl_finder_fast").joinpath("UPL_Finder.txt").read_text(encoding="utf-8").rstrip(
            "\n"
        )
    except Exception:
        return "upl_finder_fast"


def _read_fasta(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    seq: list[str] = []
    for line in text.splitlines():
        if not line or line.startswith(">"):
            if seq:
                break
            continue
        seq.append(line.strip())
    return "".join(seq)


def _species_from_str(s: str) -> Species:
    v = s.strip().lower()
    if v in {"human", "homo_sapiens", "hs"}:
        return Species.HUMAN
    if v in {"mouse", "mus_musculus", "mm"}:
        return Species.MOUSE
    raise ValueError(f"Unsupported species: {s}")


def _safe_target_label(raw: str) -> str:
    out = "".join(c if c.isalnum() or c in {"-", "_"} else "_" for c in str(raw)).strip("_")
    return out or "target"


def _load_params(path: Path) -> dict[str, Any]:
    return tomllib.loads(path.read_text(encoding="utf-8"))


def _resolve_path(base: Path, value: str | None) -> Path | None:
    if value is None:
        return None
    p = Path(value)
    if p.is_absolute():
        return p
    return (base / p).resolve()


def _build_inputs(cfg: dict[str, Any], base_dir: Path, input_type: str, raw_input: str) -> DesignInputs:
    design = cfg.get("design", {})
    spec = cfg.get("specificity", {})
    species = _species_from_str(str(design.get("species", "human")))
    tx_db = _resolve_path(base_dir, spec.get("transcriptome_blast_db") or None)
    g_db = _resolve_path(base_dir, spec.get("genome_blast_db") or None)
    return DesignInputs(
        species=species,
        input_type=input_type,
        raw_input=raw_input,
        primer_tm_opt=float(design.get("primer_tm_opt", 59.5)),
        primer_tm_min=float(design.get("primer_tm_min", 59.0)),
        primer_tm_max=float(design.get("primer_tm_max", 60.0)),
        primer_tm_diff_max=float(design.get("primer_tm_diff_max", 3.0)),
        product_size_min=int(design.get("product_size_min", 60)),
        product_size_max=int(design.get("product_size_max", 150)),
        min_probe_offset_bp=int(design.get("min_probe_offset_bp", 10)),
        max_pairs=int(design.get("max_pairs", 50)),
        selected_transcript_id=design.get("selected_transcript_id") or None,
        specificity_mode=str(spec.get("mode", "local_blast (in_silico_pcr)")),
        blastn_path=str(spec.get("blastn_path", "blastn")),
        transcriptome_blast_db=str(tx_db) if tx_db else None,
        genome_blast_db=str(g_db) if g_db else None,
        specificity_top_n=int(spec.get("top_n", 20)),
        max_offtarget_amplicons=int(spec.get("max_offtarget_amplicons", 1)),
        min_mismatches_total=int(spec.get("min_mismatches_total", 2)),
        min_mismatches_3p=int(spec.get("min_mismatches_3p", 2)),
        three_prime_window=int(spec.get("three_prime_window", 5)),
        ignore_mismatches_total_ge=int(spec.get("ignore_mismatches_total_ge", 6)),
        max_target_amplicon_size=int(spec.get("max_target_amplicon_size", 1000)),
        blastn_parallel_jobs=int(spec.get("blastn_parallel_jobs", -1)),
    )


def _resolve_rust_mode(cfg: dict[str, Any], cli_mode: str | None) -> str:
    if cli_mode is not None:
        return cli_mode
    runtime = cfg.get("runtime", {})
    mode = str(runtime.get("rust", "auto")).strip().lower()
    if mode not in {"auto", "on", "off"}:
        raise ValueError(f"Unsupported runtime.rust mode: {mode}")
    return mode


def main() -> None:
    print(_load_logo_text())
    p = argparse.ArgumentParser(description="Primer design CLI (UPL + specificity)")
    p.add_argument("--param", default="parameter.toml", help="Path to parameter.toml")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--gene", help="Gene symbol")
    g.add_argument("--transcript", help="Ensembl transcript ID")
    g.add_argument("--seq", help="Target cDNA sequence")
    g.add_argument("--fasta", help="FASTA file path")
    p.add_argument("--rust", choices=["auto", "on", "off"], help="Rust acceleration mode")
    p.add_argument("--output", help="Output markdown path (optional)")
    args = p.parse_args()

    param_path = Path(args.param)
    cfg = _load_params(param_path)
    rust_mode = _resolve_rust_mode(cfg, args.rust)
    set_rust_mode(rust_mode)
    if rust_mode == "on" and not rust_extension_available():
        raise RuntimeError(
            "Rust mode 'on' was requested but upl_finder_rust is not available. "
            "Build/install it from ./upl_finder_rust with maturin."
        )
    paths = cfg.get("paths", {})
    upl_probe_value = (
        paths.get("upl_probe_tsv")
        or paths.get("upl_probe_path")
        or paths.get("upl_probe_pkl")  # backward compat
        or "roche_upl_sequences.tsv"
    )
    upl_probe_path = _resolve_path(param_path.parent, str(upl_probe_value)) or Path("roche_upl_sequences.tsv")

    if args.gene:
        input_type = "Gene symbol"
        raw_input = args.gene.strip()
        target_label = raw_input
    elif args.transcript:
        input_type = "Ensembl Transcript ID"
        raw_input = args.transcript.strip()
        target_label = raw_input
    elif args.seq:
        input_type = "Paste cDNA sequence"
        raw_input = args.seq.strip()
        target_label = "sequence"
    else:
        fasta_path = Path(args.fasta)
        input_type = "Paste cDNA sequence"
        raw_input = _read_fasta(fasta_path)
        target_label = fasta_path.stem

    inputs = _build_inputs(cfg, param_path.parent, input_type=input_type, raw_input=raw_input)
    ensembl = EnsemblClient()
    upl_probes = load_upl_probes(str(upl_probe_path))

    last_pct: float | None = None
    last_msg: str | None = None

    def _progress(pct: float, msg: str) -> None:
        nonlocal last_pct, last_msg
        pct_i = int(max(0, min(100, pct)))
        if last_pct == pct_i and last_msg == msg:
            return
        last_pct = pct_i
        last_msg = msg
        sys.stderr.write(f"[{pct_i:3d}%] {msg}\n")
        sys.stderr.flush()

    _progress(0, "Designing primers and matching UPL probes...")
    result, attempts = run_design_with_specificity_retries(
        inputs=inputs,
        ensembl=ensembl,
        upl_probes=upl_probes,
        progress_cb=_progress,
    )

    yymmddhhmm = datetime.now().strftime("%y%m%d_%H%M")
    md = design_results_markdown(attempts)
    outputs: list[Path] = []
    if args.output:
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(md, encoding="utf-8")
        outputs.append(out_path)
    else:
        out_dir = _resolve_path(param_path.parent, cfg.get("output", {}).get("out_dir", "upl_primer_probe")) or Path(
            "upl_primer_probe"
        )
        out_dir.mkdir(parents=True, exist_ok=True)
        safe_target = _safe_target_label(target_label)
        out_path = out_dir / f"primer_upl_{safe_target}_{yymmddhhmm}.md"
        out_path.write_text(md, encoding="utf-8")
        outputs.append(out_path)
    for p in outputs:
        print(f"Wrote: {p}")


if __name__ == "__main__":
    main()
