# upl_finder_fast

高速版のCLIパッケージです（Rust拡張が利用可能な場合は自動で使用します）。Fast CLI for qPCR primer + Roche UPL probe design with optional Rust acceleration.

## Overview / 概要
- ターゲット決定: 遺伝子記号・Ensembl Transcript ID・cDNA配列（貼り付け/FASTA）を受け取り、必要なら Ensembl API から cDNA を取得して正規化します。  
  Target resolution: accepts gene symbols, Ensembl transcript IDs, pasted cDNA strings, or FASTA files and fetches/normalizes cDNA via the Ensembl API when needed.
- プライマー設計とUPL対応: primer3 で設計した左右プライマー候補を Roche UPL プローブ全体と突き合わせ、プローブの配置距離やTm/GC条件を満たす組み合わせを探索します（Rust拡張が存在すれば自動利用）。  
  Primer + UPL pairing: designs primer pairs with primer3, matches them against all Roche UPL probes, and keeps pairs that satisfy distance/Tm/GC constraints, using the Rust extension automatically when available.
- 特異性チェック（任意）: BLAST+ (`blastn`/`makeblastdb`) を使って転写産物・ゲノムDBへのオフターゲット（in silico PCR）をカウントし、許容閾値を超える組を除外します。DB が無ければ Ensembl FASTA を自動ダウンロードして `data/blastdb/` に構築します。  
  Specificity (optional): uses BLAST+ against transcriptome/genome databases to score off-target amplicons and drops pairs beyond the thresholds; missing DBs are built automatically from Ensembl FASTA files under `data/blastdb/`.
- 出力: スコア順に並べたマークダウンを生成し、プライマー/プローブ配列、座標、Tm/GC、接合部跨ぎ情報、オフターゲット件数などを記録します。  
  Output: produces a ranked Markdown report containing primer/probe sequences, positions, Tm/GC, junction-spanning flags, and off-target counts.

## インストール / Install
- 前提: Python 3.10+, BLAST+（`blastn` と `makeblastdb`）が必要です。`environment.yml` には BLAST と Python 3.12、依存パッケージが含まれています。  
  Requirements: Python 3.10+, BLAST+ (`blastn`, `makeblastdb`). The provided `environment.yml` includes BLAST, Python 3.12, and all Python deps.

## Conda 環境 / Conda Environment

```bash
conda env create -f environment.yml
conda activate upl_finder_fast
```

Conda environment setup (English): same as above.

Pip (editable) install if you already have dependencies:  
Editable install via pip when deps are available:

```bash
python -m pip install -e .
```

## 実行 / Run

```bash
upl_finder_fast --param parameter.toml --gene GAPDH
```

Run (English): same as above.

Other inputs / その他の入力例:

```bash
upl_finder_fast --param parameter.toml --transcript ENST00000323125
upl_finder_fast --param parameter.toml --fasta target.fa
upl_finder_fast --param parameter.toml --seq ACTG...   # paste cDNA
```

- `parameter.toml` の編集: `[design]` で標的種（例: GAPDH なら `human`）やTm・産物サイズを調整し、`[specificity]` でBLAST DBパスやモード（特異性不要なら `mode = "none"`）を設定します。`[paths]` の `upl_probe_pkl` は同梱の `roche_upl_sequences.pkl` か、自作の {probe_id: sequence} pickle に差し替えられます。  
  Edit `parameter.toml`: tune species/Tm/product size under `[design]`, set BLAST DB paths or `mode = "none"` under `[specificity]`, and point `[paths].upl_probe_pkl` to the bundled `roche_upl_sequences.pkl` or your own `{probe_id: sequence}` pickle.
- `--param`: `parameter.toml` で設計条件（Tm範囲、プローブ距離など）と特異性用BLAST DBのパス (`specificity.transcriptome_blast_db`, `specificity.genome_blast_db`) を設定します。  
  Configure design constraints and BLAST DB paths in `parameter.toml` (`design.*`, `specificity.*`, `paths.upl_probe_pkl`).
- `--output` を省略すると `output.out_dir`（デフォルト `upl_primer_probe/`）に `primer_upl_<target>_<YYMMDD_HHMM>.md` を保存します。  
  When `--output` is omitted, a timestamped Markdown report is written under `output.out_dir` (default `upl_primer_probe/`).
- Roche UPL プローブ配列は pickle 版の `roche_upl_sequences.pkl`（同梱）から読み込みます。  
  Roche UPL probe sequences are read from the bundled pickle `roche_upl_sequences.pkl`.
