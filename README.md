```text
  _   _ ____  _       _____ _           _           
 | | | |  _ \| |     |  ___(_)_ __   __| | ___ _ __ 
 | | | | |_) | |     | |_  | | '_ \ / _` |/ _ \ '__|
 | |_| |  __/| |___  |  _| | | | | | (_| |  __/ |   
  \___/|_|   |_____| |_|   |_|_| |_|\__,_|\___|_|   
                                                    
```

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

## Specification (Current) / 現行仕様（要点）

この節は **実装の現在の挙動** を短くまとめたものです（READMEの他の説明より“仕様”寄り）。  
This section summarizes the **current implemented behavior** (more “spec” than marketing).

### 1) Inputs / 入力
- `--gene <SYMBOL>`: species に応じて Ensembl API から転写産物候補を取得し、（必要なら）`selected_transcript_id` を優先して 1 つ選択します。  
  Fetches transcript candidates from Ensembl API and selects one (prefers `selected_transcript_id` when provided).
- `--transcript <ENST...>`: 指定transcriptのcDNAを Ensembl API から取得します。  
  Fetches the transcript cDNA from Ensembl API.
- `--refseq <NM_...>`: RefSeq transcript ID から Ensembl transcript を解決し、cDNAを Ensembl API から取得します（複数候補がある場合は先頭を採用し警告に候補一覧を出します）。  
  Resolves RefSeq transcript ID to an Ensembl transcript, then fetches cDNA from Ensembl API (if multiple candidates are found, uses the first and prints a warning with candidates).
- `--seq <cDNA>` / `--fasta <file>`: cDNA配列を直接使用します（Ensembl問い合わせ無し）。  
  Uses the provided cDNA sequence directly (no Ensembl lookup).

入力配列は A/C/G/T/N のみ残し、U→Tに正規化します。  
Input sequences are normalized (U→T; keep only A/C/G/T/N).

### 2) Primer design / プライマー設計
Primer3（`primer3-py`）で primer pair 候補を生成します。  
Primer3 (`primer3-py`) generates candidate primer pairs.

主な制約（`parameter.toml` の `[design]`）:
- 産物長（`product_size_min`〜`product_size_max`）
- Tm（`primer_tm_*` と `primer_tm_diff_max`）
- 追加の簡易フィルタ: 3'端のGC連続、poly-run などは除外  
  Additional filters: 3' GC runs, poly-runs, etc.

### 3) UPL probe matching / UPLプローブ適合
各 primer pair ごとに、amplicon（産物配列）を切り出し、Roche UPL probe（8-mer）の完全一致を探索します。  
For each primer pair, the amplicon is extracted and Roche UPL probes (8-mers) are searched by **exact match**.

- 探索は **probe配列そのもの（`+`）と逆相補（`-`）の両方** を対象にします。  
  Both **probe sequence (`+`) and its reverse-complement (`-`)** are searched.
- Rust版（`--rust auto|on` で有効）の内部実装は **Aho–Corasick法（マルチパターン文字列探索）** で、プローブ配列と逆相補配列の集合からオートマトンを一度構築し、各ampliconに対して **重なりを許す一致（overlapping match）** を列挙します。候補（primer pair）評価は `rayon` による並列処理です。  
  The Rust path (enabled by `--rust auto|on`) uses **Aho–Corasick multi-pattern search**: it builds one automaton from probe patterns (+ reverse-complements) and enumerates **overlapping** exact matches on each amplicon, evaluating candidates in parallel via `rayon`.
- 進捗表示 `[ 60%] Matching UPL probes...` はこの工程のラベルですが、Rust有効時はマッチング計算がRust側で先に完了してから表示されるため、長時間の計算中に同じ行で止まって見える場合があります（Pythonフォールバック時は候補ループ中に更新されます）。  
  The progress label `[ 60%] Matching UPL probes...` corresponds to this phase, but with Rust enabled the heavy matching work can complete inside Rust before the label is printed (Python fallback updates the label during the candidate loop).
- `min_probe_offset_bp`（現在の推奨デフォルト: **10 bp**）により、各プライマーの3'端とプローブの間の距離（間に挟まるbp数。隣接=0）を下回るものは除外します（プライマーとプローブの**重複は常に除外**）。  
  `min_probe_offset_bp` (recommended default **10 bp**) filters probes that are too close to primer 3' ends (gap in bp between them; adjacent=0). Primer/probe **overlap is always rejected**.
- 同一の primer pair に複数の probe がヒットする場合、**中央寄り・左右バランスが良い** probe を 1 つ選びます（候補爆発を抑制）。  
  When multiple probes match for the same primer pair, the tool keeps **one “more ideal”** probe (more central/balanced) to avoid candidate explosion.
- 結果には `upl_probe_strand`（`+`/`-`）を出力します。  
  The output includes `upl_probe_strand` (`+`/`-`).

### 4) Exon junction preference / exon-junction優先
Ensembl由来のexon情報がある場合、exon-exon junction を跨ぐ設計（primerまたはampliconがjunctionを跨ぐ）を優先します。  
If Ensembl exon info is available, designs spanning exon-exon junctions are preferred.

### 5) Specificity (optional) / 特異性チェック（任意）
`[specificity].mode != "none"` の場合、BLAST+（`blastn`）で転写産物DB・ゲノムDBに対するオフターゲットを評価します（上位 `top_n` 件のみ）。  
When `[specificity].mode != "none"`, BLAST+ (`blastn`) evaluates off-targets against transcriptome/genome DBs (only top `top_n` pairs).

- `local_blast (in_silico_pcr)`: 3'近傍ミスマッチも考慮して危険側ヒットを残し、産物サイズ範囲で成立するアンプリコン数を数えます。  
  Filters “risky” hits using mismatch rules and counts in-range amplicons (in silico PCR-like).
- `local_blast (exact match count)`: 完全一致ヒット数のみで簡易評価します。  
  Uses only exact-match counts for a simpler check.
- DBが無い場合は Ensembl FASTA をダウンロードして `data/blastdb/` に構築できます（`makeblastdb` が必要）。  
  Missing DBs can be built under `data/blastdb/` from Ensembl FASTA (requires `makeblastdb`).

### 6) Ranking and output / ランキングと出力
- `score` は primer3 penalty、junction優先、probe位置（中央/バランス）、特異性ペナルティ等を合わせて計算します。  
  `score` combines primer3 penalty, junction preference, probe placement, and optional specificity penalties.
- 上位候補を Markdown レポートとして保存します（`--output` 省略時は `output.out_dir/`）。  
  Writes a Markdown report (defaults to `output.out_dir/` when `--output` is omitted).

## インストール / Install
- 前提: Python 3.10+, BLAST+（`blastn` と `makeblastdb`）。  
  Requirements: Python 3.10+, BLAST+ (`blastn`, `makeblastdb`).
- **仮想環境名は `upl_finder_fast` に統一**しています。  
  **Use a single environment name: `upl_finder_fast`**.

### Step-by-step (GitHub公開後の利用想定)

```bash
git clone <YOUR_GITHUB_REPO_URL>
cd upl_finder_fast_repo
conda env create -f environment.yml
conda activate upl_finder_fast
```

注意: `environment.yml` の `pip:` に `-e .`（このリポジトリを editable install）を含むため、**必ず `upl_finder_fast_repo/` に `cd` した状態で** `conda env create -f environment.yml` を実行してください。  
If you run `conda env create` from another directory, `-e .` may point to the wrong path and you may end up running an older CLI.

（任意）Rust拡張を使う場合:

```bash
python -m pip install maturin
python -m maturin develop --release --manifest-path upl_finder_rust/Cargo.toml
```

`maturin` を実行しない場合でも CLI は動作し、Python 実装に自動フォールバックします。  
Without `maturin`, CLI still works and falls back to Python implementation.

## 実行 / Run

```bash
upl_finder_fast --param parameter.toml --gene GAPDH
```

Run (English): same as above.

Other inputs / その他の入力例:

```bash
upl_finder_fast --param parameter.toml --transcript ENST00000323125
upl_finder_fast --param parameter.toml --refseq NM_000546
upl_finder_fast --param parameter.toml --fasta target.fa
upl_finder_fast --param parameter.toml --seq ACTG...   # paste cDNA
```

- `parameter.toml` の編集: `[design]` で標的種（例: GAPDH なら `human`）やTm・産物サイズを調整し、`[specificity]` でBLAST DBパスやモード（特異性不要なら `mode = "none"`）を設定します。`[paths]` の `upl_probe_tsv`（または互換キー `upl_probe_path`）で同梱の `roche_upl_sequences.tsv` か、自作の probe 定義ファイル（`.tsv`/`.txt`/`.json`）を指定します。  
  Edit `parameter.toml`: tune species/Tm/product size under `[design]`, set BLAST DB paths or `mode = "none"` under `[specificity]`, and point `[paths].upl_probe_tsv` (or compat key `upl_probe_path`) to the bundled `roche_upl_sequences.tsv` or your own probe definition file (`.tsv`/`.txt`/`.json`).
- `--param`: `parameter.toml` で設計条件（Tm範囲、プローブ距離など）と特異性用BLAST DBのパス (`specificity.transcriptome_blast_db`, `specificity.genome_blast_db`) を設定します。  
  Configure design constraints and BLAST DB paths in `parameter.toml` (`design.*`, `specificity.*`, `paths.upl_probe_tsv`).
- `--output` を省略すると `output.out_dir`（デフォルト `upl_primer_probe/`）に `primer_upl_<target>_<YYMMDD_HHMM>.md` を保存します。  
  When `--output` is omitted, a timestamped Markdown report is written under `output.out_dir` (default `upl_primer_probe/`).
- Roche UPL プローブ配列は `roche_upl_sequences.tsv`（同梱）から読み込みます（`.txt`/`.json` 形式も利用可能）。  
  Roche UPL probe sequences are read from the bundled `roche_upl_sequences.tsv` (also supports `.txt`/`.json`).

## CLI-first Rust toggle / CLIでのRust切替

このリポジトリは **CLI実行を主用途** にしています。Rust利用有無はCLIで選択できます。  
This repository is organized for **CLI-first usage**. Rust on/off can be selected from CLI.

```bash
# auto (default): 利用可能ならRust、無ければPython
upl_finder_fast --param parameter.toml --rust auto --gene GAPDH

# Rustを必須化（拡張未インストール時はエラー）
upl_finder_fast --param parameter.toml --rust on --gene GAPDH

# Rustを使わずPython実装のみ
upl_finder_fast --param parameter.toml --rust off --gene GAPDH
```

- `parameter.toml` の `[runtime].rust = "auto|on|off"` でも設定可能です（CLI引数が優先）。
- 互換のため `UPL_FINDER_DISABLE_RUST=1` でもRustを無効化できます。
- Rust拡張のビルド/導入手順は `upl_finder_rust/README.md` を参照してください。  
  For building/installing the Rust extension, see `upl_finder_rust/README.md`.

## Repository layout / リポジトリ構成

- `upl_finder_fast/` : CLI本体（Python）
- `upl_finder_rust/` : Rust拡張（PyO3/maturin）
- `parameter.toml` : CLI設定（`[design]`, `[specificity]`, `[paths]`, `[runtime]`）
- `roche_upl_sequences.tsv` : UPLプローブ辞書（同梱。`.txt`/`.json` 形式も利用可能）
- `UPL_Finder.txt` : CLI/README で使用するASCIIロゴ

## Acknowledgements / 謝辞

- 本プロジェクトのASCIIアート（ロゴ）は `npx oh-my-logo` を用いて生成しました（project: `https://github.com/shinshin86/oh-my-logo`）。
