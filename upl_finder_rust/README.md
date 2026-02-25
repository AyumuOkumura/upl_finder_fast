# upl_finder_rust

`upl_finder_fast_repo` 同梱の Rust 拡張です。`upl_finder_fast` CLI から任意で利用できます（Linux想定）。

## Build (from repository root)

```bash
conda activate upl_finder_fast
python -m pip install maturin
python -m maturin develop --release --manifest-path upl_finder_rust/Cargo.toml
```

`upl_finder_fast --rust on ...` で Rust 使用を強制、`--rust off` で Python 実装に固定できます。

## Python usage

```python
import upl_finder_rust

matches = upl_finder_rust.find_upl_matches_from_file(
    "ACTGACTGACTG",
    "../roche_upl_sequences.tsv",  # also supports .json/.txt
)
for m in matches:
    print(m.probe_id, m.probe_seq, m.start, m.end, m.strand)
```

## Rust CLI binary

```bash
cargo run --release -- ../roche_upl_sequences.tsv ACTGACTGACTG  # also supports .json/.txt
```
