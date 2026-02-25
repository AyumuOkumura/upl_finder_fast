```text
  _   _ ____  _       _____ _           _           
 | | | |  _ \| |     |  ___(_)_ __   __| | ___ _ __ 
 | | | | |_) | |     | |_  | | '_ \ / _` |/ _ \ '__|
 | |_| |  __/| |___  |  _| | | | | | (_| |  __/ |   
  \___/|_|   |_____| |_|   |_|_| |_|\__,_|\___|_|   
                                                    
```

# upl_finder_fast

高速版のCLIパッケージです（Rust拡張が利用可能な場合は自動で使用します）。

## Conda 環境 / Conda Environment

```bash
cd upl_finder_fast_repo
conda env create -f environment.yml
conda activate upl_finder_fast
```

Conda environment setup (English):

```bash
cd upl_finder_fast_repo
conda env create -f environment.yml
conda activate upl_finder_fast
```

## 実行 / Run

```bash
upl_finder_fast --param parameter.toml --gene GAPDH
```

Rust mode examples:

```bash
upl_finder_fast --param parameter.toml --rust auto --gene GAPDH
upl_finder_fast --param parameter.toml --rust on --gene GAPDH
upl_finder_fast --param parameter.toml --rust off --gene GAPDH
```

Run (English):

```bash
upl_finder_fast --param parameter.toml --gene GAPDH
```

## Acknowledgements / 謝辞

- 本プロジェクトのASCIIアート（ロゴ）は `npx oh-my-logo` を用いて生成しました（project: `https://github.com/shinshin86/oh-my-logo`）。
