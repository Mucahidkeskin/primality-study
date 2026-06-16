# Primality Tests — code and data

Source code, datasets, and result figures for the experiments in the article
*Primality Tests: Methods and Experimental Results*. The article presents six classical
primality tests — trial division, Fermat, Solovay–Strassen, Miller–Rabin, AKS, and
Lucas–Lehmer — and reports reproducible operation-count measurements comparing their cost and
reliability.

This repository is the **code-and-data artifact** for reproducing those measurements; the paper
itself is maintained separately.

- **Tool:** [`primality.c`](primality.c) — a self-contained C implementation of all six tests
  with exact operation-count instrumentation and no external libraries.

## Repository layout
| Path | Contents |
|------|----------|
| `primality.c` | The C tool (interactive menu + non-interactive `--csv` batch mode). |
| `make_datasets.sh` | Regenerates the measurement datasets in `data/`. |
| `plot.py` | Regenerates the result figures from `data/`. |
| `figures/` | Result figures (PNG, SVG). |
| `data/` | Measured datasets (CSV) and `data/README.md`. |
| `SHA256SUMS.txt` | SHA-256 hashes of the source, scripts, and datasets. |

## Build the tool
```sh
gcc primality.c -O2 -o primality.exe
```
No external libraries are required.

## Run
Interactive (enter a number, pick a method, or "run all"):
```sh
./primality.exe
```
Batch / CSV mode (one CSV row per method to stdout):
```sh
./primality.exe --csv --method=all --k=20 --seed=88172645463325252
# reads "<number>[,<label>]" lines from stdin; methods: all|trial|fermat|solovay|miller|aks|lucas
```

## Reproduce the measurements
```sh
gcc primality.c -O2 -o primality.exe   # 1. build
bash make_datasets.sh                  # 2. regenerate data/*.csv
python plot.py                         # 3. regenerate figures/
```
Requirements: a C compiler; OpenSSL (probable-prime generation); Python 3 + matplotlib (figures).
The probabilistic experiments use the fixed seed `88172645463325252`, `k = 20` rounds, a
trial-division bound of `10^6`, and an AKS cap of 20 bits.

## Verify integrity
```sh
sha256sum -c SHA256SUMS.txt
```

## Cost metric
Cost is reported as deterministic operation counts (word multiplications, modular
multiplications/squarings, and method-specific counts) under a fixed 32-bit-limb schoolbook
arithmetic implementation, rather than wall-clock time. See the article for the exact
counter-increment rules.
