#!/usr/bin/env bash
# Generate the article datasets (CSV) for the primality-test cost study.
# Requires: ./primality.exe, openssl, python (for Mersenne decimals).
set -e
cd "$(dirname "$0")"
EXE=./primality.exe
HDR="label,bits,method,verdict,word_muls,word_div_steps,mod_muls,mod_sqrs,mod_reductions,high_label,high_value,k,seed"
mkdir -p data

echo "[Fig 1] scaling on probable primes (16..2048 bits, 3 reps each)"
IN=data/_fig1_in.txt; : > "$IN"
for bits in 16 32 64 128 256 512 1024 2048; do
  for rep in 1 2 3; do
    p=$(openssl prime -generate -bits "$bits" 2>/dev/null)
    printf '%s,b%04d_r%d\n' "$p" "$bits" "$rep" >> "$IN"
  done
done
$EXE --csv --method=all --aks-cap=0 --k=20 < "$IN" > data/fig1_scaling.csv

echo "[Fig 2] reliability: Carmichael numbers + strong pseudoprimes (k=20)"
cat > data/_fig2_in.txt <<'EOF'
561,carmichael
1105,carmichael
1729,carmichael
2465,carmichael
2821,carmichael
6601,carmichael
8911,carmichael
41041,carmichael
294409,carmichael_bigfactors
56052361,carmichael_bigfactors
118901521,carmichael_bigfactors
2047,strong_psp2
3277,strong_psp2
4033,strong_psp2
4681,strong_psp2
8321,strong_psp2
3215031751,strong_psp_2357
EOF
$EXE --csv --method=all --aks-cap=0 --k=20 < data/_fig2_in.txt > data/fig2_reliability.csv

echo "[Fig 2b] Fermat fooled vs Miller across 10 seeds (large-factor Carmichaels, k=5)"
OUT=data/fig2b_fermat_seeds.csv; echo "$HDR" > "$OUT"
for seed in 1 2 3 4 5 6 7 8 9 10; do
  for n in 294409 56052361 118901521; do
    for meth in fermat miller; do
      printf '%s,bigCarmichael\n' "$n" \
        | $EXE --csv --method=$meth --k=5 --seed=$seed | tail -n +2 >> "$OUT"
    done
  done
done

echo "[Fig 3] equal-k vs equal-confidence on one 1024-bit prime (k=5..80)"
P1024=$(openssl prime -generate -bits 1024 2>/dev/null)
OUT=data/fig3_confidence.csv; echo "$HDR" > "$OUT"
for k in 5 10 20 40 80; do
  printf '%s,p1024\n' "$P1024" \
    | $EXE --csv --method=all --aks-cap=0 --k=$k | tail -n +2 >> "$OUT"
done

echo "[Fig 4] AKS blowup on small primes (<=17 bits)"
cat > data/_fig4_in.txt <<'EOF'
31,p05bit
127,p07bit
1009,p10bit
8191,p13bit
32749,p15bit
131071,p17bit
EOF
$EXE --csv --method=all --aks-cap=20 --k=20 < data/_fig4_in.txt > data/fig4_aks.csv

echo "[Fig 5] Lucas-Lehmer niche on Mersenne numbers (vs Miller-Rabin)"
IN=data/_fig5_in.txt; : > "$IN"
for p in 31 61 89 107 127 521 607 1279; do
  m=$(python -c "print(2**$p-1)")
  printf '%s,Mp%04d_primeExp\n' "$m" "$p" >> "$IN"
done
for p in 67 257 1061; do
  m=$(python -c "print(2**$p-1)")
  printf '%s,Mp%04d_composite\n' "$m" "$p" >> "$IN"
done
$EXE --csv --method=all --aks-cap=0 --k=10 --trial-bound=1000 < "$IN" > data/fig5_lucas.csv

rm -f data/_fig*_in.txt
echo "done -> data/*.csv"
