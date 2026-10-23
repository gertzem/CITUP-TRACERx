#! /bin/bash

set -eu

for ground_truth in data/sims/*; do
    case=$(basename $ground_truth)
    mach2_out=data/conipher_mach2_runs/threshold_0p025/$case/tree_1/mach2_out
    if [ ! -d "$mach2_out" ]; then
        continue
    fi
    python scripts/calc_f1_scores.py --type mach2 "$ground_truth/G.tree" "$mach2_out"/P-G-*.graph
done