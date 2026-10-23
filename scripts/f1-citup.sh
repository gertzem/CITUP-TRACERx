#! /bin/bash

set -eu

for ground_truth in data/sims/*; do
    case=$(basename $ground_truth)
    citup_out=results/min_trees/sims/$case/trees
    if [ ! -d "$citup_out" ]; then
        continue
    fi
    python scripts/calc_f1_scores.py --type citupm "$ground_truth/G.tree" $(fd 'mike_edges_sol0(_refined_\d+)?\.tsv' "$citup_out")
done