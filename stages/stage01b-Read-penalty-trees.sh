#! /bin/bash

set -euf


python scripts/read_citup_obj.py \
    $(find results/penalty_trees -name '*.tree_0.txt') > results/penalty_trees/tracerx/optima.tsv

python scripts/plot_by_tree_path.py \
    $(find results/penalty_trees -name trees)
